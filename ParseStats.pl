#!/usr/bin/perl

# Perl script written by Joseph Hughes, University of Glasgow
# use this to count the overall errors, in the forward and the reverse and create the matrix for circos 

use strict;
use Getopt::Long; 

my ($stats,$out,$help);
&GetOptions(
	    'stats:s'      => \$stats,#the text-tab delimited stats file
	    'out:s'  =>  \$out, #for output  
           );

if (($help)&&!($help)||!($out)||!($stats)){
 print "Usage : perl ParseStats.pl -stats  CD_T1_L1802_stats.txt -out CD_T1_L1802_mat.txt \n";
 print " -stats <txt> - the input text tab delimited stats file\n";
 print " -out <txt>  - the output matrix\n";
 print " -help        - Get this help\n";
 exit();
 }

# first go through the stats file and put all the sequence names that have the same forward and reverse into a hash
# @M01569:4:000000000-A383C:1:2109:4418:15978 1:N:0:2	Jseq_1-6	33	forward	0
# @M01569:4:000000000-A383C:1:2109:4418:15978 1:N:0:2	TRBV(12-3-5)	24	reverse	0
# @M01569:4:000000000-A383C:1:2109:21508:15979 1:N:0:2	TRBV2	31	forward	0

my (%info,%matrix,%xlab,%ylab);
my ($totalbases,$totalerror,$revbases,$reverror,$forbases,$forerror)=0;
open(STATS,$stats)||die "Can't open $stats\n";
while(<STATS>){
  chomp($_);
  my @cells=split(/\t/,$_);
  my $type=$cells[1];
  $type=~s/(\w{3,4}).+/$1/;
  my $nobracket=$cells[1];
  $nobracket=~s/\(//g;
  $nobracket=~s/\)//g;
  if (!$info{$cells[0]}){
    $info{$cells[0]}=$nobracket;
  }elsif ($info{$cells[0]}){
    if ($info{$cells[0]}<$nobracket){
      $matrix{$info{$cells[0]}}{$nobracket}++;
      $ylab{$info{$cells[0]}}++;
      $xlab{$nobracket}++;
    }else{
      $matrix{$cells[1]}{$info{$cells[0]}}++;
      $ylab{$nobracket}++;
      $xlab{$info{$cells[0]}}++;
    }
  }
  if ($cells[3]=~/forward/){
    $totalbases=$totalbases+$cells[2];
    $totalerror=$totalerror+$cells[4];
    $forbases=$forbases+$cells[2];
    $forerror=$forerror+$cells[4];
  }
  if ($cells[3]=~/reverse/){
    $totalbases=$totalbases+$cells[2];
    $totalerror=$totalerror+$cells[4];
    $revbases=$revbases+$cells[2];
    $reverror=$reverror+$cells[4];
  }  
}

print "Total bases = $totalbases\tTotal errors = $totalerror\tRatio = ".$totalerror/$totalbases."\n";
print "Forward bases = $forbases\tForward errors = $forerror\tRatio = ".$forerror/$forbases."\n";
print "Reverse bases = $revbases\tReverse errors = $reverror\tRatio = ".$reverror/$revbases."\n";


open (OUT,">$out")||die "Can't open $out\n";
my $firstline="labels\t";
foreach my $xlab (keys %xlab){
  $firstline=$firstline."$xlab\t";
}
$firstline=~s/\t$/\n/;
print OUT "$firstline";
foreach my $ylab (keys %ylab){
  my $newline="$ylab\t";
  foreach my $xlab (keys %xlab){
    if ($matrix{$ylab}{$xlab}=~/\d+/){
      $newline=$newline."$matrix{$ylab}{$xlab}\t";
    }else{
      # print "Nothing\n";
      $newline=$newline."-\t";
    }
  }
  $newline=~s/\t$/\n/;
  print OUT "$newline";
}