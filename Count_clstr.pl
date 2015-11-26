#!/usr/bin/perl -w

# Perl script written by Joseph Hughes, University of Glasgow
# use this to compare the number of sequences in different clusters for different samples
# create a table with the relative counts for each, the cluster name and the representative from the cluster

use strict;
use Getopt::Long; 
use Bio::SeqIO;

my ($infile,$outfile,$clstr,%count,%clusters,$help,$ids);
&GetOptions(
	    'i:s'      => \$infile,#fasta file
	    'clstr:s'  =>\$clstr, #a cd-hit generated cluster file
	    'ids:s'    =>\$ids, #a list of identifiers comma separated 
	    'o:s'         => \$outfile, # use this to output a text tab-delimited file with the counts from representative sequences
           );

if (($help)||!$infile||!$outfile){
  print "usage: Count_clstr.pl [-h] [-i INFILE] [-o OUTFILE] [-clstr <txt>] \n";
  print "by J. Hughes (joseph(dot)hughes(at)glasgow(dot)ac(dot)uk\n\n";

  print "   [-h]         = This helpful help screen.\n";
  print "   [-clstr <txt>]       = the cd-hit cluster file .clstr\n";
  print "   [-ids <txt>]       = a list of comma separated identifiers\n";
  print "   [-i INFILE]  = FASTA input file.\n";
  print "   [-o OUTFILE] = the name of the test file to output\n";
  exit();
}

print "$clstr\n";
my @ids=split(/,/,$ids);
print "@ids\n";
my $seq_in  = Bio::SeqIO->new(-format => 'fasta',-file   => $infile);
my $seqcnt=0;
while(my $seq = $seq_in->next_seq() ) {
    my $desc = $seq->desc();
    my $id = $seq->display_id();
    my $seq_str= $seq->seq();
    $count{$id}=$desc;
    #print "$id\n";
    $seqcnt++;
}
open(CLUSTER,"<$clstr")||die "Can't open $clstr\n";
my $clusterid;
my $longest;
my %longest;

open(OUT,">$outfile")||die "Can't open $outfile\n";
while(<CLUSTER>){
  if ($_=~/^>(.+)$/){
    $clusterid=$1;
    $clusterid=~s/\s//g;
   # print "$clusterid\n";
  }elsif ($_=~/\d+.+\>(.+)\.\.\..+\*$/){
    my $id=$1;
    $longest=$id;
    my $sample_name=$id;
    $sample_name=~s/_\d+$//;
    #print "$sample_name\n";
    $clusters{$clusterid}{$sample_name}=$longest;
    #print "$clusterid\t$id\t$longest\n";
    $longest{$clusterid}=$longest;
  }
  elsif ($_=~/\d+.+\>(.+)\.\.\..+$/){
    my $id=$1;
    my $sample_name=$id;
    $sample_name=~s/_\d+$//;
   # print "$sample_name\n";
    $clusters{$clusterid}{$sample_name}=$id;
    #print "$clusterid\t$id\t$longest\n";
  }
}

foreach my $clusterid (keys %clusters){
   my $nbseqs=keys %{$clusters{$clusterid}};
   print OUT "$clusterid\t$longest{$clusterid}\t$nbseqs\t";
   my $results="";
   foreach my $sample (@ids){
     my $repid=$clusters{$clusterid}{$sample};
     if ($repid){
      # print "$repid\n";
       $results=$results."$clusters{$clusterid}{$sample}\t$count{$repid}\t";
     }else{
       $results=$results."NA\tNA\t";
     }
   }
   $results=~s/\t$//;
   print OUT "$results\n";
}
