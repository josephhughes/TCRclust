#!/usr/bin/perl -w

# Perl script written by Joseph Hughes, University of Glasgow
# use this to sort through the output from cd-hit for example remove singletons
# or pull out representative sequences


use strict;
use Getopt::Long; 
use Bio::SeqIO;

my ($infile,$outfile,$clstr,$freq,$singleton,%sequences,%clusters,$singletons,$rep,$help);
&GetOptions(
	    'i:s'      => \$infile,#fasta file
	    'clstr:s'  =>\$clstr, #a cd-hit generated cluster file
	    'o:s'         => \$outfile, # use this to output a fasta formatted file
	    'freq:s'   => \$freq, #the frequency cut-off
	    'singleton' => \$singletons, #remove singletons
	    'rep' => \$rep, #get rep sequence
           );

if (($help)||!$infile||!$outfile){
  print "usage: sort-cdhit.pl [-h] [-i INFILE] [-o OUTFILE] [-clstr <txt>] [-freq <txt>] [-singleton <txt>] [-rep <txt>]\n";
  print "by J. Hughes (joseph(dot)hughes(at)glasgow(dot)ac(dot)uk\n\n";

  print "   [-h]         = This helpful help screen.\n";
  print "   [-clstr <txt>]       = the cd-hit cluster file .clstr\n";
  print "   [-freq <txt>]       = the frequency cut-off for a cluster to be kept\n";
  print "   [-singleton <txt>]       = to sort out the singletons versus non-singletons\n";
  print "   [-rep <txt>]       = to pull out the representative sequence\n";
  print "   [-i INFILE]  = FASTA input file.\n";
  print "   [-o OUTFILE] = stub for FASTA output file.\n";
  exit();
}

print "$clstr\n";
my $seq_in  = Bio::SeqIO->new(-format => 'fasta',-file   => $infile);
my $seqcnt=0;
while(my $seq = $seq_in->next_seq() ) {
    my $desc = $seq->desc();
    my $id = $seq->display_id();
    my $seq_str= $seq->seq();
    $sequences{$id}=$seq_str;
    #print "$id\n";
    $seqcnt++;
}
open(CLUSTER,"<$clstr")||die "Can't open $clstr\n";
my $clusterid;
my $longest;
my %longest;

while(<CLUSTER>){
  if ($_=~/^>(.+)$/){
    $clusterid=$1;
    $clusterid=~s/\s//g;
   # print "$clusterid\n";
  }elsif ($_=~/\d+.+\>(.+)\.\.\..+\*$/){
    my $id=$1;
    $longest=$id;
    $clusters{$clusterid}{$id}=$sequences{$longest};
    #print "$clusterid\t$id\t$longest\n";
    $longest{$clusterid}=$longest;
  }
  elsif ($_=~/\d+.+\>(.+)\.\.\..+$/){
    my $id=$1;
    $clusters{$clusterid}{$id}=$sequences{$id};
    #print "$clusterid\t$id\t$longest\n";
  }
}

if ($singletons){
  my $nosingletons = Bio::SeqIO->new(-file => ">$outfile.nosingletons.fa" , '-format' => 'fasta');
  my $singletons = Bio::SeqIO->new(-file => ">$outfile.singletons.fa" , '-format' => 'fasta');
  my $repfile = Bio::SeqIO->new(-file => ">$outfile.nosingletons.rep.fa", '-format' => 'fasta');
  foreach my $cluster (keys %clusters){
    my $nbseqs=keys %{$clusters{$cluster}};
    if ($nbseqs>1){
    #  print "$cluster $nbseqs\n";
      foreach my $id (keys %{$clusters{$cluster}}){
        my $seq_obj = Bio::Seq->new(-seq => $clusters{$cluster}{$id}, 
                          -display_id => $id);#-desc => $nbseqs, 
        $nosingletons->write_seq($seq_obj);
      }
     if ($rep){
        my $nbseqs=keys %{$clusters{$cluster}};
        my $rep_obj = Bio::Seq->new(-seq=> $sequences{$longest{$cluster}},
                                    -display_id=> $longest{$cluster}, -desc => "$nbseqs");
        $repfile->write_seq($rep_obj);
      }
    }elsif ($nbseqs<=1){
      foreach my $id (keys %{$clusters{$cluster}}){
         my $seq_obj = Bio::Seq->new(-seq => $clusters{$cluster}{$id},
                          -display_id => $id);# -desc => $nbseqs, 
        $singletons->write_seq($seq_obj);
      }
    }
  }
}if ($freq){
  my $above = Bio::SeqIO->new(-file => ">$outfile.above$freq.fa" , '-format' => 'fasta');
  my $below = Bio::SeqIO->new(-file => ">$outfile.below$freq.fa" , '-format' => 'fasta');
  my $repfile = Bio::SeqIO->new(-file => ">$outfile.rep$freq.fa", '-format' => 'fasta');
  foreach my $cluster (keys %clusters){
    my $nbseqs=keys %{$clusters{$cluster}};
    my $seqfreq=$nbseqs/$seqcnt;
    #print "$nbseqs $seqfreq\n";
    if ($seqfreq>$freq){
      foreach my $id (keys %{$clusters{$cluster}}){
         my $seq_obj = Bio::Seq->new(-seq => $clusters{$cluster}{$id}, 
                          -display_id => $id);#-desc => $nbseqs, 
        $above->write_seq($seq_obj);
      }
      if ($rep){
        my $nbseqs=keys %{$clusters{$cluster}};
        my $rep_obj = Bio::Seq->new(-seq=> $sequences{$longest{$cluster}},
                                    -display_id=> $longest{$cluster}, -desc => "$nbseqs");
        $repfile->write_seq($rep_obj);
      }
    }elsif ($seqfreq<=$freq){
      foreach my $id (keys %{$clusters{$cluster}}){
         my $seq_obj = Bio::Seq->new(-seq => $clusters{$cluster}{$id},
                          -display_id => $id);# -desc => $nbseqs, 
        $below->write_seq($seq_obj);
      }
    }
  }
}if ($rep && !$freq && !$singletons){
my $repseqs = Bio::SeqIO->new(-file => ">$outfile.rep.fa" , '-format' => 'fasta');
  foreach my $clusterid (keys %longest){
	my $nbseqs=keys %{$clusters{$clusterid}};
	my $seq_obj = Bio::Seq->new(-seq => $sequences{$longest{$clusterid}}, 
							      -display_id => $longest{$clusterid}, -desc => "$nbseqs"); 
	$repseqs->write_seq($seq_obj);
  }
}




