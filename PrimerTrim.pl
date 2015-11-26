#!/usr/bin/perl

# Perl script written by Joseph Hughes, University of Glasgow
# use this script to provide count differences to primer and  
# trim the primer
# and sort out no primer matches

use strict;
use Getopt::Long; 
use POSIX;
use IO::Handle;
use Carp;
use Bio::SeqIO;
# subroutines 
sub open_and_detect_input_format;
sub create_output_files;
sub read_record;
sub write_record($);
sub match_sequences ;
# Global variables 
# Populated by 'create_output_files'
my $input_file_io;
my $fastq_format = 1;
my $debug = 1 ;
my $file;
my $progress=0;
my $alength=13;
# The Four lines per record in FASTQ format.
# (when using FASTA format, only the first two are used)
my $seq_name;
my $seq_bases;
my $seq_name2;
my $seq_qualities;

my ($in, $filename, $primerfile,$help, $i,%primer,$stats,%files);

&GetOptions(
	    'in:s'      => \$in,#fastafile or fastq
	    'out:s'   => \$filename,#output fasta file or fastq
	    'primer:s'   => \$primerfile,#file with primers	
	    'stats:s'   => \$stats,
	    'alength:i'  =>  \$alength, #length of the adaptor  
           );

if (($help)&&!($help)||!($filename)||!($in)){
 print "Usage : cat 2108_join.fq | perl PrimerTrim.pl -in 2108_join.fq -out 2108_join_trim.fq -primer primers_58.fa -alength 13 -stats 2108_stats.txt \n";
 print " -in <txt> - the input fasta or fastq file\n";
 print " -out <txt> - the name of your output fasta of fastq file\n";
 print " -primer <txt> - the primer sequences in fasta\n";
 print " -stats <txt> - output statistics\n";
 print " -alength <integer>  - length in nucleotides of the adaptor (default = 13)\n";
 print " -help        - Get this help\n";
 exit();
 }

# grabs the FASTa parser
# my $primerseq = Bio::SeqIO->new(-format    => 'fasta',
#                            -file      => $primerfile);
# 
# while (my $pseq = $primerseq->next_seq) {
#     my $pid=$pseq->id();
#      my $pseqstr=$pseq->seq();
#      $primer{$pid}=$pseqstr;
# }
open (PRIMER,"$primerfile")||die "Can't open $primerfile\n";
while(<PRIMER>){
  chomp($_);
  if ($_=~/>(\S+)/){
    print $1,"\n";
    my $nextline = <PRIMER>;
    chomp($nextline);
    print "$nextline\n";
    $primer{$1}=$nextline;
  }
}


create_output_files;
open_and_detect_input_format;
my $total;
if ($fastq_format){
  $total = `grep -c "^+" $in 2>&1`;
  $total++;
  print "Total number of fastq sequences: $total\n";
}else{
  $total = `grep -c "^>" $in 2>&1`;
  $total++;
  print "Total number of fasta sequences: $total\n";

}
open(OUT,">$stats")||die "Can't open output $stats\n";

match_sequences ;


sub match_sequences {

while ( read_record ) {
  my $percent=$progress*100/$total;
  printf STDOUT "Progress %03d%%\r",$percent;
  chomp($seq_name);
  chomp($seq_bases);
  chomp($seq_name2);
  chomp($seq_qualities);
  #use the unmatched substrings to look for best barcode match
  #first check for exact match on the forward within the first 10 bases
  #save the no matches to unmatched
  #then on the reverse use the match function allowing for 1 mismatch to get the longest substring
  my $best_fprimer_mmcount = 30;
  my $best_fprimer_ident = undef;
  my $best_rprimer_mmcount = length($seq_name);
  my $best_rprimer_ident = undef;
   foreach my $ident (keys %primer) {
     my $primer_length=length($primer{$ident});
     my $rprimer =	revcompl($primer{$ident});
     my $bsequence_fragment = substr $seq_bases, 0, ($primer_length+$alength+1);
     my $esequence_fragment = substr $seq_bases, -($primer_length+$alength);
     my ($align1,$align2,$begseq,$endseq)=&SmithWaterman($primer{$ident},$bsequence_fragment);#use length of $raling2 and $begseq to trim the begining
     my ($ralign1,$ralign2,$rbegseq,$rendseq)=&SmithWaterman($rprimer,$esequence_fragment);#use length of $ralign2 and $rendseq to trim the end
     #print "New\n$esequence_fragment\n$ralign1\n$ralign2\n$rbegseq\n$rendseq\n";
     #print "New\n$bsequence_fragment\n$align1\n$align2\n$begseq\n$endseq\n";
     my $mm = mismatch_count($align1, $align2) ;
     my $rmm = mismatch_count($ralign1, $ralign2) ;
     #print "Mismatch $mm\n"; 
     if ($mm<$best_fprimer_mmcount){
        $best_fprimer_mmcount = $mm;
        $best_fprimer_ident = "$ident,$primer{$ident},$align1,$align2,$begseq,$endseq,$best_fprimer_mmcount";
        #print "$best_fprimer_ident\n";
     }
    if ($rmm<$best_rprimer_mmcount){
        $best_rprimer_mmcount = $rmm;
        $best_rprimer_ident = "$ident,$primer{$ident},$ralign1,$ralign2,$rbegseq,$rendseq,$best_rprimer_mmcount";
        #print "$best_fprimer_ident\n";
     }     
     
     
   }
  #print "BestF: $best_fprimer_ident\nBestR: $best_rprimer_ident\n";
  my @fres=split(/,/,$best_fprimer_ident);
  my @rres=split(/,/,$best_rprimer_ident);
  if (($fres[6]/(length($fres[3])))<0.7 && ($rres[6]/(length($rres[3])))<0.7){ 
	print OUT "$seq_name\t$fres[0]\t".length($fres[3])."\tforward\t$fres[6]\n";
	print OUT "$seq_name\t$rres[0]\t".length($rres[3])."\treverse\t$rres[6]\n";
	#print "$seq_bases\n$seq_qualities\n";
	my $begtrim= length($fres[3]) + length ($fres[4]);
	my $endtrim=  length($rres[3]) + length ($rres[5]);
	$seq_bases=substr $seq_bases, $begtrim, - $endtrim;
	$seq_qualities=substr $seq_qualities, $begtrim, - $endtrim;
	my $file=$files{$filename};
	write_record($file);
  }else{
	my $file=$files{"error"};
	write_record($file);
  }
}
}

sub revcompl {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

sub create_output_files {
    # %labels should be a unique list of identifiers
	my @exceptions=qw/error/;
	
   foreach my $except (@exceptions){
     my $new_filename = $filename;
     $new_filename=~s/\.fq/_error.fq/; 
	 open my $file, ">$new_filename" or die "Error: failed to create output file ($new_filename)\n"; 
	 $files{$except} = $file ;
    }
    
	my $new_filename = $filename; 
	open my $file, ">$new_filename" or die "Error: failed to create output file ($new_filename)\n"; 
	$files{$filename} = $file ;
}

sub read_record
{
	$seq_name = $input_file_io->getline();
    $progress++;
	return undef unless defined $seq_name; # End of file?

	$seq_bases = $input_file_io->getline();
	die "Error: bad input file, expecting line with sequences\n" unless defined $seq_bases;

	# If using FASTQ format, read two more lines
	if ($fastq_format) {
		$seq_name2  = $input_file_io->getline();
		die "Error: bad input file, expecting line with sequence name2\n" unless defined $seq_name2;

		$seq_qualities = $input_file_io->getline();
		die "Error: bad input file, expecting line with quality scores\n" unless defined $seq_qualities;
	}
	return 1;
}

sub write_record($)
{
	my $file = shift;

	croak "Bad file handle" unless defined $file;

	print $file $seq_name,"\n";
	print $file $seq_bases,"\n";
	#if using FASTQ format, write two more lines
	print $file $seq_name2,"\n";
    print $file $seq_qualities, "\n";
}

sub open_and_detect_input_format
{
	$input_file_io  = new IO::Handle;
	die "Failed to open STDIN " unless $input_file_io->fdopen(fileno(STDIN),"r");

	# Get the first characeter, and push it back
	my $first_char = $input_file_io->getc();
	$input_file_io->ungetc(ord $first_char);

	if ($first_char eq '>') {
		# FASTA format
		$fastq_format = 0 ;
		print STDERR "Detected FASTA format\n" if $debug;
	} elsif ($first_char eq '@') {
		# FASTQ format
		$fastq_format = 1;
		print STDERR "Detected FASTQ format\n" if $debug;
	} else {
		die "Error: unknown file format. First character = '$first_char' (expecting > or \@)\n";
	}
}

#Quickly calculate hamming distance between two strings
#
#NOTE: Strings must be same length.
#      returns number of different characters.
#see  http://www.perlmonks.org/?node_id=500235
sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }

sub SmithWaterman{
my ($seq1,$seq2)=@_;
#print "IN $seq1 $seq2\n";
# scoring scheme
my $MATCH     =  1; # +1 for letters that match
my $MISMATCH = 0; # -1 for letters that mismatch
my $GAP       = -1; # -1 for any gap

# initialization
my @matrix;
$matrix[0][0]{score}   = 0;
$matrix[0][0]{pointer} = "none";
for(my $j = 1; $j <= length($seq1); $j++) {
     $matrix[0][$j]{score}   = 0;
     $matrix[0][$j]{pointer} = "none";
}
for (my $i = 1; $i <= length($seq2); $i++) {
     $matrix[$i][0]{score}   = 0;
     $matrix[$i][0]{pointer} = "none";
}

# fill
 my $max_i     = 0;
 my $max_j     = 0;
 my $max_score = 0;

 for(my $i = 1; $i <= length($seq2); $i++) {
     for(my $j = 1; $j <= length($seq1); $j++) {
         my ($diagonal_score, $left_score, $up_score);
         
         # calculate match score
         my $letter1 = substr($seq1, $j-1, 1);
         my $letter2 = substr($seq2, $i-1, 1);      
         if ($letter1 eq $letter2) {
             $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
          }
         else {
             $diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
          }
         
         # calculate gap scores
         $up_score   = $matrix[$i-1][$j]{score} + $GAP;
         $left_score = $matrix[$i][$j-1]{score} + $GAP;
         
         if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
             $matrix[$i][$j]{score}   = 0;
             $matrix[$i][$j]{pointer} = "none";
             next; # terminate this iteration of the loop
          }
         
         # choose best score
         if ($diagonal_score >= $up_score) {
             if ($diagonal_score >= $left_score) {
                 $matrix[$i][$j]{score}   = $diagonal_score;
                 $matrix[$i][$j]{pointer} = "diagonal";
              }
             else {
                 $matrix[$i][$j]{score}   = $left_score;
                 $matrix[$i][$j]{pointer} = "left";
              }
          } else {
             if ($up_score >= $left_score) {
                 $matrix[$i][$j]{score}   = $up_score;
                 $matrix[$i][$j]{pointer} = "up";
              }
             else {
                 $matrix[$i][$j]{score}   = $left_score;
                 $matrix[$i][$j]{pointer} = "left";
              }
          }
         
       # set maximum score
         if ($matrix[$i][$j]{score} > $max_score) {
             $max_i     = $i;
             $max_j     = $j;
             $max_score = $matrix[$i][$j]{score};
          }
      }
 }

 # trace-back

 my $align1 = "";
 my $align2 = "";
 my $foverhang = "";
 my $roverhang = "";
#print "$max_j\n$max_i\n";
 my $j = $max_j;
 my $i = $max_i;
if ($j<$i){
  ###$align1 .= reverse(substr($seq1, $j, length($seq1)-$j));#uncomment ### if you want full alignment with overhangs
  ###$align2 .= reverse(substr($seq2, $i, length($seq2)-$i));
  $roverhang = substr($seq2, $i, length($seq2)-$i);
#  print "Number of gaps to add " , (length($align2)-length($align1)), "\n";
  my $nbgaps=(length($align2)-length($align1));
  if ($nbgaps>0){
    ###$align1 .= "-" x (length($align2)-length($align1));
  }elsif ($nbgaps<0){
    ###$align2 = "-" x (length($align1)-length($align2));
    ###$align2 .= reverse(substr($seq2, $i, length($seq2)-$i));
  }
}

 while (1) {
     if ($matrix[$i][$j]{pointer} eq "none"){
       #print "Min $i, $j\n";
       if ($j<$i){
         ###$align1 .= "-" x $i;
         ###$align2 .= reverse(substr($seq2, $j, $i));
         $foverhang = substr($seq2, $j, $i);
       }elsif ($i<$j){# this hasn't been checked with data yet
         ###$align1 .= reverse(substr($seq2, $i, $j));
         ###$align2 .= "-" x $j;
       }  
       last;
     }  
     if ($matrix[$i][$j]{pointer} eq "diagonal") {#rather than trim substring, I want to do padding
         $align1 .= substr($seq1, $j-1, 1);
         $align2 .= substr($seq2, $i-1, 1);
         $i--; $j--;
      }
     elsif ($matrix[$i][$j]{pointer} eq "left") {
         $align1 .= substr($seq1, $j-1, 1);
         $align2 .= "-";
         $j--;
      }
     elsif ($matrix[$i][$j]{pointer} eq "up") {
         $align1 .= "-";
         $align2 .= substr($seq2, $i-1, 1);
         $i--;
      }  
 }

 $align1 = reverse $align1;
 $align2 = reverse $align2;
 return ($align1,$align2,$foverhang,$roverhang);
 }