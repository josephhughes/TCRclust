#!/usr/bin/perl

# Perl script written by Joseph Hughes, University of Glasgow
# use this to remove the fastq sequences that have the same forward and reverse primer 
# provide a fastq and stats input, give a subset fastq and stats output

use strict;
use Getopt::Long; 
use IO::Handle;
use Carp;
# subroutines 
sub open_and_detect_input_format;
sub create_output_files;
sub read_record;
sub write_record($);
# global variables
my $input_file_io;
my $fastq_format = 1;
my $debug = 1 ;
my %files;
my $newfile_suffix='.fq';
my $newfile_prefix  ;
my %filenames;
# the 4 lines per record in FASTQ format.
# (in FASTA format, only the first 2 are used)
my $seq_name;
my $seq_bases;
my $seq_name2;
my $seq_qualities;

my ($stats,$stub,$help);
&GetOptions(
	    'stats:s'      => \$stats,#the text-tab delimited stats file
	    'stub:s'  =>  \$stub, #stub for output  
           );

if (($help)&&!($help)||!($stub)||!($stats)){
 print "Usage : cat CD_T1_L1802_trim.fq | perl RemoveSamePrimer.pl -stats  CD_T1_L1802_stats.txt -stub CD_T1_L1802_new_ \n";
 print " -stats <txt> - the input text tab delimited stats file\n";
 print " -stub <txt>  - the output stub for the new .fq and stats file\n";
 print " -help        - Get this help\n";
 exit();
 }

# first go through the stats file and put all the sequence names that have the same forward and reverse into a hash
# @M01569:4:000000000-A383C:1:2109:4418:15978 1:N:0:2	Jseq_1-6	33	forward	0
# @M01569:4:000000000-A383C:1:2109:4418:15978 1:N:0:2	TRBV(12-3-5)	24	reverse	0
# @M01569:4:000000000-A383C:1:2109:21508:15979 1:N:0:2	TRBV2	31	forward	0

my (%info,%types);
open(STATS,$stats)||die "Can't open $stats\n";
while(<STATS>){
  chomp($_);
  my @cells=split(/\t/,$_);
  my $type=$cells[1];
  $type=~s/(\w{3,4}).+/$1/;
  $info{$cells[0]}{$type}=$_;
  $types{$type}++;
}

# go through the fastq and only printout those that are not found in the hash
# print out a new stats file with the subset
open_and_detect_input_format;
create_output_files;

open (OUT,">$stub\stats.txt")||die "Can't open $stub\stats.txt\n";
my $cnt=0;
my $kept=0;
my $seqs=0;
while ( read_record ) {
  chomp($seq_name);
  chomp($seq_bases);
  chomp($seq_name2);
  chomp($seq_qualities);
  if (!%{$info{$seq_name}}){
    print "$seq_name doesn't exist in stats\n";
  }else{
	  my $nb = keys %{$info{$seq_name}};
	  if ($nb<=1){
		# print $seq_name;
		$cnt++;
	  }elsif ($nb>1){
		my $file = $files{"trim"};
		write_record($file);
		$seqs++;
		foreach my $type (keys %{$info{$seq_name}}){
		  print OUT "$info{$seq_name}{$type}\n";
		  $kept++;
		}
	  }
  }
}
print "Total number of fastq removed: $cnt\n"; 
print "Total number kept: ".($kept/2)."\n";  
# print "Total number of sequences: ".$seqs."\n"; 
########### subroutines #################

sub create_output_files {
    # %labels should be a unique list of identifiers
	my @exceptions=qw/trim/;
	
   foreach my $except (@exceptions){
     my $new_filename = $stub . $except . $newfile_suffix; 
	 $filenames{$except} = $new_filename;
	 open my $file, ">$new_filename" or die "Error: failed to create output file ($new_filename)\n"; 
	 $files{$except} = $file ;
    }
}

sub read_record
{
	$seq_name = $input_file_io->getline();
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



