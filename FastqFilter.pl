#!/usr/bin/perl

# Perl script written by Joseph Hughes, University of Glasgow
# use this sorts the sequences according to specified min and max sequence length cut-off  
# DOES NOT WORK FOR MULTILINE FASTQ

use strict;
use Getopt::Long; 
use POSIX;
use IO::Handle;
use Carp;
# subroutines 
sub open_and_detect_input_format;
sub create_output_files;
sub read_record;
sub write_record($);
sub filter_sequences ;
# Global variables 
# Populated by 'create_output_files'
my $input_file_io;
my $fastq_format = 1;
my $debug = 1 ;
my $progress=0;
# The Four lines per record in FASTQ format.
# (when using FASTA format, only the first two are used)
my $seq_name;
my $seq_bases;
my $seq_name2;
my $seq_qualities;

my ($filename,$help, $i,%files,$min,$max);

&GetOptions(
	    'stub:s'   => \$filename,#output fasta file or fastq
	    'min:i'  =>  \$min, #minimum length of sequences to keep
	    'max:i'  =>  \$max, #max length of sequences to keep  
           );

if (($help)&&!($help)||!($filename)){
 print "Usage : cat 2108_join.fq | perl FastqFilter.pl -stub 2108_sort -min 40 -max 90 \n";
 print " -stub <txt> - the stub name to use for the output fasta of fastq file: _ok, _above90 _below40\n";
 print " -min <integer> - minimum length in bases\n";
 print " -max <txt> - maximum length in bases\n";
 print " -help        - Get this help\n";
 exit();
 }

create_output_files;
open_and_detect_input_format;
filter_sequences ;

my $mincnt=0;
my $maxcnt=0;
my $okcnt=0;
sub filter_sequences {

	while ( read_record ) {
	  chomp($seq_name);
	  chomp($seq_bases);
	  chomp($seq_name2);
	  chomp($seq_qualities);
	  if (length($seq_bases)>$max){
		my $file=$files{"_above$max"};
		write_record($file);
		$maxcnt++;
	  }elsif (length($seq_bases)<$min){
		my $file=$files{"_below$min"};
		write_record($file);
		$mincnt++;
	  }else{
		my $file=$files{"_ok"};
		write_record($file);
		$okcnt++;
	  }
	}
	print "Final summary:\nNumber of sequences below $min: $mincnt\n";
    print "Number of sequences above $max: $maxcnt\n";
    print "Number meeting the criteria: $okcnt\n";
}



sub create_output_files {
    # %labels should be a unique list of identifiers
	my @exceptions=("_ok","_below$min","_above$max");
	
   foreach my $except (@exceptions){
     my $new_filename = $filename.$except.".fq";
	 open my $file, ">$new_filename" or die "Error: failed to create output file ($new_filename)\n"; 
	 $files{$except} = $file ;
    }
    
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

