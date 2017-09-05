#!/usr/bin/perl -w

# fasta2merge.pl
# Merge all sequences from a multiFasta file into a single fasta record
# required for the abacas.pl -r input 
#
# Stephane Plaisance (VIB-NC+BITS) 2017/09/05; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

use warnings;
use strict;
use Bio::SeqIO;
use File::Basename;
use Getopt::Std;

############################
# handle command parameters
############################
getopts('i:o:h');
our ( $opt_i, $opt_o, $opt_h );
my $version="1.0, 2017_09_05";

my $usage = "## Usage: fasta2merge.pl <-i fasta-file>
# script version ".$version."
# Additional optional parameters are:
# <-o prefix for output (default to merged-<inputname>.fasta)>
# <-h to display this help>";

my $infile = $opt_i || die $usage . "\n";
my $outfile = defined($opt_o)?$opt_o:"merged-".$infile;
defined($opt_h) && die $usage . "\n";

# handle IO
my $inpath = dirname($infile);
my @sufx = ( ".fa", ".fasta", ".fsa" );
my $name = basename( $infile, @sufx );
my $outpath = $inpath."/".$outfile;

# IO handles
my $seqIN = Bio::SeqIO->new(-file => $infile, -format=>"Fasta");
my $seqOUT = Bio::SeqIO->new(-file => ">$outpath", -format=>"Fasta");

# loop through sequences and search motifs
my $mergedname="merged_sequences_from_".$name;
my $mergedseq="";

while ( my $seq = $seqIN->next_seq() ) {
	my $title = $seq->id;
	my $chrlen = $seq->length;
	print STDOUT "# adding ".$title."\t".$chrlen."\n";
	$mergedseq .= $seq->seq();
	}

my $seq_obj = Bio::Seq->new(-seq => $mergedseq,
                         -alphabet => 'dna',
                         -display_id => $mergedname );
                         
$seqOUT -> write_seq($seq_obj);

exit 0;
