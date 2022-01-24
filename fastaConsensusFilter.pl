#!/usr/bin/perl -w

# script: fastaConsensusFilter.pl
#
# filter records based on the proportion of capital letters
# more capital letters means better consensus
# save only filtered subset to a new file
#
# Stephane Plaisance (VIB-NC+BITS) 2020/06/11; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

use warnings;
use strict;
use Bio::SeqIO;
use File::Basename;
use Getopt::Std;

# disable buffering to get output during long process (loop)
$|=1; 

############################
# handle command parameters
############################
getopts('i:m:o:h');
our ( $opt_i, $opt_m, $opt_o, $opt_h );
my $version="1.0, 2021_12_30";

my $usage = "## Usage: fastaConsensusFilter.pl <-i fasta-file> <-m minimum % capital Bases>
# script version ".$version."
# Additional optional parameters are:
# <-o prefix for output (default to <inputname>_gt<min>pc.fa)>
# <-h to display this help>";

defined($opt_h) && die $usage . "\n";

# handle IO
my $inpath = $opt_i || die $usage . "\n";
my @sufx = ( ".fa", ".fasta", ".fsa" );
my ($name,$path,$suffix) = fileparse($inpath,@sufx);
my $min = $opt_m || die $usage . "\n";

my $defout = $path."/".$name."_gt".$min."pc.fa";
my $outname = defined($opt_o)?$opt_o:$defout;
my $outpath = $path."/".$outname;
(my $outpathbad = $outpath) =~ s/\.fa$/_bad\.fa/;

# IO handles
my $seqIN = Bio::SeqIO->new(-file => "$inpath", -format=>"Fasta");
my $seqOUT = Bio::SeqIO->new(-file => ">$outpath", -format=>"Fasta");
my $seqOUTbad = Bio::SeqIO->new(-file => ">$outpathbad", -format=>"Fasta");

# count records and show progress
my $counter=0;

while (my $seq = $seqIN->next_seq) {

$counter ++;
# tell every 10000
$counter =~ m/0000$/ and print STDERR "# processed $counter sequences\r";

my $len = $seq->length;
my $uc = () = $seq->seq =~ m/\p{Uppercase}/g;
my $frac = 100*$uc/$len;

# test proportion
my $result = ( $frac >= $min ) ? $seqOUT -> write_seq($seq) : $seqOUTbad -> write_seq($seq); 

}

exit 0;
