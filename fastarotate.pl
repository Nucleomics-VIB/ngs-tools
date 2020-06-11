#!/usr/bin/perl -w

# fastarotate.pl
# rotate a circular map with new +1 as user-provided position
# cut a fasta sequence at position <pos>
# paste the second part followed by the first
#
# Stephane Plaisance (VIB-NC+BITS) 2020/06/11; v1.0
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
getopts('i:b:o:h');
our ( $opt_i, $opt_b, $opt_o, $opt_h );
my $version="1.0, 2020_06_11";

my $usage = "## Usage: fastarotate.pl <-i fasta-file> <-b break position>
# script version ".$version."
# Additional optional parameters are:
# <-o prefix for output (default to <inputname>_rotate-<brk>.fasta)>
# <-h to display this help>";

defined($opt_h) && die $usage . "\n";

# handle IO
my $break = $opt_b || die $usage . "\n";
my $inpath = $opt_i || die $usage . "\n";
my @sufx = ( ".fa", ".fasta", ".fsa" );
my ($name,$path,$suffix) = fileparse($inpath,@sufx);
my $defout = $path."/".$name."_rotate-$break.fa";
my $outname = defined($opt_o)?$opt_o:$defout;
my $outpath = $path."/".$outname;

# IO handles
my $seqOUT = Bio::SeqIO->new(-file => ">$outpath", -format=>"Fasta");

my $seqIN = Bio::SeqIO->new(-file => "$inpath", -format=>"Fasta");
my $seq = $seqIN->next_seq();
my $len = $seq->length;

if ( $break gt $len ){
die "# break is further than the sequence end";
}

my $newseq = $seq->subseq($break+1,$len).$seq->subseq(1,$break);

my $seq_obj = Bio::Seq->new(-seq => $newseq,
                         -alphabet => 'dna',
                         -display_id => $seq->display_id(),
                         -desc => $seq->desc()." rotated_at_$break"
                         );
                         
                         
$seqOUT -> write_seq($seq_obj);

exit 0;
