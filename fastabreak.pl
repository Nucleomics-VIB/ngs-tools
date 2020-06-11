#!/usr/bin/perl -w

# fastabreak.pl
# prepare file for circlator trimming
# cut a fasta sequence at position <pos>
# insert a new header '>Break'
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

my $usage = "## Usage: fastabreak.pl <-i fasta-file> <-b break position>
# script version ".$version."
# Additional optional parameters are:
# <-o prefix for output (default to <inputname>_break-<brk>.fasta)>
# <-h to display this help>";

defined($opt_h) && die $usage . "\n";

# handle IO
my $break = $opt_b || die $usage . "\n";
my $inpath = $opt_i || die $usage . "\n";
my @sufx = ( ".fa", ".fasta", ".fsa" );
my ($name,$path,$suffix) = fileparse($inpath,@sufx);
my $defout = $path."/".$name."_break-$break.fa";
my $outname = defined($opt_o)?$opt_o:$defout;
my $outpath = $path."/".$outname;

# IO handles
my $seqOUT = Bio::SeqIO->new(-file => ">$outpath", -format=>"Fasta");

my $seqIN = Bio::SeqIO->new(-file => "$inpath", -format=>"Fasta");
my $seq = $seqIN->next_seq();
my $len = $seq->length;

# test end
if ( $break > $len ){
die "# break is further than the sequence end";
}

my $seq_obj1 = Bio::Seq->new(-seq => $seq->subseq(1,$break),
                         -alphabet => 'dna',
                         -display_id => $seq->display_id(),
                         -desc => $seq->desc()." break_at_$break"
                         );
                         
                         
$seqOUT -> write_seq($seq_obj1);

my $seq_obj2 = Bio::Seq->new(-seq => $seq->subseq($break+1,$len),
                         -alphabet => 'dna',
                         -display_id => "Break" );
                         
                         
$seqOUT -> write_seq($seq_obj2);

exit 0;
