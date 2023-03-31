#!/usr/bin/perl -w

# joinFastqPair.pl
# Bookend join reads in pairs to produce single reads for FastQC
#
# Stephane Plaisance - VIB Nucleomics Core
# version 1.0, 2023-03-31
# visit our Git: https://github.com/Nucleomics-VIB

use warnings;
use strict;
use FASTX::Reader;
use File::Basename;
use Getopt::Std;

############################
# handle command parameters
############################
getopts('i:j:o:h');
our ( $opt_i, $opt_j, $opt_o, $opt_h );
my $version="1.0, 2023_03_31";

my $usage = "## Usage: joinFastqPair.pl -i <first-read> -j <second-read>
# script version ".$version."
# Additional optional parameters are:
# <-o prefix for output (default to merged-<first-read>)>
# <-h to display this help>";

my $first = $opt_i || die $usage . "\n";
my $second = $opt_j || die $usage . "\n";
defined($opt_h) && die $usage . "\n";

# handle IO
my $inpath = dirname($first);
my @sufx = ( ".fq", ".fastq", ".fq.gz", ".fastq.gz" );
my $name = basename( $first, @sufx );
my $outfile = defined($opt_o)?$opt_o:"merged-".$name.".fq.gz";
my $outpath = $inpath."/".$outfile;

# IO handles
my $seqIN1 = FASTX::Reader->new({ filename => "$first" });
my $seqIN2 = FASTX::Reader->new({ filename => "$second" });
open (my $seqOUT, "| gzip -c > $outpath") or die "error starting gzip $!";

# loop through sequences and search motifs
my $mergedseq="";
my $mergedqual="";
my $counter=0;
my $saved=0;

while ( my $seq1 = $seqIN1->getRead() ) {

        # tell every 100000
        $counter++;
        $counter =~ m/00000$/ and print STDERR ".\r";

        my $seq2 = $seqIN2->getRead();
        my $title = "@".$seq1->{name}." ".$seq1->{comment};
        my $newseq = $seq1->{seq}.$seq2->{seq};
        my $newlen = length($newseq);
        my $newqual = $seq1->{qual}.$seq2->{qual};
        
        # only keep full length pairs to avoid async in FastQC
        if ( length($seq1->{seq}) == length($seq2->{seq}) ){
                print $seqOUT $title, "\t", "merged-pair of length ".$newlen, "\n", $newseq, "\n+\n", $newqual, "\n";
                $saved++;
                }
        }

print STDERR "\n# ".$saved." joined read pairs, out of ".$counter." input pairs, were written to ".$outpath."\n";

close $seqOUT;
