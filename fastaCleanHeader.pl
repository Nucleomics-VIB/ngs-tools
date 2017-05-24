#!/usr/bin/perl -w

# Remove spaces from fasta headers in a  multifasta
#
# Stephane Plaisance (VIB-NC+BITS) 2017/05/19; v1.0
# supports compressed files (zip, gzip, bgzip)
# requires BioPerl and bgzip to save compressed
# v1.1 further clean after delimiter
#
# visit our Git: https://github.com/Nucleomics-VIB

use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Std;

my $usage="## Usage: fastaCleanHeader.pl <-i fasta_file (required)>
# <-o output file name (default to \"cleaned_\"<infile>; optional)>
# <-c keep only the leftmost word (display_id field; optional)>
# <-s further delete after separator (default to <space>); optional>
# <-d output delimiter (default to \'|\'; optional)>
# <-z to compress results with bgzip>
# <-h to display this help>";

####################
# declare variables
####################
getopts('i:o:d:s:czh');
our ($opt_i, $opt_o, $opt_s, $opt_d, $opt_c, $opt_z, $opt_h);

my $infile = $opt_i || die $usage."\n";
(my $base = $infile) =~ s/\.f[nasta]+(.*z.*)?$//;
my $outfile = $opt_o || $base."_cleaned";
my $indelim = $opt_s || undef;
my $outdelim = $opt_d || "|";
my $clean = $opt_c || undef;
my $zipit = defined($opt_z) || undef;
defined($opt_h) && die $usage."\n";

# load keys from keyfile
my @header = "";

# process multifasta
my $seq_in = OpenArchiveFile($infile);
my $seq_out;
if ( defined($zipit) ) {
	$seq_out = Bio::SeqIO -> new( -format => 'Fasta', -file => "|bgzip -c >$outfile.fa.gz" );
} else {
	$seq_out = Bio::SeqIO -> new( -format => 'Fasta', -file => ">$outfile.fa" );
}

while ( my $seq = $seq_in->next_seq() ) {
	my $curname;
	if ( defined($clean) ) {
		$curname = $seq->display_id();
	} else {
		$curname = $seq->display_id()." ".$seq->accession_number()." ".$seq->desc();
	}
	
	my $newname = $curname;
	
	# further clean after delimiter
	if (defined($indelim)) {
		$newname =~ s/$indelim.*$//g;
	}
	
	# further clean spaces
	$newname =~ s/\ +/$outdelim/g;

	# echo changes
	print STDERR $curname." => ".$newname."\n";
	
	$seq->display_id($newname);
	$seq->accession_number("");
	$seq->desc("");
	$seq_out->write_seq($seq); 
	}

undef $seq_in;

exit 0;

#### Subs ####
sub OpenArchiveFile {
    my $infile = shift;
    my $FH;
    if ($infile =~ /.fa$|.fasta$|.fna$/i) {
    $FH = Bio::SeqIO -> new(-file => "$infile", -format => 'Fasta');
    }
    elsif ($infile =~ /.fa.bz2$|.fasta.bz2$|.fna.bz2$/i) {
    $FH = Bio::SeqIO -> new(-file => "bgzip -c $infile | ", -format => 'Fasta');
    }
    elsif ($infile =~ /.fa.gz$|.fasta.gz|.fna.gz/i) {
    $FH = Bio::SeqIO -> new(-file => "gzip -cd $infile | ", -format => 'Fasta');
    }
    elsif ($infile =~ /.fa.zip$|.fasta.zip$|.fna.zip$/i) {
    $FH = Bio::SeqIO -> new(-file => "unzip -p $infile | ", -format => 'Fasta');
    } else {
	die ("$!: do not recognise file type $infile");
	# if this happens add, the file type with correct opening proc
    }
    return $FH;
}
