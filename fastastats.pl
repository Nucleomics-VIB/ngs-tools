#!/usr/bin/perl -w

# fastastats.pl
# initial version: 2017-05-06
# parse a multi-Fasta file and compute simple stats

# Stephane Plaisance (VIB-NC+BITS) 2017/05/05; v1.0
# added NXX/tot metrics (2017/05/08)
# visit our Git: https://github.com/BITS-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;
use Bio::SeqIO;
use Statistics::Descriptive;
use POSIX qw(strftime);

my $version = 1.1;
my $date = strftime "%m/%d/%Y", localtime;

############################
# handle command parameters
############################

getopts('i:m:x:p:h');
our ( $opt_i, $opt_m, $opt_x, $opt_h );

my $usage = "You must provide a Fasta file with -i
## Usage: fastastats.pl <-i fasta-file (can be compressed)>
# Additional optional filtering parameters are:
# <-m minsize in kb (0)>
# <-x maxsize in kb (1E+6 | 1Gb)>
# <-h to display this help>";

####################
# declare variables
####################

my $infile = $opt_i || die $usage . "\n";
my $minlen = $opt_m || 0;
my $maxlen = $opt_x || 1E+6;
defined($opt_h) && die $usage . "\n";

# handle IO
my $inpath = dirname($infile);
my @sufx = (
	".fa", ".fasta", ".fsa", ".fna",
	".fa.gz", ".fasta.gz", ".fsa.gz", ".fna.gz",
	".fa.zip", ".fasta.zip", ".fsa.zip", ".fna.zip",
	".fa.bz2", ".fasta.bz2", ".fsa.bz2", ".fna.bz2", 
	);
my $outbase = basename( $infile, @sufx );
my $outfile = $inpath . "/" . $outbase . "_stats.txt";

# open stream 
my $in = OpenArchiveFile($infile);
open OUT, "> $outfile" || die $!;

# working variables
my $cnttot   = 0;
my $cntflt   = 0;
my $lentot   = 0;
my $lenflt   = 0;
my $ntot     = 0;
my $nflt     = 0;

# summary variables
my @BigArray = ();
my $stat     = Statistics::Descriptive::Full->new();
my @result   = ();

################################
# parse data and store in array
################################

while ( my $seq = $in->next_seq() ) {
	my $title = $seq->id;
	my $seqlen = $seq->length;
	my $countN = $seq->seq =~ tr/Nn/Nn/;
	push( @BigArray, ( [ $seqlen, $countN ] ) );
}

###################################
# Create filtered molecule subsets
###################################

my @filtered = grep { $_->[0] >= $minlen*1000 && $_->[1] <= $maxlen*1000 } @BigArray;
my $filteredpc = sprintf( "%.1f", 100 * scalar(@filtered) / scalar(@BigArray) );

##################
# prepare summary
##################

my @sizeall = map $_->[0], @BigArray;
# total length in data
$lentot += $_ for @sizeall;
my $lentotf = reverse join ",", (reverse $lentot) =~ /(\d{1,3})/g;

my @sizedistall = map { sprintf( "%d", $_ ) } ( get_stats(@sizeall) );

my ( $n50_length_all, $n50_rank_all ) = get_NXX(50, @sizeall);
my ( $n90_length_all, $n90_rank_all ) = get_NXX(90, @sizeall);

my @nbaseall = map $_->[1], @BigArray;
$ntot += $_ for @nbaseall;
$ntot = reverse join ",", (reverse $ntot) =~ /(\d{1,3})/g;

my @sizeflt = map $_->[0], @filtered;
# total length in filtered data
$lenflt += $_ for @sizeflt;
my $lenfltf = reverse join ",", (reverse $lenflt) =~ /(\d{1,3})/g;
my @sizedistflt = map { sprintf( "%d", $_ ) } ( get_stats(@sizeflt) );

my ( $n50_length_flt, $n50_rank_flt ) = get_NXX(50, @sizeflt);
my ( $n90_length_flt, $n90_rank_flt ) = get_NXX(90, @sizeflt);

my @nbaseflt = map $_->[1], @filtered;
$nflt += $_ for @nbaseflt;
$nflt = reverse join ",", (reverse $nflt) =~ /(\d{1,3})/g;

my $topline = "# fasta length statistics report (v" . $version . "), " . $date;
push( @result, $topline );
push( @result, "# input file: " . basename($infile) );

# print stats for all molecules
push( @result, "\n## Stats for All Molecules");
push( @result, "#-------------------------");
push( @result, "# Molecule count: " . scalar(@BigArray) );
push( @result, "# Total-bases: " . $lentotf );
push( @result, "# N-bases: " . $ntot );

# molecule lengths
push( @result,
      "# molecules-length-N50 [kb]: "
        . ( sprintf( "%.3f", $n50_length_all / 1000 ) ." (".$n50_rank_all." records)" ) );
push( @result,
      "# N50/genome length: " . ( sprintf( "%.3f", 1000 * $n50_length_all / $lentot ) ) );
push( @result,
      "# molecules-length-N90 [kb]: "
        . ( sprintf( "%.3f", $n90_length_all / 1000 ) ." (".$n90_rank_all.")" ) );
push( @result,
      "# N90/genome length: " . ( sprintf( "%.3f", 1000 * $n90_length_all / $lentot ) ) );        
push( @result,
      "\n# quantile distribution for All molecules (" . scalar(@BigArray) . ")" );

my $header = "variable\t"
  . join( "\t",
          ( "min", "25%", "median", "75%", "90%", "max", "mean") );
push( @result, $header );

my $length = "length [kb]\t" . join( "\t", @sizedistall );
push( @result, $length );

####################################
# print stats for filtered molecules

push( @result, "\n## Stats for filtered Molecules");
push( @result, "#------------------------------");
push( @result, "# minimum length [kb]: " . $minlen );
push( @result, "# maximum length [kb]: " . $maxlen );
push( @result, "# Molecule count: " . scalar(@filtered) );
push( @result, "# Filtered-bases: " . $lenfltf );
push( @result, "# N-bases: " . $nflt );

# molecule lengths
push( @result,
      "# molecules-length-N50 [kb]: "
        . ( sprintf( "%.3f", $n50_length_flt / 1000 ) ." (".$n50_rank_flt." records)" ) );
push( @result,
      "# N50/genome length: " . ( sprintf( "%.3f", 1000 * $n50_length_flt / $lenflt ) ) );
push( @result,
      "# molecules-length-N90 [kb]: "
        . ( sprintf( "%.3f", $n90_length_flt / 1000 ) ." (".$n90_rank_flt." records)" ) );
push( @result,
      "# N90/genome length: " . ( sprintf( "%.3f", 1000 * $n90_length_flt / $lenflt ) ) );        
push( @result,
      "\n# quantile distribution for Filtered molecules (" . scalar(@filtered) . ")" );

push( @result, $header );

my $lengthflt = "length [kb]\t" . join( "\t", @sizedistflt );
push( @result, $lengthflt );

# print results to STDOUT and to file
print STDOUT join( "\n", @result ) . "\n";
# echo to file
print OUT join( "\n", @result ) . "\n";
close OUT;

exit 0;

##################
sub OpenArchiveFile {
    my $infile = shift;
    my $FH;
    if ($infile =~ /.fa$|.fasta$|.fna$/) {
    $FH = Bio::SeqIO -> new(-file => "$infile", -format => 'Fasta');
    }
    elsif ($infile =~ /.gz$/) {
    $FH = Bio::SeqIO -> new(-file => "bgzip -cd $infile| ", -format => 'Fasta');
    }
    elsif ($infile =~ /.bz2$/) {
    $FH = Bio::SeqIO -> new(-file => "bzip2 -c $infile| ", -format => 'Fasta');
    }
    elsif ($infile =~ /.zip$/) {
    $FH = Bio::SeqIO -> new(-file => "unzip -p $infile| ", -format => 'Fasta');
    } else {
	die ("$!: do not recognise file type $infile");
	# if this happens, add the file type with correct opening proc
    }
    return $FH;
}

sub get_N50 {
    my @sort = sort { $b <=> $a } @_;
    my $totsum;
    map { $totsum += $_ } @_;
    my $cumsum = 0;     # cumulative sum
    my $ranknum = 0;    # sorted rank
    foreach my $curlength (@sort) {
        $cumsum += $curlength;
        $ranknum++;
        if ( $cumsum >= $totsum / 2 ) {
            return ($curlength, $ranknum);
            last;
        }
    }
}

sub get_NXX {
	my ( $lim, @data ) = @_;
	my $xx = 100/$lim; # 50 => 2
    my @sort = sort { $b <=> $a } @data;
    my $totsum;
    map { $totsum += $_ } @_;
    my $cumsum = 0;     # cumulative sum
    my $ranknum = 0;    # sorted rank
    foreach my $curlength (@sort) {
        $cumsum += $curlength;
        $ranknum++;
        if ( $cumsum >= $totsum / $xx ) {
            return ($curlength, $ranknum);
            last;
        }
    }
}

sub get_stats {

    # uses https://metacpan.org/module/Statistics::Descriptive
    # test input
    return undef unless ( scalar(@_) );

    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@_);
    my @result = (
                   $stat->quantile(0),
                   $stat->quantile(1),
                   $stat->quantile(2),
                   $stat->quantile(3),
                   scalar( $stat->percentile(90) ),
                   $stat->quantile(4),
                   $stat->mean()
                 );
    return @result;
}

