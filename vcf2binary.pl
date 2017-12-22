#!/usr/bin/perl -w

# vcf2binary.pl
# convert multisample VCF data to binary tsv format for counting and clustering (GT field)
# !! NA genotypes '.' are counted as 0 (ref)
#
# Stephane Plaisance (VIB-NC) 2017/13/19; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;
use POSIX qw(strftime);
use List::Util qw( sum );

my $version = 1.0;
my $date = strftime "%m/%d/%Y", localtime;

# autoflush
$|=1;

#############################
# handle command parameters #
#############################

getopts('i:m:M:c:bh');
our ( $opt_i, $opt_m, $opt_M, $opt_b, $opt_c, $opt_h );

my $usage = "## Usage: vcf2binary.pl <-i vcf-file (can be compressed)>
# <-m x (minimal number of variant samples, default to 1)>
# <-M y (minimal number of variant samples, default to sample#)>
# optional <-b write data as one title column (chr_start_ref_alt) + binary matrix>
# ! if -b is unset, writes the first 8 VCF columns + variant count + matrix>
# optional <-c prefix for chr (defrault to '')>
# <-h to display this help>
version: ".$version;

my $vcf_file = $opt_i || die $usage . "\n";
my $min_var = $opt_m || 1;
my $binmat = $opt_b || 0;
my $pfx = $opt_c || '';
defined($opt_h) && die $usage . "\n";

####################
# declare variables
####################

# make outname
my @suffixlist=(".vcf",".vcf.gz");
my $name = basename($vcf_file, @suffixlist);
my $outfile;

# collect sample names from last header line
my @samples = ();

# limits and counters
my $max_var = 0;
my $cntln = 0;
my $found_header = 0;
my $counter = 0;
my $sample_cnt = 0;
my $header;

##########################################
# read first line and check if valid VCF #
##########################################

my $VCF = OpenArchiveFile($vcf_file) or die $!;
my $firstline=<$VCF>;
$firstline =~ /##fileformat=VCFv/ || die "##fileformat=VCFv not found in first row, check this is a VCF file";

##################
# parse VCF file #
##################

while ( my $line = <$VCF> ) {

	# check header
	$cntln++;
	($cntln>1000 && $found_header==0) && die "# header row not found in first 1'000 lines!";
	next if ($line !~ /^#CHROM/ && $found_header==0);
	
	# one-time create output and print header
	if ( $line =~ /^#CHROM/ ) {
		$found_header=1;
    	chomp($line);
    	@samples = split("\t", $line);
    	@samples = @samples[9..$#samples];
    
		# check valid VCF
		$sample_cnt = scalar(@samples);
		$sample_cnt == 0 && die "No sample column found, check you provided a valid VCF file!";
		$max_var = $opt_M || $sample_cnt;

		# echo 
		print STDOUT ($#samples+1)." samples found in VCF data\n";
		print STDOUT "# filtering variants in at least $min_var, and at most $max_var samples\n";

		if ($binmat == 0) {
			$outfile = $name."_ge".$min_var."_le".$max_var."_full_binary.txt";
			$header= join("\t", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "VARCNT_SMPL", @samples);
		} else {
			$outfile = $name."_ge".$min_var."_le".$max_var."_binary.txt";
			$header= join("\t", ($pfx.join("_", "CHROM", "POS", "REF", "ALT", "VARCNT_SMPL")), @samples);
		}

		# create output file
		open (OUT, " > $outfile") or die "Error: cannot write to $outfile: $!\n";
			
		# print header
		print OUT $header."\n";
		#print STDERR $header."\n";
	}

	# count data lines and echo every 10000
	$counter++;
	$counter =~/0000$/ && print STDOUT ".";
	 
	# ignore header rows
	next if $line =~ /^#/;
	my @col = split("\t", $line);
	my @genotypes = @col[9..$#col];
	my @format = $col[8];
	my ( $gtidx )= grep { $format[$_] =~ /GT/ } 0..$#format;

	# convert genotypes to binary
	my @binaries = ();
	my $sum_all = 0;
	
	foreach my $sample (@genotypes) {
		my $gt = (split(":", $sample))[$gtidx];
		# remove all non numeric characters './|'
		$gt =~ s/\D/0/g;
		my $sum = sum (split(//, $gt));
		my $bin = ( $sum > 0 ) ? 1 : 0;
		$sum_all += $bin;
		push @binaries, $bin;
		}

	my @info = ();
	# print only when at least $min_var genome(s) is/are variant
	if ($sum_all >= $min_var && $sum_all <= $max_var) {
		if ($binmat == 0) {
			@info = @col[0..6];
			push @info, $sum_all."_".$sample_cnt;
		} else {
			my $rowname = $pfx.(join("_",  @col[1,2,4,5], $sum_all, $sample_cnt));
			push @info, $rowname;
		}	

		print OUT join("\t", @info, @binaries)."\n"; 
		#print STDERR join("\t", @info, @binaries)."\n";
		}
		
	# EOF
	}
	
undef $VCF;
close OUT;
exit 0;

##############
#### Subs ####

sub OpenArchiveFile {
	# $Filename passed in, handle to file returned
	my $File = shift;    # filename
	my $FH;              # file handle

	if ( $File =~ /.vcf$/ ) {
		open( $FH, "cat $File | " )
		  or die("$!: can't open file $File");
	} elsif ( $File =~ /.zip$/ ) {
		open( $FH, "unzip -p $File | " )
		  or die("$!: can't open file $File");
	} elsif ( $File =~ /(.gzip|.gz)$/ ) {
		open( $FH, "gzip -dc $File | " )
		  or die("$!: can't open file $File");
	} else {
		die("$!: the file $File does seem to be a 'vcf' file");
	}
	return $FH;
}
