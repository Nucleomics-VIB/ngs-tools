#!/usr/bin/perl
use strict;
use File::Basename;
use Switch;
use DateTime qw( );
$| = 1;

# parse an augustus GFF prediction and matching blastP output
# replace augustus arbitrary gene names by blastP best match
# usage: gffRenameGM.pl <augustus.gff> <blastp.txt>
#
# Stéphane Plaisance VIB-BITS May-06-2022 v1.0

my $script = basename($0);
my $datetag = scalar localtime ();
my $asm = "R64-3-1_20210421";

@ARGV == 2 || die ("usage: $script augustus.gff blastp.txt");

my($gff, $blastp) = @ARGV;
chomp($gff);
chomp($blastp);

# load blastp name table in name hash
open BLAST, "<", $blastp or die "Can't read file '$blastp' [$!]\n";
my %names=();

foreach my $row (<BLAST>){
    chomp $row;
    my @cols = split /\t/, $row;
    $names{$cols[0]}=$cols[1];
}

# parse results (to be removed)
for(keys %names){
	print STDERR "name for $_ is $names{$_}\n";
}

# parse GFF and rename GM
open GFF, "<", $gff or die "Can't read file '$gff' [$!]\n";

# print header
print STDOUT "##gff-version 3\n";
print STDOUT "#!date-produced $datetag\n";
print STDOUT "#!data-source AUGUSTUS (3.4.0)\n";
print STDOUT "#!assembly $asm\n";

foreach my $line (<GFF>){

    # print through lines starting with #
    if ($line =~ /^#/) {
    #print STDOUT $line;
    } else {
    my @field = split /\t/, $line;

    # replace name in field #9 (column #8)
    my $type = $field[2];
    my $annot = $field[8];
    my $augid;
	chomp $annot;
	
    switch($type) {
		case "gene" {
			# g3
			#print STDERR "# $type $annot\n";
			$annot =~ s/$annot/$names{$annot}/g;
			my $newline = join( "\t", @field[0..$#field-1], $annot );
			print STDOUT "$newline\n";
			}
		case ["intron","start_codon","stop_codon","CDS"] {
			# transcript_id "g3.t1"; gene_id "g3";
			#print STDERR "# $type $annot\n";
			if ($annot =~ /transcript_id.*gene_id \"(.*)\";/){
				$augid = $1;
				$annot =~ s/$augid/$names{$augid}/g;
				my $newline = join( "\t", @field[0..$#field-1], $annot );
				print STDOUT "$newline\n";
				}
			}
		case "transcript" {
			# g6.t1
			#print STDERR "# $type $annot\n";
			if ($annot =~ /(.*)\..*/){
				$augid = $1;
				$annot =~ s/$augid/$names{$augid}/g;
				my $newline = join( "\t", @field[0..$#field-1], $annot );
				print STDOUT "$newline\n";
				}
			}
		else {
			die "unrecognized annotation type"
			}
		}
    
    }
 
 }
