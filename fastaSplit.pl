#!/usr/bin/perl -w

# fastaSplit.pl input.fa prefix bins
# adapted from https://www.biostars.org/p/2226/
# also available as 'faSplit sequence input.fa 100 part_' (Kent tools)

use strict;
use Bio::SeqIO;

my $from = shift;
my $toprefix = shift;
my $seqnum = shift;

# count records in infile
my $tot = `grep -c "^>" $from`;
chomp($tot);
my $dig = length(int($tot / $seqnum))+1;

# process input
my $in  = new Bio::SeqIO(-file  => $from);

my $count = 0;
my $fcount = 1;
my $suff = sprintf ("%0${dig}d", $fcount );
# create first output
my $out = new Bio::SeqIO(-file => ">".$toprefix."_".$suff.".fasta", -format=>'fasta');

while (my $seq = $in->next_seq) {
        if ($count > 0 && $count % $seqnum == 0) {
                $fcount++;
                $suff = sprintf ("%0${dig}d", $fcount );
                print STDOUT $suff;
                $out = new Bio::SeqIO(-file => ">".$toprefix."_".$suff.".fasta", -format=>'fasta');
        }
        $out->write_seq($seq);
        $count++;
}