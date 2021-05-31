#!/usr/bin/perl -w

# addflancs2circ.pl -in=<input_file.fa> -out=<output_file.fa> -addflank=1000 
# in order to better map reads crossing the sequence origin
# take a circular sequence (plasmid fasta)
# Using the current sequence
#   write addFlank (default=1000) last bases at front
#   write current sequence
#   write addFlank (default=1000) first bases at back
# if flank >= total length take only half length as flank

use strict;
use Bio::SeqIO;
use Getopt::Long;

# variables
my $defaultFlank=1000;
my $maxPos=undef;

GetOptions(
  'in=s'=>\my $inFile,
  'out=s'=>\my $outFile,
  'flank=i'=>\(my $addFlank = 1000),
  ) or die "syntax: addflancs2circ.pl -in=<input_file.fa> -out=<output_file.fa> -addflank=<length to add at either side (default 1000bps)>\n";

print "$inFile - $outFile - $addFlank\n";

# connect to input file for reading
my $rdr=new Bio::SeqIO(-file=>$inFile);

# connect for output
my $writer=new Bio::SeqIO(-format=>'fasta',-file=>">$outFile");

# process each input sequence (accepts multifasta)
while (my $rec=$rdr->next_seq)
 {  
     my $fullSeq=$rec->seq();
     $maxPos=$rec->length;
     # flank cannot be larger than seq/2
     $addFlank=$defaultFlank unless (defined $addFlank);
     # if flank >= total length take only half length as flank
     if ( $addFlank >= $maxPos )
      { 
        $addFlank=$maxPos/2
      }
     my $newId=$rec->id.".with-".$addFlank."bps-flanks";
     my $headSeq=$rec->subseq(1,$addFlank);
     my $tailSeq=$rec->subseq($maxPos-$addFlank,$maxPos);
     my $newSeq=new Bio::Seq(-id=>$newId,-seq=>$tailSeq.$fullSeq.$headSeq);
     $writer->write_seq($newSeq);
 }
