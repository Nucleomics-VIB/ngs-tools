#!/usr/bin/env perl -w

# http://omicsomics.blogspot.com/2013/12/assembly-could-benefit-from-more.html
# circ-perm.pl split=200223 max=4561237 in=orig.fa out=perm.fa 

use strict;  
use Bio::SeqIO;  
use Getopt::Long;  

my $splitPoint=undef;  
my $inFile=undef;  
my $reverse=0;  
my $outFile=undef;  
my $prefix="";  
my $maxPos=undef;  
my $regexp=".";  
my $targetId=undef;  

&GetOptions(
  'split=i'=>\$splitPoint,
  'in=s'=>\$inFile,
  'out=s'=>\$outFile,
  'rev'=>\$reverse,
  'reg=s'=>\$regexp,
  'pre=s'=>\$prefix,
  'max=i'=>\$maxPos,
  'id=s'=>\$targetId
  );  
 
$inFile=shift if (!defined $inFile && scalar(@ARGV)==1);  

unless ( defined $splitPoint && defined $inFile && defined $outFile )  
 {  
   die "circ-perm.pl -in inFile -out outFile -split spsplitPoint (-r)\n";  
 }  

my $rdr=new Bio::SeqIO(-file=>$inFile);  
my $writer=new Bio::SeqIO(-format=>'fasta',-file=>">$outFile");  

while (my $rec=$rdr->next_seq)  
 {  
   if ($rec->id=~/$regexp/ || (defined $targetId && $rec->id eq $targetId))  
   {  

     my $headSeq=$rec->subseq(1,$splitPoint);  
     $maxPos=$rec->length unless (defined $maxPos);  
     my $tailSeq=$rec->subseq($splitPoint+1,$maxPos);  
     my $revFlag=""; $revFlag=".rc" if ($reverse);  
     my $newSeq=new Bio::Seq(
       -id=>$prefix.$rec->id.".cp.$splitPoint$revFlag",
       -seq=>$tailSeq.$headSeq
       );  
     $newSeq=$newSeq->revcom if ($reverse);  
     $writer->write_seq($newSeq);  
   }  
   else  
   {  
     $writer->write_seq($rec);  
   }  
 }
