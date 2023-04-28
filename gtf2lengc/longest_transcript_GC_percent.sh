#!/bin/bash

# script: longest_transcript_GC_percent.sh
# merges all exons per gene based on a GTF file
# extracts the exon sequences from the reference fasta
# computes the length and GC% for each gene model
# Stephane Plaisance - VIB-Nucleomics Core - 2019-12-23 v1.0

# requires:
# GTF_to_exons.R custom script (SP@NC)
# bedtools
# gawk

# edit below to point to your own files
reffa=RefBeet-1.2.2.fa
gtf=Beta_vulgaris.RefBeet-1.2.2.56.chr.gtf

# end of edits ##############################################

# extract merged exons using R

GTF_to_exons.R -f ${reffa}.fai -g ${gtf}

# extract exon sequences using bedtools
bedtools getfasta \
 -fi ${reffa} \
 -bed merged_exons.bed \
 -fo merged_exons_seqs.txt \
 -tab \
 -nameOnly

# proceed each sequence
mkdir sequences
gawk 'BEGIN{FS="\t"; OFS="\t"}
  {outfile="sequences/"$1".seq"; print $2 >> outfile}' merged_exons_seqs.txt

# merge exons sequences per gene using awk
for seq in sequences/*.seq; do
pfx=$(basename ${seq%.seq});
cat ${seq} | tr -d '\n' \
  | awk -v pfx="${pfx}" '{len=length($0); gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");}
    END{ printf "%s\t%d\t%.2f\n", pfx, len, (gc*100)/(gc+at) }' \
    >> len_gc_percent.txt;
done
