## Produce data required for RNASeq normalisation

script: longest_transcript_GC_percent.sh

# requires:
# GTF_to_exons.R custom script (SP@NC)
# bedtools
# gawk

* takes annotations and genome fasta downloaded from ensembl
* isolates all exons and merge them per gene into a exonic gene model to get the full exon footprint (not a real transcript)
* measures the sequence and GC% of each merged sequence
* prints the results in a text file
  geneID, length, GC%

edit the input files and run the shell script **longest_transcript_GC_percent.sh** will call the R script **GTF_to_exons.R** (found in the local folder or in your PATH)

example run:

SP@NC 2023-04-27 v1.0
