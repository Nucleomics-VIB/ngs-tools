#!/usr/bin/env Rscript

# script: GTF_to_exons.R
# Aim: create a BED file from a GTF annotation file and genome index
# requires "GenomicFeatures", "GenomicRanges", "rtracklayer", "argparse"
#
# Stéphane Plaisance - VIB-Nucleomics Core - 2023-04-27 v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# Load libraries and check if they are installed
suppressPackageStartupMessages(suppressWarnings(library("GenomicFeatures")))
suppressPackageStartupMessages(suppressWarnings(library("GenomicRanges")))
suppressPackageStartupMessages(suppressWarnings(library("rtracklayer")))
suppressPackageStartupMessages(suppressWarnings(library("argparse")))

# Define command line arguments
parser <- ArgumentParser()
parser$add_argument("-f", "--fai", dest="fai_file", help="Path to fasta index file", required=TRUE)
parser$add_argument("-g", "--gtf", dest="gtf_file", help="Path to GTF annotation file", required=TRUE)
parser$add_argument("-p", "--prefix", dest="prefix", help="Prefix for the output BED files", required=FALSE)
args <- parser$parse_args()

# Check if both arguments were given
if (is.null(args$fai_file) || is.null(args$gtf_file)) {
  stop("Both arguments --fai and --gtf must be given.")
}

# Set prefix to empty string if not provided
if (is.null(args$prefix)) {
  args$prefix <- "merged"
}

############
# the code
############

chrominfo <- read.table(file = args$fai_file, sep = '\t')
chrominfo <- chrominfo[,1:2]
colnames(chrominfo) <- c("chrom", "length")

# import
txdb <- GenomicFeatures::makeTxDbFromGFF(organism = NA, 
                                         format = "gtf",
                                         file =  args$gtf_file, 
                                         chrominfo = chrominfo
                                         )

# get exons and save as bed
exons <- exonsBy(txdb, by = c("gene"))
exons <- unlist(exons)
bed1 <- paste0(args$prefix,"_exons.bed")
rtracklayer::export(exons, bed1)
