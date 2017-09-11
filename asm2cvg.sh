#!/bin/bash

# usage: asm2cvg.sh -r <reads> -a <assembly> -w <bin size (1000)>
# align reads to a provided fasta assembly with BWA mem
# computer coverage in bins using bedtools
# plot coverage for each chromosome using R
#
# Stephane Plaisance - VIB-Nucleomics Core - September-11-2017 v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# requirements
# functional bedtools 2.25+ in your PATH
# R installed with optparse and ggplot2
# a functional Rscript in your PATH
# the plotting script btcvg2plots.R in your PATH

version="1.0, 2017_09_11"

usage='# Usage: asm2cvg.sh -r <reads> -a <assembly>
# -x <pacbio or ont2d (preset for long reads)>
# -t <threads to be used for mapping>
# -w < window width for coverage stats>
# script version '${version}'
# [optional: -w <window size|1000>]'

while getopts "r:a:x:t:w:h" opt; do
  case $opt in
    r) reads=${OPTARG} ;;
    a) assembly=${OPTARG} ;;
    x) preset=$OPTARG ;;
    t) thread_opt=${OPTARG} ;;
    w) window_opt=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# test if minimal arguments were provided
if [ -z "${reads}" ]; then
   echo "# no read file provided!"
   echo "${usage}"
   exit 1
fi

if [ ! -f "${reads}" ]; then
	echo "${reads} file not found!"
	exit 1
fi

if [ -z "${assembly}" ]; then
   echo "# no assembly file provided!"
   echo "${usage}"
   exit 1
fi

if [ ! -f "${assembly}" ]; then
	echo "${assembly} file not found!"
	exit 1
fi

# bwa mapping preset pacbio or nanopore
if [ "${preset}" == "pacbio" ] || [ "${preset}" == "ont2d" ] ; then
bwapreset="-x ${preset}"
else
bwapreset=""
fi

# defaults
binwidth=${window_opt:-1000}
threads=${thread_opt:-8}

ref_index=${assembly%%.f*}

# create genome size list from .fai indices
samtools faidx ${assembly}
cut -f1,2 ${assembly}.fai > ${assembly}.sizes

# create title lists for plots from fasta headers
grep "^>" ${assembly} > ${assembly}.titles

# create bwa index if does not exists
if [ ! -f "${ref_index}.bwt" ]; then
cmd="bwa index -p ${ref_index} ${ref_index}.fasta"

echo "# ${cmd}"
eval ${cmd}
fi

ref_mappings="${reads%%.fastq*}_to_${assembly%%.f*}-alignments.bam"

# if not done yet, map reads to index, sort results, and index bam
if [ ! -f "${ref_mappings}" ]; then
cmd="bwa mem ${bwapreset} -t ${threads} ${ref_index} ${reads} \
	| samtools view -Sb - \
	| samtools sort - \
	> ${ref_mappings} && samtools index ${ref_mappings}"
	
echo "# ${cmd}"
eval ${cmd}
fi

# create windows
if [ ! -f "${assembly}.${binwidth}bp-bin.bed" ]; then
cmd="bedtools makewindows -g ${assembly}.sizes -w ${binwidth} \
	> ${assembly}.${binwidth}bp-bin.bed"

echo "# ${cmd}"
eval ${cmd}
fi

# bedtools genome wide coverage analysis
if [ ! -f "bedtools_coverage.${binwidth}bp_${assembly}-hist.txt" ]; then
cmd="bedtools coverage -b ${ref_mappings}\
  -a ${assembly}.${binwidth}bp-bin.bed \
  -hist \
  > bedtools_coverage.${binwidth}bp_${assembly}-hist.txt"

echo "# ${cmd}"
eval ${cmd}
fi

# merge and generate stats for each bin (ignore 'all' rows)
if [ ! -f "bedtools_coverage.${binwidth}bp_${assembly}-stats.txt" ]; then
cmd="bedtools groupby \
  -i <(grep -v "^all" bedtools_coverage.${binwidth}bp_${assembly}-hist.txt) \
  -g 1,2,3 \
  -c 5 -o min,median,mean,max \
  > bedtools_coverage.${binwidth}bp_${assembly}-stats.txt"

echo "# ${cmd}"
eval ${cmd}
fi

# run R plotting script
cmd="btcvg2plots.R -b bedtools_coverage.${binwidth}bp_${assembly}-stats.txt \
	-t ${assembly}.titles"
	
echo "# ${cmd}"
eval ${cmd}

