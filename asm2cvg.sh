#!/bin/bash

# usage: asm2cvg.sh -r <reads> -a <assembly> -w <bin size (1000)>
# align reads to a provided fasta assembly with BWA mem
# computer coverage in bins using bedtools
# plot coverage for each chromosome using R
# create mosaic using imagemagic
#
# Stephane Plaisance - VIB-Nucleomics Core - September-11-2017 v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# requirements
# functional bedtools 2.25+ in your PATH
# R installed with optparse and ggplot2
# a functional Rscript in your PATH
# the plotting script btcvg2plots.R in your PATH

version="1.2, 2017_10_10"

usage='# Usage: asm2cvg.sh -r <reads> -a <assembly>
# -x <pacbio or ont2d (preset for long reads)>
# -s <small genome, use "-a is" for bwa index (default undef)>
# -m <existing bam data - no mapping will be done (default undef)>
# -t <threads to be used for mapping>
# -w <window width for coverage stats>
# -p <plot stat [min,mean,median,max,all] (default median)]>
# script version '${version}'
# [optional: -w <window size|1000>]'

while getopts "r:a:m:x:t:w:p:sh" opt; do
  case $opt in
    r) reads=${OPTARG} ;;
    a) assembly=${OPTARG} ;;
    m) mappings=${OPTARG} ;;
    s) small=1 ;;
    x) preset=$OPTARG ;;
    t) thread_opt=${OPTARG} ;;
    w) window_opt=${OPTARG} ;;
    p) plot_opt=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# check if executables are present
declare -a arr=("samtools" "bwa" "bedtools" "btcvg2plots.R" "montage")
for prog in "${arr[@]}"; do
eval $( hash "${prog}" 2>/dev/null ) || ( echo "# required ${prog} not found in PATH"; exit 1 )
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

# is genome small or large?
if [ ${small}==1 ]; then
indexpreset="-a is"
else
indexpreset="-a bwtsw"
fi

# metric(s) to plot?
if [ -n "${plot_opt}" ]; then
plot=${plot_opt}
else
plot="median"
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

ref_mappings=${mappings:-"$(basename ${reads%.f*})_to_$(basename ${assembly%.f*})-alignments.bam"}

# if not done yet, map reads to index, sort results, and index bam
if [ ! -f "${ref_mappings}" ]; then

# create bwa index if does not exists
if [ ! -f "${ref_index}.bwt" ]; then
cmd="bwa index ${indexpreset} -p ${ref_index} ${ref_index}.fasta"

echo "# ${cmd}"
eval ${cmd}
fi

# map reads to reference index
cmd="bwa mem ${bwapreset} -t ${threads} ${ref_index} ${reads} | \
  samtools sort -o ${ref_mappings} - && \
    samtools index ${ref_mappings}"

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

# run R plotting script for all contigs / chromosomes
cmd="btcvg2plots.R -b bedtools_coverage.${binwidth}bp_${assembly}-stats.txt \
	-t ${assembly}.titles \
	-s ${plot}"
	
echo "# ${cmd}"
eval ${cmd}

# create mosaic
montage -mode concatenate coverage_plots-${assembly%%.*}/*.pdf cvg_mosaic-${assembly}.pdf