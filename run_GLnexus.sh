#!/bin/bash

# run_GLnexus on all "*.g.vcf.gz" files in the input folder
#
# Stephane Plaisance (VIB-NC) 2023-08-17; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# version number
version="1.1, 2024_03_01"

usage='# Usage: $0 
  -d <gVCF folder> 
  -o <outfolder>
  -c <config (default to DeepVariantWGS>
# script version '${version}'
# [optional: -b <BED file to restrict to chr or region>
# [optional: -t <max threads|4>]
# [optional: -m <max memory|16>]'

while getopts "d:o:c:b:t:m:h" opt; do
  case $opt in
    d) opt_gvcfdir=${OPTARG} ;;
    o) opt_outdir=${OPTARG} ;;
    c) opt_config=${OPTARG} ;;
    b) opt_bed=${OPTARG} ;;
    t) opt_thr=${OPTARG} ;;
    m) opt_mem=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2;
       echo "${usage}" >&2; 
       exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; 
       echo "${usage}" >&2; 
       exit 1 ;;
  esac
done

# GLnexus docker
IMAGE="ghcr.io/dnanexus-rnd/glnexus"
BIN_VERSION="v1.4.3"
img="${IMAGE}:${BIN_VERSION}"

# check for the docker image
res=$(docker images | grep ${img} | cut -d " " -f 1)
if [ ! ${res} == "${IMAGE}" ]; then
  echo "docker image not found"
  exit 1
fi

# check minimal arguments
if [ -z "${opt_gvcfdir}" ]; then
   echo "# no gVCF input folder provided!"
   echo "${usage}"
   exit 1
fi

if [ -z "${opt_outdir}" ]; then
   echo "# no outfolder name provided!"
   echo "${usage}"
   exit 1
fi

# arguments with default
config=${opt_config:-"DeepVariantWGS"}
thr=${opt_thr:-4}
ram=${opt_mem:-16}

# set in and out folders
infolder=${opt_gvcfdir}
outfolder=$(realpath ${opt_outdir})

# prepare analysis, the docker image can only access local files
# copy original inputs to the outfolder
mkdir -p ${outfolder}

# copy gvcf files to outfolder
find ${infolder} -type f -name "*.g.vcf.gz*" -exec cp {} ${outfolder}/ \;

# copy bed if defined to outfolder
if [ -n "${opt_bed}" ]; then
  bed=$(basename "${opt_bed}")
  cp ${opt_bed} ${outfolder}/
  bedarg="--bed /data/${bed}"
else
  bedarg=""
fi

###########################
# merge gVCF using glnexus
###########################

# create a space delimited list of inputs
gvcflist=$(find ${outfolder} -type f -name "*.g.vcf.gz" \
           | tr " " "\n" \
           | sort -k 1V,1 \
           | xargs -n1 basename \
           | gawk -v p="/data/" '{printf "%s%s ",p, $1}' \
           | tr -d "\n")

# bcftools reheader -s sample_names.txt -o reordered.vcf merged.vcf
find ${outfolder} -type f -name "*.g.vcf.gz" -exec basename {} \; \
           | sort -k 1V,1 \
           > ${outfolder}/sample_order.txt

echo "# processing ${gvcflist}"

outfile="merged.vcf.gz"

{ time docker run \
  --rm \
  -v "${outfolder}":"/data" \
  ${img} \
  /usr/local/bin/glnexus_cli \
  --config ${config} \
  ${bedarg} \
  --mem-gbytes ${ram} \
  --threads ${thr} \
  ${gvcflist} \
  | bcftools view - \
  | bcftools reheader -s ${outfolder}/sample_order.txt \
  | bgzip -@ 4 -c \
  > ${outfolder}/${outfile} && \
  tabix -p vcf ${outfolder}/${outfile}; } \
  |& tee ${outfolder}/${outfile%.vcf.gz}_glnexus_log.txt

exit 0

# [1] [2023-07-07 06:47:57.294] [GLnexus] [info] glnexus_cli release v1.4.3-0-gcecf42e Sep 20 2021
# [1] [2023-07-07 06:47:57.295] [GLnexus] [info] detected jemalloc 5.2.1-0-gea6b3e973b477b8061e0076bb257dbd7f3faa756
# Usage: /usr/local/bin/glnexus_cli [options] /vcf/file/1 .. /vcf/file/N
# Merge and joint-call input gVCF files, emitting multi-sample BCF on standard output.
# 
# Options:
#   --dir DIR, -d DIR              scratch directory path (mustn't already exist; default: ./GLnexus.DB)
#   --config X, -c X               configuration preset name or .yml filename (default: gatk)
#   --bed FILE, -b FILE            three-column BED file with ranges to analyze (if neither --range nor --bed: use full length of all contigs)
#   --list, -l                     expect given files to contain lists of gVCF filenames, one per line
# 
#   --more-PL, -P                  include PL from reference bands and other cases omitted by default
#   --squeeze, -S                  reduce pVCF size by suppressing detail in cells derived from reference bands
#   --trim-uncalled-alleles, -a    remove alleles with no output GT calls in postprocessing
# 
#   --mem-gbytes X, -m X           memory budget, in gbytes (default: most of system memory)
#   --threads X, -t X              thread budget (default: all hardware threads)
# 
#   --help, -h                     print this help message#
