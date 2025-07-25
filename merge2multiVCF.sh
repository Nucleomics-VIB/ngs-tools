#!/bin/bash

# script: merge2multiVCF.sh
# merge multiple single genome & single chromosome VCF files
# to a single multi-sample VCF file
# requirements: same reference used for all files and vcf indexed
#
# Stephane Plaisance (VIB-NC) 2023/09/28; v1.0
# visit our Git: https://github.com/Nucleomics-VIB

version="1.0"

usage="# Usage: merge2multiVCF.sh -i <folder with vcf files> -c <chromosome>
# script version '${version}'
# [optional: -o <output file|default to ${chr}.vcf.gz>]
# [optional: -n <threads|1>]"

while getopts ":i:c:o:n:h" opt; do
  case $opt in
    i)
      infolder="$OPTARG"
      ;;
    c)
      chr="$OPTARG"
      ;;
    o)
      output_vcf="$OPTARG"
      ;;
    n)
      nthr="$OPTARG"
      ;;
    h) 
      echo "${usage}" >&2; exit 0
      ;;  
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Check if the required arguments are provided
if [[ -z "$infolder" || -z "$chr" ]]; then
  echo "${usage}"
  exit 1
fi

# List of individual VCF files for your samples
vcf_files=( $(ls "${infolder}"/*_"${chr}"_*.vcf.gz) )

output_vcf=${output_vcf:-"${chr}.vcf.gz"}

nthr=${nthr:-1}

# Use bcftools to merge the VCF files without merging variants
bcftools merge \
  --threads "${nthr}" \
  -m none \
  -o "${output_vcf}" \
  "${vcf_files[@]}" && \
    bcftools index \
      -t "${output_vcf}" \
      --threads "${nthr}"
