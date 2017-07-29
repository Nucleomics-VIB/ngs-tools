#!/bin/bash

# script vcf2index.sh
# sort and compress vcf data and add tabix index
#
# Stephane Plaisance (VIB-NC) 2017/07/29; v1.0

usage="Usage: vcf2index.sh -i <vcf.file>"

while getopts "i:" opt; do
  case $opt in
    i)
      infile=${OPTARG}
      ;;
    \?)
      echo ${usage}
      exit 1
      ;;
    *)
      echo ${usage} >&2
      exit 1
      ;;
  esac
done

# test if minimal arguments were provided
if [ -z "${infile}" ]
then
   echo "# no input provided!"
   echo "${usage}"
   exit 1
fi

# filter string
filter="-k 1V,1 -k 2n,2"

( grep ^"#" ${infile}; grep -v ^"#" ${infile} | sort ${filter} ) | \
	bgzip -c > ${infile}".gz" && \
	tabix -f -p vcf ${infile}".gz"