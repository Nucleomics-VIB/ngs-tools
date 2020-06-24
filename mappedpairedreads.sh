#!/bin/bash

# mappedpairedreads.sh
# BWA mem commands to collect reads mapping to a reference for denovo asm
# converts the results back to compressed fastQ for de-novo asm

# Stephane Plaisance (VIB-NC) 2020/06/23; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# check parameters for your system
version="1.0, 2020_06_23"

usage='# Usage: mappedlongreads.sh -i <fasta assembly> -1 <PE reads1> -2 <PE reads2>
# [optional: -t <threads|1>]
# [optional: -T <threads2|1>]
# [optional: -h <this help text>]
# script version '${version}

while getopts "i:1:2:t:T:h" opt; do
  case $opt in
    i) ref=${OPTARG} ;;
    1) reads1=${OPTARG} ;;
    2) reads2=${OPTARG} ;;
    t) opt_thr=${OPTARG} ;;
    T) opt_thr2=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

startts=$(date +%s)

# threads
thr=${opt_thr:-1}
thr2=${opt_thr2:-1}

# type can be map-pb or map-ont
# for PB, add -H use homopolymer-compressed k-mer (preferrable for PacBio)
type=${opt_plf:-"map-pb -H"}

if [[ -z "${ref}" ]]; then
  echo "# no reference assembly provided!"
  echo "${usage}"
  exit 1
fi

if [[ ! -f "${ref}" ]]; then
	echo "${ref} file not found!"
	exit 1
fi

if [[ -z "${reads1}" ]] || [[ -z "${reads2}" ]]; then
	echo "# requires two matched PE read files !"
	echo "${usage}"
	exit 1
fi

if [[ ! -f "${reads1}" ]] ; then
    echo "reads file ${reads1} not found!";
    exit 1
fi

if [[ ! -f "${reads2}" ]] ; then
    echo "reads file ${reads2} not found!";
    exit 1
fi


#################################
# BWA mem align to ref assembly #
#################################

mkdir -p bwaidx
bwa index ${ref} -p bwaidx/$(basename ${ref%.fa*})

bwa mem \
  -t ${thr} \
  bwaidx/$(basename ${ref%.fa*}) \
  ${reads1} \
  ${reads2} \
| samtools fastq -@${thr2} \
  -N \
  -F 4 \
  -1 >(bgzip -c > $(basename ${ref%.fa*})_mapped_$(basename ${reads1%.f*}).fq.gz) \
  -2 >(bgzip -c > $(basename ${ref%.fa*})_mapped_$(basename ${reads2%.f*}).fq.gz) \
  -0 /dev/null \
  -s /dev/null \
  - \
  && rm -rf bwaidx

endts=$(date +%s)
dur=$(echo "${endts}-${startts}" | bc)
echo "Done in ${dur} sec"