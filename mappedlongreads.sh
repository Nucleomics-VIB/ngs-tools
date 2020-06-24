#!/bin/bash

# mappedlongreads.sh
# minimap2 commands to collect reads mapping to a reference for denovo asm
# converts the results back to compressed fastQ for de-novo asm

# Stephane Plaisance (VIB-NC) 2020/06/23; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# check parameters for your system
version="1.0, 2020_06_23"

usage='# Usage: mappedlongreads.sh -i <fasta assembly> -r <long reads>
# [optional: -p <platform map-pb|map-ont (default to "map-pb -H")>]
# [optional: -t <threads|1>]
# [optional: -T <threads2|1>]
# [optional: -h <this help text>]
# script version '${version}

while getopts "i:r:p:t:T:h" opt; do
  case $opt in
    i) ref=${OPTARG} ;;
    r) reads=${OPTARG} ;;
    p) opt_plf=${OPTARG} ;;
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

if [[ -z "${reads}" ]]; then
	echo "# requires long reads !"
	echo "${usage}"
	exit 1
fi

if [[ ! -f "${reads}" ]] ; then
    echo "reads file ${reads} not found!";
    exit 1
fi

##################################
# map long reads to raw assembly #
##################################

minimap2 \
  -ax ${type} \
  -t ${thr} \
  ${ref} \
  ${pbreads} \
| samtools view -h -F 4 -O BAM - \
| samtools fastq -@${thr2} - \
| bgzip -c > $(basename ${ref%.fa*})_mapped_$(basename ${pbreads%.fq*}).fq.gz

endts=$(date +%s)
dur=$(echo "${endts}-${startts}" | bc)
echo "Done in ${dur} sec"