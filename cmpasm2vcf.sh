#!/bin/bash

# cmpasm2vcf.sh
# compare assembly to reference and produce VCF
# requires minimap2, k8, and paftools.js
# https://github.com/lh3/minimap2/issues/109
#
# Stephane Plaisance (VIB-NC+BITS) 2020/06/23; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# check parameters for your system
version="1.0, 2020/06/23"

usage='# Usage: cmpasm2vcf.sh -r <reference assembly> -q <query assembly>
# [optional: -L <minimal alignment length|20000>]
# [optional: -t <threads|4>]
# [optional: -h <this help text>]
# [from: https://github.com/lh3/minimap2/issues/109]
# script version '${version}

while getopts "r:q:L:t:h" opt; do
  case $opt in
    r) ref=${OPTARG} ;;
    q) query=${OPTARG} ;;
    L) opt_minl=${OPTARG} ;;
    t) threads=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# defaults
startts=$(date +%s)

# test if minimal arguments were provided
if [ -z "${ref}" ]
then
   echo "# no reference assembly provided!"
   echo "${usage}"
   exit 1
fi

if [ ! -f "${query}" ]; then
	echo "${query} file not found!"
	exit 1
fi

# required
# check executables present
declare -a arr=( "k8" "paftools.js" "minimap2" )
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || \
    ( echo "# required ${prog} not found in PATH"; exit 1 )
done

k8=$(which k8)
paftools=$(which paftools.js)
minimap2=$(which minimap2)

minl=${opt_minl:-20000}

${minimap2} -c --cs ${ref} ${query} \
| sort -k6,6 -k8,8n \
| ${k8} ${paftools} call \
  -f ${ref} \
  -L ${minl} - \
  > ${query%.fa*}_vs_${ref%.fa*}.vcf \
2> ${query%.fa*}_vs_${ref%.fa*}.log

endts=$(date +%s)
dur=$(echo "${endts}-${startts}" | bc)
echo "Done in ${dur} sec"
