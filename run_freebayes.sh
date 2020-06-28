#!/bin/bash

# call variants using freebayes at CLI
#
# Requirements:
# run on a unix computer installed with freebayes
# reference genome fasta present and indexed
#
# this script can be used for a batch of BNX files like:
# for b in *.bam; do
#	run_freebayes.sh -i $b -r reference.fa;
# done
#
# Stephane Plaisance (VIB-NC+BITS) 2015/09/28; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# check parameters for your system
TOOLS=$BIOTOOLS
version="1.0, 2016_09_28"

usage='# Usage: run_freebayes.sh -i <bam file> -r <fasta reference>
# script version '${version}'
# [optional: -m <minmapq|20>]
# [optional: -q <minbaseq|20>]
# [optional: -F <minaltfrac|0.01>]
# [optional: -C <minaltcnt|10>]
# [optional: -n <add parameters between quotes>]'

while getopts "i:r:m:q:F:C:n:h" opt; do
  case $opt in
    i) bamfile=${OPTARG} ;;
    r) reffile=${OPTARG} ;;
    m) minmapqual=${OPTARG} ;;
    q) minbasequal=${OPTARG} ;;
    F) minaltfraction=${OPTARG} ;;
    C) minaltcount=${OPTARG} ;;
    n) moreoptions=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# defaults
startts=$(date +%s)

# test if minimal arguments were provided
if [ -z "${bamfile}" ]
then
   echo "# no bam provided!"
   echo "${usage}"
   exit 1
fi

if [ ! -f "${bamfile}" ]; then
	echo "${bamfile} file not found!"
	exit 1
fi

# define outfile
outfile="${bamfile%%.bam$}.vcf"

if [ -z "${reffile}" ]
then
	echo "# no reference cmap provided"
	echo "${usage}"
	exit 1
fi

if [ ! -f "${reffile}" ]; then
    echo "${reffile} file not found!";
    exit 1
fi

# other parameters or defaults
minmapq=${minmapqual:-20}
minbaseq=${minbasequal:-20}
minaltfrac=${minaltfraction:-0.01}
minaltcnt=${minaltcount:-10}

if [ -z "${moreoptions}" ]
then
	moreopt="${moreoptions}"
else
	moreopt=''
fi

# build the command
cmd="freebayes \
	-f ${reffile} \
	--min-mapping-quality ${minmapq} \
	--min-base-quality ${minbaseq} \
	--min-alternate-fraction ${minaltfrac} \
	--min-alternate-count ${minaltcnt} \
	${moreopt} \
	${bamfile} > ${outfile}"

# show and execute
echo "# ${cmd}"
eval ${cmd}

endts=$(date +%s)
dur=$(echo "${endts}-${startts}" | bc)
echo "Done in ${dur} sec"
