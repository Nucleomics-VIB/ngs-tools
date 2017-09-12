#!/bin/bash

# nucmer_biplot.sh: produce a pairwise plot from two fasta sequences
#
# Requirements:
# run on a unix computer installed with mummer3 (mummer apps in $PATH)
# two related fasta references to be compared
# this code was designed for DNA
#
# Stephane Plaisance (VIB-NC+BITS) 2017/09/22; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# check parameters for your system
version="1.0, 2017_09_22"

usage='# Usage: nucmer_biplot.sh -x <first assembly> -y <second assembly>
# script version '${version}'
# [optional: -o <result folder>]
# [optional: -c <min-cluster|100>]
# [optional: -h <this help text>]'

while getopts "x:y:o:c:p:h" opt; do
  case $opt in
    x) assembly1=${OPTARG} ;;
    y) assembly2=${OPTARG} ;;
    p) mummerpath=${OPTARG} ;;
    o) outpathopt=${OPTARG} ;;
    c) clust=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# defaults
cluster=${clust:-100}

# test if minimal arguments were provided
if [ -z "${assembly1}" ]
then
   echo "# no first assembly provided!"
   echo "${usage}"
   exit 1
fi

if [ ! -f "${assembly1}" ]; then
	echo "${assembly1} file not found!"
	exit 1
fi

if [ -z "${assembly2}" ]
then
	echo "#  no second assembly provided!"
	echo "${usage}"
	exit 1
fi

if [ ! -f "${assembly2}" ]; then
    echo "${assembly2} file not found!";
    exit 1
fi

# check if mummer requirements are present
$( hash nucmer 2>/dev/null ) || ( echo "# nucmer not found in PATH"; exit 1 )
$( hash show-coords 2>/dev/null ) || ( echo "# show-coords not found in PATH"; exit 1 )
$( hash mummerplot 2>/dev/null ) || ( echo "# mummer-plot not found in PATH"; exit 1 )

# other parameters or defaults
outpath=${outpathopt:-"mummer_results"}
mkdir -p ${outpath}

result="${outpath}/mummer_plot-${assembly1%.f*}_vs_${assembly2%.f*}"

# build the command
cmd="nucmer --maxmatch \
	-c ${cluster} \
	-p ${result} \
	${assembly1} ${assembly2}
	> mummer3-log.txt 2>&1"

# show and execute	
echo "# ${cmd}"
eval ${cmd}
 
# after success create alignment file and plot
if [ $? -eq 0 ]; then
	cmd="(show-coords -r -c -l ${result}.delta > ${result}-all_coords.txt && \
		mummerplot --fat --filter --png --large -p ${result} ${result}.delta) \
		>> mummer3-log.txt 2>&1"
	echo "# ${cmd}"
	eval ${cmd}
else
    echo "Mummer analysis seems to have failed, please check nucmer-log.txt!"
fi
