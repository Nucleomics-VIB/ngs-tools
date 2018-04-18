#!/bin/bash

# assemblytics.sh: produce a assembly QC based on mummer comparison to a reference
#
# Requirements:
# run on a unix computer
# nucmer (mummer3 apps) and Assemblytics both in $PATH
# two related fasta references to be compared
#
# Stephane Plaisance (VIB-NC+BITS) 2017/11/27; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# check parameters for your system
version="1.0, 2017_11_27"

# path to the Assemblytics scripts
default_path_to_scripts="/opt/biotools/Assemblytics/"

usage='# Usage: assemblytics.sh -x <reference assembly> -y <de-novo assembly>
# script version '${version}'
# [optional: -o <result folder|Assemblytics_results>]
# [optional: -w <uniqseqlen|10000>]
# [optional: -p <path to scripts|default set in the code>]
# [optional: -h <this help text>]'

while getopts "x:y:o:w:p:h" opt; do
  case $opt in
    x) assembly1=${OPTARG} ;;
    y) assembly2=${OPTARG} ;;
    o) outpathopt=${OPTARG} ;;
    w) uniqseqlen=${OPTARG} ;;
    p) p2script=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# defaults
path_to_scripts=${p2script:-${default_path_to_scripts}}
uniq_seq_len=${uniqseqlen:-"10000"}

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

# check if nucmer requirements are present
prog="Assemblytics"
$( hash ${prog} 2>/dev/null ) || ( echo "# ${prog} not found in PATH (nucmer or promer?)"; exit 1 )
$( hash nucmer 2>/dev/null ) || ( echo "# mummer not found in PATH"; exit 1 )

# other parameters or defaults
outpath=${outpathopt:-"Assemblytics_results"}
asm1=$(basename ${assembly1%.f*})
asm2=$(basename ${assembly2%.f*})
deltabase=${asm2}_vs_${asm1}
mkdir -p ${outpath}

# build the nucmer command
nucmercmd="nucmer -maxmatch -l 100 -c 500 \
	${assembly1} \
	${assembly2} \
	-prefix ${outpath}/${deltabase} \
	> ${outpath}/assemblytics-log.txt 2>&1"

# show and execute	
echo "# pairwise alignments with nucmer: ${nucmercmd}"
eval ${nucmercmd}

# check for failure
if [ $? -ne 0 ]; then
	echo "# the nucmer command failed, please check your inputs"
	exit 0
fi

# build the Assemblytics command
cmd="(Assemblytics ${outpath}/${deltabase}.delta \
	${outpath}/${asm2}_contigs \
	${uniq_seq_len} \
	${path_to_scripts}) \
	>> ${outpath}/assemblytics-log.txt 2>&1"

# show and execute	
echo "# Assemblytics analysis: ${cmd}"
eval ${cmd}
 
exit 0

########################################################################################
# man pages for the main executables used above

# Usage:
# Assemblytics delta output_prefix unique_length_required path_to_R_scripts
