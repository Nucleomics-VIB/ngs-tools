#!/bin/bash

# run_assemblytics.sh: produce a assembly QC based on mummer comparison to a reference
#
# Requirements:
# run on a unix computer
# nucmer (mummer4 apps) and Assemblytics both in $PATH
# R packages for plotting
# two related fasta references to be compared
#
# Stephane Plaisance (VIB-NC+BITS) 2017/11/27; v1.0
#
# v1.1: mummer4 nucmer is now multithreaded
#  updated Assemblytics 1.2 on 2020_03_27 from github
#  https://github.com/MariaNattestad/Assemblytics
#
# visit our Git: https://github.com/Nucleomics-VIB

# check parameters for your system
version="1.1.1, 2023_06_13"

# path to the Assemblytics scripts should be in PATH
# default_path_to_scripts="/opt/biotools/Assemblytics/scripts"

usage='# Usage: run_assemblytics.sh -x <reference (fasta)> -y <query-asm fasta)>
# script version '${version}'
# [optional: -o <result folder|Assemblytics_results>]
# [optional: -w <uniqseqlen|10000>]
# [optional: -m <min variant length|50>]
# [optional: -M <max variant length|10000>]
# [optional: -t <threads for alignment (4)>]
# [optional: -p <path to scripts|default set in the code>]
# [optional: -h <this help text>]'

# Initialize variables with defaults
outpathopt="Assemblytics_results"
uniqseqlen="10000"
minl="50"
maxl="10000"
threads="4"
path_to_scripts=""

# Parse command line arguments
while getopts "x:y:o:w:m:M:t:p:h" opt; do
  case $opt in
    x) assembly1=${OPTARG} ;;
    y) assembly2=${OPTARG} ;;
    o) outpathopt=${OPTARG} ;;
    w) uniqseqlen=${OPTARG} ;;
    m) minl=${OPTARG} ;;
    M) maxl=${OPTARG} ;;
    t) threads=${OPTARG} ;;
    p) path_to_scripts=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; echo "${usage}" >&2; exit 1 ;;
    :) echo "Option -${OPTARG} requires an argument." >&2; echo "${usage}" >&2; exit 1 ;;
  esac
done

# check executables present
declare -a arr=( "nucmer"  "R" "Assemblytics" )
for prog in "${arr[@]}"; do
hash ${prog} 2>/dev/null || ( echo "# required ${prog} not found in PATH"; exit 1 )
done

# defaults parameters
uniq_seq_len=${uniqseqlen}
minlen=${minl}
maxlen=${maxl}
thr=${threads}

# Validate numeric parameters
if ! [[ "$uniq_seq_len" =~ ^[0-9]+$ ]] || ! [[ "$minlen" =~ ^[0-9]+$ ]] || ! [[ "$maxlen" =~ ^[0-9]+$ ]] || ! [[ "$thr" =~ ^[0-9]+$ ]]; then
  echo "Error: numeric parameters must be positive integers" >&2
  exit 1
fi

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

# other parameters or defaults
outpath=${outpathopt:-"Assemblytics_results"}
asm1=$(basename "${assembly1%.f*}")
asm2=$(basename "${assembly2%.f*}")
deltabase="${asm2}_vs_${asm1}"
mkdir -p "${outpath}"

# build the nucmer4 command
nucmercmd="nucmer \
    --maxmatch \
    --threads=${thr} \
    --minmatch=100 \
    --mincluster=500 \
    \"${assembly1}\" \
    \"${assembly2}\" \
    --prefix=\"${outpath}/${deltabase}\" \
    > \"${outpath}/assemblytics-log.txt\" 2>&1"

# show and execute
echo "## pairwise alignments with nucmer:"
echo "# ${nucmercmd}"
eval ${nucmercmd}

# check for failure
if [ $? -ne 0 ]; then
    echo "# the nucmer command failed, please check your inputs"
    exit 1  # Changed from exit 0 to exit 1 to indicate error
fi

# Check if delta file was created
if [ ! -f "${outpath}/${deltabase}.delta" ]; then
    echo "# Error: Delta file was not created by nucmer"
    exit 1
fi

# build the Assemblytics command
cmd="Assemblytics \"${outpath}/${deltabase}.delta\" \
    \"${outpath}/${deltabase}\" \
    ${uniq_seq_len} \
    ${minlen} \
    ${maxlen} \
    >> \"${outpath}/assemblytics-log.txt\" 2>&1"

# show and execute
echo
echo "## Assemblytics analysis:"
echo "# ${cmd}"
eval ${cmd}

# Check Assemblytics exit status
if [ $? -ne 0 ]; then
    echo "# Error: Assemblytics command failed"
    exit 1
fi

echo "# Assemblytics analysis completed successfully"
exit 0

########################################################################################
# man pages for the main executables used above

# Usage:
# Assemblytics delta output_prefix unique_length_required min_size max_size
