#!/bin/bash

# gepard_plot.sh: create a diagalign DNA plot from two assemblies / references
# get the code from https://github.com/univieCUBE/gepard
# ref: Krumsiek J, Arnold R, Rattei T, Bioinformatics 2007; 23(8): 1026-8. PMID: 17309896
#
# Stephane Plaisance (VIB-NC+BITS) 2017/04/21; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# check parameters for your system
# Post-install
# create an alias named 'gepard.jar' (pointing to the actual jar file) next to a matrices folder containing edna.mat and other matrices 
# refer with -p to that place
version="1.0, 2017_04_21"

usage='# Usage: gepard_plot.sh -x <reference assembly> -y <draft assembly> -p <path to gepard.jar and matrices>
# script version '${version}'
# [optional: -o <result path/prefix>]
# [optional: -w <word size:10>]
# [optional: -W <window size:0>]
# [optional: -J <java extra parameters (eg -Xmx1G, put between double quotes if it contains spaces)>
# [optional: -h <this help text>]'

while getopts "x:y:w:W:o:p:J:h" opt; do
  case $opt in
    x) reference=${OPTARG} ;;
    y) draftassembly=${OPTARG} ;;
    w) wordopt=${OPTARG} ;;
    W) windowopt=${OPTARG} ;;
    o) outfile=${OPTARG} ;;
    p) gepardpath=${OPTARG} ;;
    J) javaargs=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# defaults
startts=$(date +%s)

# test if minimal arguments were provided
if [ -z "${reference}" ]
then
	echo "# no reference assembly provided!"
	echo "${usage}"
	exit 1
fi

if [ ! -f "${reference}" ]; then
    echo "${reference} file not found!";
    exit 1
fi

if [ -z "${draftassembly}" ]
then
   echo "# no draft assembly provided!"
   echo "${usage}"
   exit 1
fi

if [ ! -f "${draftassembly}" ]; then
	echo "${draftassembly} file not found!"
	exit 1
fi

if [ -z "${gepardpath}" ]
then
	echo "# no path to gepard.jar provided!"
	echo "${usage}"
	exit 1
fi

if [ ! -f "${gepardpath}/gepard.jar" ]; then
    echo "gepard.jar file not found at ${gepardpath}!";
    exit 1
fi

# required matrix file
if [ ! -f "${gepardpath}/matrices/edna.mat" ]; then
    echo "edna.mat file not found at ${gepardpath}/matrices!";
    exit 1
fi

# other parameters or defaults
destfile=${outfile:-"gepard-$(basename "${draftassembly}" | cut -d. -f1)_vs_$(basename "${reference}" | cut -d. -f1)"}_w${wordopt:-10}_W${windowopt:-0}.png

# build the command
cmd="java ${javaargs:-""} -cp ${gepardpath}/gepard.jar org.gepard.client.cmdline.CommandLine \
	-seq1 ${reference} \
	-seq2 ${draftassembly} \
	-matrix ${gepardpath}/matrices/edna.mat \
	-word ${wordopt:-10} \
	-window ${windowopt:-0} \
	-format png \
	-outfile ${destfile}"

# show and execute	
echo "# ${cmd}"
eval ${cmd}

endts=$(date +%s)
dur=$(echo "${endts}-${startts}" | bc)
echo "Done in ${dur} sec"
