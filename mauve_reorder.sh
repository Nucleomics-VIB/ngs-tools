#!/bin/bash

# mauve_reorder.sh: reorder contigs based on a reference assembly using Mauve
#
# Requirements:
# run on a unix computer installed with Mauve
# samtools installed
# draft assembly and reference genome fasta present
#
# Stephane Plaisance (VIB-NC+BITS) 2017/04/21; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# check parameters for your system
version="1.0, 2017_04_21"

usage='# Usage: mauve_reorder.sh -i <draft assembly> -r <reference assembly> -p <mauve path>
# script version '${version}'
# [optional: -o <result folder>]
# [optional: -m <available RAM|1G>]
# [optional: -h <this help text>]'

while getopts "i:r:o:p:m:h" opt; do
  case $opt in
    i) draftassembly=${OPTARG} ;;
    r) reference=${OPTARG} ;;
    p) mauvepath=${OPTARG} ;;
    o) outpath=${OPTARG} ;;
    m) memory=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# defaults
startts=$(date +%s)

# test if minimal arguments were provided
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

if [ -z "${mauvepath}" ]
then
	echo "# no path to Mauve.jar provided!"
	echo "${usage}"
	exit 1
fi

if [ ! -f "${mauvepath}/Mauve.jar" ]; then
    echo "Mauve.jar file not found at ${mauvepath}!";
    exit 1
fi

# check if requirements are present
$( hash samtools 2>/dev/null ) || ( echo "# samtools not found in PATH"; exit 1 )

# other parameters or defaults
destfolder=${outpath:-"mauve_ordered-${draftassembly}"}
mem=${memory:-"1G"}

# build the command
cmd="java -Xmx${mem} -cp ${mauvepath}/Mauve.jar \
  org.gel.mauve.contigs.ContigOrderer \
  -output ${destfolder} \
  -ref ${reference} \
  -draft ${draftassembly} \
  >>mauve-reorder_${draftassembly}-log.txt 2>&1"

# show and execute	
echo "# ${cmd}"
eval ${cmd}
 
# after completion, copy the final ordered assembly to <destfolder> and index it
if [ $? -eq 0 ]; then
    cp ${destfolder}/$(ls ${destfolder}/ | sort -r | head -1)/${draftassembly} \
  ${destfolder}/ordered-${draftassembly} && \
  samtools faidx ${destfolder}/ordered-${draftassembly}
else
    echo "Mauve ordering seems to have failed, please check the mauve-reorder_${draftassembly}-log.txt!"
fi

endts=$(date +%s)
dur=$(echo "${endts}-${startts}" | bc)
echo "Done in ${dur} sec"
