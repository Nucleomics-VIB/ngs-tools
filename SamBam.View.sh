#!/bin/bash

# script SamBam.view.sh
# view [top N-lines (default to 50)] or full SAM or BAM data 
#
# Stephane Plaisance (VIB-NC) 2017/07/29; v1.0

usage="Usage: SamBam.view.sh -i <SAM¬BAM.file>
 -t <only first N lines (default to 50; 0 for full data)>
 -h <include header>
 -H <show only header>"

while getopts "i:t:hH" opt; do
  case $opt in
    i)
      infile=${OPTARG}
      ;;
    h)
      showheader=" -h"
      ;;    
    H)
      headeronly=" -H"
      ;;  
    t)
      top=${OPTARG}
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
if [ -z "${infile}" ]; then
echo "# no input provided!"
echo "${usage}"
exit 1
fi

if [ ${top} == 0 ]; then
samtools view ${showheader} ${headeronly} ${infile}
else
samtools view ${showheader} ${headeronly} ${infile} | head -${top}
fi