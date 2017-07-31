#!/bin/bash

# script tabix_search.sh
# search interval from a tabix indexed bgzip-compressed filebordercolor=0 .5 .5
#
# Stephane Plaisance (VIB-NC) 2017/07/29; v1.0

usage="Usage: tabix_search.sh -i <tabular.file.gz> -q <query interval (chr:start-end)>"

while getopts "i:q:h" opt; do
  case $opt in
    i)
      infile=${OPTARG}
      ;;
    q)
      query=${OPTARG}
      ;;
    h)
      echo ${usage}
      exit 0
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

if [[ ! ${infile} = *.gz ]]; then
echo "Error: input file should be compressed with bgzip!"
exit 1
fi

if [ -z "${infile}.tbi" ]; then
echo "# no tabix index found!"
echo "${usage}"
exit 1
fi

if [ -z "${query}" ]; then
echo "# no query interval provided!"
echo "${usage}"
exit 1
fi

tabix ${infile} ${query}