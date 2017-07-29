#!/bin/bash

# script sort_tabular.sh
# sort tabular files
# accept any GNU sort options passed from -f between quotes
# by default removes all lines starting with '#' unless -h is provided
#
# Stephane Plaisance (VIB-NC) 2017/07/28; v1.0

usage="Usage: sort_tabular.sh -i <input.file> -f <filter expression> -h <keep header (default remove #-lines)>"

while getopts "i:f:h:" opt; do
  case $opt in
    i)
      infile=${OPTARG}
      ;;
    f)
      filterstring=${OPTARG}
      ;;
    h)
      keepheader=${OPTARG}
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
if [ -z "${infile}" ]
then
   echo "# no input provided!"
   echo "${usage}"
   exit 1
fi

# filter string or default
filter=${filterstring:-"-k 1V,1 -k 2n,2 -k 3n,3"}

if [ "${keepheader}" == "Yes" ]; then
grep "^#" ${infile}; grep -v "^#" ${infile} | sort ${filter}
else
grep -v "^#" ${infile} | sort ${filter}
fi
