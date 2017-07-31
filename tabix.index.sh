#!/bin/bash

# script tabix.index.sh
# create tabix index from a bgzip compressed tabular file
#
# Stephane Plaisance (VIB-NC) 2017/07/29; v1.0

usage="Usage: tabix.index.sh -i <tabular.file.gz> 
	-b <begin-coordinate (2 for BED|GFF, 4 for SAM)>
	-c <comment char (#)>
	-e <end coordinate (3 for BED)>
	-p <preset (gff|bed|sam|vcf; do not apply when any of [-s, -b, -e, -c and -0] is defined> 
	-s <sequence name (1 for most file types)
	-0 <data is zero-based>"

while getopts "i:b:c:e:p:s:h0" opt; do
  case $opt in
    i)
      infile=${OPTARG}
      ;;
    0)
      zerob=" -0"
      ;;
    b)
      beginc=" -b "${OPTARG}
      ;;
    c)
      comment=" -c "${OPTARG}
      ;;
    e)
      endc=" -e "${OPTARG}
      ;;    
    p)
      preset=" -p "${OPTARG}
      ;;
    s)
      seqc=" -s "${OPTARG}
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

# handle -0
if [ -n "${zerob}" ]; then
zerobased=" -0"
fi

##################################
# ensure data is bgzip-compressed
##################################

if [[ ! ${infile} = *.gz ]]; then
echo "The input file will first be compressed with bgzip!"
bgzip -c ${infile} > ${infile}.gz && infile=${infile}".gz"
fi

##################
# run with preset
##################

if [ -n "${preset}" ]; then
# test incompatibility
if [[ -n ${seqc} || -n ${beginc} || -n ${endc} || -n ${comment} || -n ${zerob} ]]
then
echo "Argument -p conflicts with -s, -b, -e, -c, and -0"
exit 1
fi
# run with -p
tabix -f ${preset} ${infile}
fi

###################
# run with details
###################

if [ -z "${preset}" ]; then
# run with all available details
tabix -f ${zerobased} ${seqc} ${beginc} ${endc} ${comment} ${infile}
fi
