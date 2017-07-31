#!/bin/bash

# script tabix_indexer.sh
# create tabix index from a bgzip compressed tabular file
#
# Stephane Plaisance (VIB-NC) 2017/07/29; v1.0

read -d '' usage <<- EOF
Usage: tabix_indexer.sh -i <tabular.file.gz> 
#	-b <begin-coordinate (2 for BED|GFF, 4 for SAM)>
#	-c <comment char (#)>
#	-e <end coordinate (3 for BED)>
#	-p <preset (gff|bed|sam|vcf; do not apply when any of [-s, -b, -e, -c and -0] is defined> 
#	-s <sequence name (1 for most file types)
#	-0 <data is zero-based (Undef=No as default)>
EOF

while getopts "i:b:c:e:p:s:h0" opt; do
  case $opt in
    i)
      infile=${OPTARG}
      ;;
    0)
      zerob=${OPTARG}
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
      echo "${usage}"
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

# comment char
comment_char="${comment:-"#"}"

# handle -0
if [ -n "${zerob}" ]; then
zerobased=" -0"
fi

##################
# run with preset
##################

if [ -n "${preset}" ]; then

#######################
# test incompatibility
#######################

if [[ -n ${seqc} || -n ${beginc} || -n ${endc} || -n ${comment} || -n ${zerob} ]]
then
echo "Argument -p conflicts with -s, -b, -e, -c, and -0"
exit 1
fi

##############
# run with -p
##############

##################################
# ensure data is bgzip-compressed
##################################

if [[ ! ${infile} = *.gz ]]; then
echo "The input file will first be sorted and compressed with bgzip!"

# handle sorting for each type
case $preset in
  gff)
  filter="-k 1V,1 -k 4n,4 -k 5n,5"
  ;;
  bed)
  filter="-k 1V,1 -k 2n,2 -k 3n,3"
  ;;
  sam)
  filter="-k 3V,3 -k 4n,4"
  ;;
  vcf)
  filter="-k 1V,1 -k 2n,2"
  ;;
esac

# create sorted archive
( grep ^${comment_char} ${infile}; grep -v ^${comment_char} ${infile} | sort ${filter} ) | \
	bgzip -c > ${infile}".gz"

# adapt infile name
infile=${infile}".gz"
fi

tabix -f ${preset} ${infile}
fi

###################
# run with details
###################

if [ -z "${preset}" ]; then

##################################
# ensure data is bgzip-compressed
##################################

if [[ ! ${infile} = *.gz ]]; then
echo "REM: The input file will first be sorted and compressed with bgzip!"

# filter based on coordinate columns
if [ -z "${seqc}" ]; then
echo "You need to specify at least sequence-name and start columns to sort the data and index it!"
exit 1
fi

# build filter string for sorting
filter="-k ${seqc: -1}V,${seqc: -1}"

if [ -n "${beginc}" ]; then
filter=$filter" -k ${beginc: -1}n,${beginc: -1}"
fi

if [ -n "${endc}" ]; then
filter=$filter" -k ${endc: -1}n,${endc: -1}"
fi

( grep ^["${comment_char}"] ${infile}; grep ^[^"${comment_char}"] ${infile} | sort ${filter} ) | \
	bgzip -c > ${infile}".gz"

# adapt infile name
infile=${infile}".gz"
fi

# run with all available details
indexcmd="tabix -f ${zerobased} ${seqc} ${beginc} ${endc} -c \"${comment_char}\" ${infile}"
eval ${indexcmd}
fi
# end run with details
