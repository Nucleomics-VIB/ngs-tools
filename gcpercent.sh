#!/bin/bash
# script name: gcpercent.sh
# create GC percent track for a given bin width
# convert to TDF format for IGV
#
## Requirements:
# A strong unix computer (more threads=>shorter execution time)
# bedtools & igv
#
# Stephane Plaisance (VIB-NC) 2018/12/03; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

version="1.0, 2018_12_03"

usage='# Usage: gcpercent.sh
# -i <reference assembly fasta (not gzipped!)>)> 
# -w <window width for binning (default to 100bps)>
# -p <prefix for output data (default to input file prefix)>
# script version '${version}'
# [-h for this help]'

while getopts "i:w:p:h" opt; do
  case $opt in
    i) opti=${OPTARG} ;;
    w) optw=${OPTARG} ;;
    p) optp=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

#############################
# check executables present
declare -a arr=( "bioawk" "bedtools" "igvtools" )
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || ( echo "# required ${prog} not found in PATH"; exit 1 )
done

# declare variables
infile=${opti}
width=${optw:-100}
pref=${optp:-"$(basename ${infile%.fa*})"}

# test if minimal arguments were provided
if [ -z "${infile}" ]
then
   echo "# no fasta provided!"
   echo "${usage}"
   exit 1
fi

# check if exists or die
[ -f ${infile} ] || ( echo "## ERROR! Fasta input not found" ; exit 1 )

# create results folder
basedir=$(dirname ${infile})
refname=$(basename ${infile%.fa*})
outdir="${basedir}/${pref}-gcpercent"
mkdir -p ${outdir}

bioawk -c fastx '{ printf("%s\t%s\n",$name,length($seq)) }' ${infile} \
	> ${outdir}/${refname}.sizes

bedtools makewindows -g ${outdir}/${refname}.sizes \
  -w $width \
  > ${outdir}/${refname}.${width}bps.bed

bedtools nuc -fi ${infile} \
  -bed ${outdir}/${refname}.${width}bps.bed \
  > ${outdir}/${refname}_nuc_${width}bps.txt

gawk -v w=${width} 'BEGIN{FS="\t"; OFS="\t"}
  {
  if (FNR>1) {print $1,$2,$3,"GCpc_"w"bps",$5}
  }' ${outdir}/${refname}_nuc_${width}bps.txt \
  > ${outdir}/${refname}_GC_${width}.igv

# create binary version
igvtools toTDF \
	-z 5 \
	-f min,max,mean \
	${outdir}/${refname}_GC_${width}.igv \
	${outdir}/${refname}_GC_${width}.tdf \
	${infile}

echo "# all done."
