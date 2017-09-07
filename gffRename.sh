#!/bin/bash
# script name: gffRename.sh
# Add prefix to ID and Parent subfields of the 9th GFF field
# required when merging several GFF files with shared ID or Parent names
#
# Stephane Plaisance (VIB-NC+BITS) 2017/09/07; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# get user input
if [[ $# != 2 ]]
then
echo "Usage: ${0##*/} <input-gff> <prefix>";
exit 1
fi

cat $1 | \
gawk -v s=${2} 'BEGIN{ FS="\t";OFS="\t"}	{ 
	if (/^#/) {print} else {
	gsub(/ID=/, "ID="s"_", $9);
	gsub(/Parent=/, "Parent="s"_", $9);
	print $1,$2,$3,$4,$5,$6,$7,$8,$9 } 
	}'
