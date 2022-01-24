#!/bin/bash

# script name: fasta2consensus.sh
# create consensus sequences from UMI clusters
# (single or multi-Fasta file of type '<barcode>_reads.fa')
# copy single Fasta records directly to output file
# derive consensus for multi-Fasta records
# rename output records from the barcode and supporting read count
#
# Stephane Plaisance (VIB-NC) 2021/12/30; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

version="1.0, 2022_12_30"

usage='# Usage: fasta2consensus.sh
# -i <multifasta NNNN_reads.fa file (1 file = 1 cluster)>
# -o <output consensus file (default to: "consensus_seqs.fa")>
# -w <workdir (location for intermediate files default to "/tmp")>
#
# NOTE: can be run in parallel on a folder of Fasta files with:
# ls ${f} | parallel -j 80 \
#   fasta2consensus.sh -i ${f}/{} -o $(basename ${f}_consensus_seqs.fa) -w temp
#
# script version '${version}'
# [-h for this help]'

while getopts "i:o:w:h" opt; do
  case $opt in
    i) opti=${OPTARG} ;;
    o) opto=${OPTARG} ;;
    w) optw=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# check executable present
declare -a arr=( "bioawk" "fold" "clustalo" "cons" )
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || ( echo "# required ${prog} not found in PATH"; exit 1 )
done

# check INPUTS
if [ -z "${opti}" ]
then
  echo "# provide a folder with UMI-cluster fasta files to be parsed"
  echo "${usage}" >&2
  exit 0
fi

# infile
infile=${opti}

# output file
outfile=${opto:-"consensus_seqs.fa"}

# workdir
workdir=${optw:-"/tmp"}
[ -d "${workdir}" ] || mkdir -p ${workdir}

# extract barcode from Fasta file name
barcode=$(basename ${infile%_reads.fa})

# count records in Fasta file
rcnt=$(bioawk -c fastx 'END{print NR}' ${infile})

# test if file contains one more sequences
if [[ "${rcnt}" == 1 ]]; then 

# single sequence Fasta file
# copy as is (60 char per line) and rename
bioawk -c fastx -v nm="${barcode}_${rcnt}" '{
print ">"nm"\n"$seq
}' ${infile} | fold -w60 >> ${outfile}

else

# multiple sequence Fasta file
# align and derive consensus
alifile=${workdir}/$(basename ${infile%.fa}.ali)
cnsfile=${workdir}/$(basename ${infile%.fa}_cons.fa)

clustalo --infile ${infile} \
  --outfmt clustal \
  --outfile ${alifile} \
  --force \
  --seqtype dna && \
cons -sequence ${alifile} \
  -outseq ${cnsfile} \
  -name ${barcode}_${rcnt} && \
cat ${cnsfile} >> ${outfile} 

# clean intermediate file
# && rm ${alifile} ${cnsfile}

fi
