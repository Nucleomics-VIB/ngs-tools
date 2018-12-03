#!/bin/bash
# script name: mappability.sh
# create GEM indices from a fasta reference
# create mappability track for a given read length
# convert to bigwig format for IGV
#
## Requirements:
# A strong unix computer (more threads=>shorter execution time)
# Jim Kent's Utils (https://github.com/ENCODE-DCC/kentUtils)
# GEM (http://algorithms.cnag.cat/wiki/GEM:Installation_instructions)
# samtools for indexing and extracting reference sizes
# picard tools for producing the dict file
# a recent bash for capitalisation
#
# Stephane Plaisance (VIB-NC+BITS) 2017/05/19; v1.0
# add --complement argument
#
# visit our Git: https://github.com/Nucleomics-VIB

version="1.1, 2018_12_03"

usage='# Usage: mappability.sh
# -i <reference assembly fasta (not gzipped!)>)> 
# -l <single read length for prediction (default to 100bps)>
# -p <prefix for output data (default to input file prefix)>
# -t <number of threads (default to 1)>
# -b <index both strands <emulate|yes|no> (default=emulate)>
# script version '${version}'
# [-h for this help]'

while getopts "i:l:p:t:b:h" opt; do
  case $opt in
    i) opti=${OPTARG} ;;
    l) optl=${OPTARG} ;;
    p) optp=${OPTARG} ;;
    t) optt=${OPTARG} ;;
    b) optb=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

#############################
# check executables present
declare -a arr=( "samtools" "gem-indexer" "gem-mappability" "gem-2-wig" "wigToBigWig" )
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || ( echo "# required ${prog} not found in PATH"; exit 1 )
done

# test provided option
containsElement () {
  local e match="$1"
  shift
  for e; do [[ "$e" == "$match" ]] && return 0; done
  return 1
}

# declare variables
infile=${opti}

# test if minimal arguments were provided
if [ -z "${infile}" ]
then
   echo "# no fasta provided!"
   echo "${usage}"
   exit 1
fi

# check if exists or die
[ -f ${infile} ] || ( echo "## ERROR! Fasta input not found" ; exit 1 )

pref=${optp:-"$(basename ${infile%.*})"}
idxpref="${pref}_index"
kmer=${optl:-100}
thr=${optt:-1}

# test complement option
comp=${optb:-"emulate"}
ARRAY=( "emulate" "yes" "no" )
containsElement "${comp}" "${ARRAY[@]}"
if [ $? = 1 ]; then 
  echo "## ERROR! invalid -b argument: ${comp}"
  exit 1 
fi

# create new folder
basedir=$(dirname ${infile})
outdir="${basedir}/${pref}-mappability"
mkdir -p ${outdir}

# create reference size list
# faSize ${infile} -detailed > ${outdir}/${pref}.sizes
# faSize does not handle correctly long sequence names with spaces
# replaced by samtools and cut
echo "# creating reference size list"
samtools faidx ${infile} \
&& cut -f 1,2 ${infile}.fai > ${outdir}/${pref}.sizes

# create index
echo "# creating GEM index ... be patient"
gem-indexer -i ${infile} \
	-o ${outdir}/${idxpref} \
	--complement ${comp} \
	-c 'dna' \
	-T ${thr}

# create mappability track
echo "# creating GEM mappability data ... be patient"
gem-mappability -T ${thr} \
	-I ${outdir}/${idxpref}.gem \
	-l ${kmer} \
	-o ${outdir}/${pref}_${kmer}

# convert to wig
echo "# converting to WIG"
gem-2-wig -I ${outdir}/${idxpref}.gem \
	-i ${outdir}/${pref}_${kmer}.mappability \
	-o ${outdir}/${pref}_${kmer}

# convert to BigWig
echo "# converting to BigWIG"
wigToBigWig ${outdir}/${pref}_${kmer}.wig \
	${outdir}/${pref}.sizes \
	${outdir}/${pref}_${kmer}.bw

echo "# all done."
