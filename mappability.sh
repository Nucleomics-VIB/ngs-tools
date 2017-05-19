#!/bin/bash
# script name: mappability.sh
# create GEM indices from a fasta reference
# create mappability track for a given read length
# convert to bigwig format for IGV
#
## Requirements:
# Fasta sequence names should not contain spaces !!!
# A strong unix computer (more threads=>shorter execution time)
# Jim Kent's Utils (https://github.com/ENCODE-DCC/kentUtils)
# GEM (http://algorithms.cnag.cat/wiki/GEM:Installation_instructions)
# samtools for indexing
# picard tools for producing the dict file
# a recent bash for capitalisation
#
# Stephane Plaisance (VIB-NC+BITS) 2017/05/19; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

version="1.0, 2017_05_19"

usage='# Usage: mappability.sh
# -i <reference assembly fasta>)> 
# -l <single read length for prediction (default to 100bps)>
# -p <prefix for output data (default to input file prefix)>
# -t <number of threads (default to 1)> 
# script version '${version}'
# [-h for this help]'

while getopts "i:l:p:t:h" opt; do
  case $opt in
    i) opti=${OPTARG} ;;
    l) optl=${OPTARG} ;;
    p) optp={OPTARGS} ;;
    t) optt={OPTARGS} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

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

# create new folder
basedir=$(dirname ${infile})
outdir="${basedir}/${pref}-mappability"
mkdir -p ${outdir}

# create reference size list
$( hash faSize 2>/dev/null ) || ( echo "## ERROR! faSize not found in PATH"; exit 1 )
echo "# creating reference size list"
faSize ${infile} -detailed > ${outdir}/${pref}.sizes

# create index
$( hash gem-indexer 2>/dev/null ) || ( echo "## ERROR! gem-indexer not found in PATH"; exit 1 )
echo "# creating GEM index ... be patient"
gem-indexer -i ${infile} \
	-o ${outdir}/${idxpref} \
	-c 'dna' \
	-T ${thr}

# create mappability track
echo "# creating GEM mappability data ... be patient"
$( hash gem-mappability 2>/dev/null ) || ( echo "## ERROR!  gem-mappability not found in PATH"; exit 1 )
gem-mappability -T ${thr} \
	-I ${outdir}/${idxpref}.gem \
	-l ${kmer} \
	-o ${outdir}/${pref}_${kmer}

# convert to wig
echo "# converting to WIG"
$( hash gem-2-wig 2>/dev/null ) || ( echo "## ERROR!  gem-2-wig not found in PATH"; exit 1 )
gem-2-wig -I ${outdir}/${idxpref}.gem \
	-i ${outdir}/${pref}_${kmer}.mappability \
	-o ${outdir}/${pref}_${kmer}

# convert to BigWig
echo "# converting to BigWIG"
$( hash wigToBigWig 2>/dev/null ) || ( echo "## ERROR!  wigToBigWig not found in PATH"; exit 1 )
wigToBigWig ${outdir}/${pref}_${kmer}.wig \
	${outdir}/${pref}.sizes \
	${outdir}/${pref}_${kmer}.bw

echo "# all done."
