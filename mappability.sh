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
# bioawk for reference cleaning (https://github.com/lh3/bioawk)
# a recent bash for capitalisation
#
# Stephane Plaisance (VIB-NC+BITS) 2017/05/19; v1.0
# add --complement argument
# add bedgraph and plot data
#
# visit our Git: https://github.com/Nucleomics-VIB

version="1.3, 2020_04_20"

usage='# Usage: mappability.sh
# -i <reference assembly fasta (not gzipped!)>)> 
# -l <single read length for prediction (default to 100bps)>
# -p <prefix for output data (default to input file prefix)>
# -t <number of threads (default to 1)>
# -b <index both strands <emulate|yes|no> (default=emulate)>
# -m <min mappability cutoff for final filtering (default to 0.2)>
# script version '${version}'
# [-h for this help]'

while getopts "i:l:p:t:b:m:h" opt; do
  case $opt in
    i) opti=${OPTARG} ;;
    l) optl=${OPTARG} ;;
    p) optp=${OPTARG} ;;
    t) optt=${OPTARG} ;;
    b) optb=${OPTARG} ;;
    m) optm=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

#############################
# check executables present
declare -a arr=( "bioawk" "gem-indexer" "gem-mappability" "gem-2-wig" "wigToBigWig" "bigWigToBedGraph")
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
refname=$(basename ${infile%.fa*})
outdir="${basedir}/${pref}-mappability"
mkdir -p ${outdir}

# clean reference fasta file (names with space or comments)
echo "# creating a clean reference and size list"
bioawk -c fastx '{printf(">%s\n%s\n",$name,$seq)}' ${infile} \
	> ${outdir}/${refname}.fa

bioawk -c fastx '{ printf("%s\t%s\n",$name,length($seq)) }' ${infile} \
	> ${outdir}/${refname}.sizes

# create index
echo "# creating GEM index ... be patient"
gem-indexer -i ${outdir}/${refname}.fa \
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
	${outdir}/${refname}.sizes \
	${outdir}/${pref}_${kmer}.bw

# convert to bedGraph & BED
echo "# converting to BED"
bigWigToBedGraph ${outdir}/${pref}_${kmer}.bw \
	${outdir}/${pref}_${kmer}.bedGraph
	
# keep only regions with mappability > min
min=${optm:-0.2}
gawk -v min=${min} 'BEGIN{FS="\t"; OFS="\t"}
	{if ($4>=min) print $1 $2 $3 "mappability_gt_"min $4}' \
	${outdir}/${pref}_${kmer}.bedGraph \
	> ${outdir}/${pref}_${kmer}_gt_${min}.bed

# produce counts for plotting
# count bases above mappability threshold in range 0:100% by 5% steps
lims=( $(seq 0 5 100) )
echo "min-mappability_pc,bases" > ${outdir}/counts.txt 
# loop and add counts to file
for thr in ${lims[*]}; do
gawk -v min=${thr} 'BEGIN{FS="\t"; OFS="\t"; tot=0} \
	{if ($4>=min/100) tot=tot+$3-$2}\
	END{print min","tot}' \
	${outdir}/${pref}_${kmer}.bedGraph
done >> ${outdir}/counts.txt

# plot it to png
gnuplot <<-EOFMarker
set datafile separator ","
set key autotitle columnhead
set style line 1 \
    linecolor rgb '#0060ad' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.5
set yrange [*:100]
set xlabel "mappability %"
set ylabel "% above threshold"
max(a,b)=a>b?a:b
a=-99999999.
set terminal pngcairo
set output '${outdir}/mappability_plot.png'
plot '${outdir}/counts.txt' using 1:(a=max(a,column(2)),0/0) notitle, \
    '' using 1:(100*column(2)/a) notitle with linespoints linestyle 1
EOFMarker

echo "# all done."
