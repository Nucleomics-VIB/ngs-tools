#!/bin/bash
# script name: genmap_mappability.sh
# create genmap indices from a fasta reference
# create mappability track for a given read length
#
## Requirements:
# genmap from https://github.com/cpockrandt/genmap
# bioawk for reference cleaning (https://github.com/lh3/bioawk)
# a recent bash for capitalisation
#
# Stephane Plaisance (VIB-NC+BITS) 2020/06/05; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

version="1.01, 2020_06_05"

usage='# Usage: genmap_mappability.sh
# -i <reference assembly fasta>
# -l <single read length for prediction (default to 100bps)>
# -p <prefix for output data (default to "genmap_mappability")>
# -e <max errors (0-4, default to 4)>
# -t <add detailed AND LARGE text output (off by default)>
# -f <count only for contigs starting with prefix (default to none, eg "chr")>
# script version '${version}'
# [-h for this help]'

while getopts "i:l:p:e:f:th" opt; do
  case $opt in
    i) opti=${OPTARG} ;;
    l) optl=${OPTARG} ;;
    p) optp=${OPTARG} ;;
    e) opte=${OPTARG} ;;
    t) optt=" -t ";;
    f) optf=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# test if minimal arguments were provided
infile=${opti}
if [ -z "${infile}" ]
then
  echo "# no fasta provided!"
  echo "${usage}"
  exit 1
fi

# check if exists or die
[ -f ${infile} ] || ( echo "## ERROR! Fasta reference file not found" ; exit 1 )

ref=$(basename ${infile%.fa*})
kmer=${optl:-100}
err=${opte:-4}
filter=${optf:-""}

####################
# custom functions #
####################

function countmappability() {
gawk -v min=$1 'BEGIN{FS="\t"; OFS=","; tot=0} \
  {if ($4>=min/10000) tot=tot+$3-$2} \
  END{printf "%d,%d\n", min/100, tot}' $2
  }
export -f countmappability

ranok () {
if [ $? -ne 0 ]
then
echo "!! something went wrong, quitting."
exit 1
fi
}
export -f ranok

#####################
# create outfolders #
#####################

outdir=${optp:-"$(pwd)/genmap_mappability"}
mkdir -p ${outdir}/tmp
export TMPDIR=${outdir}/tmp

#######################
# create index (once) #
#######################

if [ ! -f ${outdir}/genmap_${ref}_idx_created ]; then

echo "# creating genmap index for ${ref}"

genmap index \
-F ${infile} \
-I ${outdir}/genmap_${ref}_idx

touch ${outdir}/genmap_${ref}_idx_created
else
echo "# genmap index already present; delete genmap_${ref}_idx and re-run to reindex"
fi

#######################
# compute mappability #
#######################

if [ ! -f ${outdir}/genmap_${ref}_${kmer}_${err}_computed ]; then

echo "# computing mappability on ${ref} for kmer:${kmer} and max ${err} errors"
genmap map \
-I ${outdir}/genmap_${ref}_idx \
-O ${outdir}/genmap_${ref}_${kmer}_${err} \
-K ${kmer} \
-E ${err} \
-w \
-bg ${optt}

ranok

touch ${outdir}/genmap_${ref}_${kmer}_${err}_computed
else
echo "# mappability already computed for kmer=${kmer} and ${err} error(s), please delete ${outdir}/genmap_${kmer}_${err}_computed to recompute"
fi

#################
# filter counts #
#################

echo "# plotting mappability on ${ref} for kmer:${kmer} and max ${err} errors"

bedg=${outdir}/genmap_${ref}_${kmer}_${err}.bedgraph

if [ -n "${optf}" ]; then
echo "# counting only contigs starting by: ${optf}"
filter="_${optf}"
data=${outdir}/genmap_${ref}_${kmer}_${err}${filter}.bedgraph
# filter bedgraph data
gawk 'BEGIN{FS="\t"; OFS="\t"}{if ($1~/'${optf}'/) print $0}' \
${bedg} \
> ${data}
else
echo "counting all contigs in reference"
filter=""
data=${bedg}
fi

ranok

####################
# summaryze counts #
####################

echo "# summarizing base counts from ${data}"

counts=${outdir}/genmap_${ref}_${kmer}_${err}${filter}_counts.txt
echo "min-mappability_pc,bases" > ${counts}
(for (( i=1; $i<=10; i++ )); do echo 10000/$i |bc; done; echo 0) \
|parallel 'countmappability {} '${data} \
|sort --field-separator=',' -k 1n,1 \
>> ${counts}

ranok

####################
# plot from counts #
####################

echo "# plotting from summarised counts"

plot=${outdir}/mappability_${ref}_${kmer}_${err}${filter}_plot.png

gnuplot <<-EOFMarker
set datafile separator ","
set key autotitle columnhead
set style line 1 \
  linecolor rgb '#0060ad' \
  linetype 1 linewidth 2 \
  pointtype 7 pointsize 1.5
set yrange [*:100]
set xlabel "mappability % (kmer: '${kmer}')"
set ylabel "% of genome width above map%"
max(a,b)=a>b?a:b
a=-99999999.
set terminal pngcairo
set output '${plot}'
plot '${counts}' using 1:(a=max(a,column(2)),0/0) notitle, \
  '' using 1:(100*column(2)/a) notitle with linespoints linestyle 1
EOFMarker

echo "# results have been written to: mappability_${ref}_${kmer}_${err}_plot.png"
