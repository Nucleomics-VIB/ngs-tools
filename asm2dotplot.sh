#!/bin/bash
# script name: asm2dotplot.sh
# align two assemblies using minimap2
# create a dot-plot from the PAF results
# requirements: 
#   https://github.com/lh3/minimap2
# plot with pafCoordsDotPlotly_SP.R adpted from dotplotty
#   https://github.com/tpoorten/dotPlotly
# Rpackages for pafCoordsDotPlotly.R:
#   "optparse", "ggplot2", "plotly"
# Stephane Plaisance (VIB-NC) 2010/03/27; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

version="1.0, 2020-03-27"

usage='# Usage: asm2dotplot.sh
# -R <reference asm (.fa* or fa*.gz)>
# -Q <query asm (.fa* or fa*.gz)>
# -a <label for x-axis (default to R)>
# -b <label for y-axis (default to Q)>
# -o <output folder ("asm2dotplot")>
# -r <comma-separated ordered-list of reference records (undef)>
# -k <alt | plot only the first k reference records (undef)>
# -q <min query-length (1000)>
# -m <min alignment length (1000)>
# -l <horizontal lines on plot for separating scaffolds (FALSE)>
# -s <color alignments by % identity (FALSE)>
# -t <calculation of % identity for on-target alignments only (FALSE)>
# -T <threads for minimap2 job (4)>
# -h <this help>
# script version '${version}'
# [-h for this help]'

# default to no-color
colored=""
ontarget=""
drawlines=""

while getopts "R:Q:o:r:k:q:m:T:a:b:tslh" opt; do
  case $opt in
    R) reference=${OPTARG} ;;
    Q) query=${OPTARG} ;;
    o) optoutf=${OPTARG} ;;
    r) optrlst=${OPTARG} ;;
    k) optrlim=${OPTARG} ;;
    q) optmqlen=${OPTARG} ;;
    m) optmalen=${OPTARG} ;;
    s) colored=" -s" ;;
    t) ontarget=" -t" ;;
    l) drawlines=" -l" ;;
    a) optxlabel=${OPTARG} ;;
    b) optylabel=${OPTARG} ;;
    T) optthr=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# check executables present
declare -a arr=( "minimap2" "R" "pafCoordsDotPlotly_SP.R" )
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || \
    ( echo "# required ${prog} not found in PATH"; exit 1 )
done

# test if minimal arguments were provided
if [ -z "${reference}" ]
then
   echo "# no first assembly provided!"
   echo "${usage}"
   exit 1
fi

if [ ! -f "${reference}" ]; then
    echo "${reference} file not found!"
    exit 1
fi

if [ -z "${query}" ]
then
    echo "#  no second assembly provided!"
    echo "${usage}"
    exit 1
fi

if [ ! -f "${query}" ]; then
    echo "${query} file not found!";
    exit 1
fi

outfolder=${optoutf:-"asm2dotplot"}
mkdir -p ${outfolder}

refpfx=$(basename ${reference%%.fa*})
qrypfx=$(basename ${query%%.fa*})
outpfx=${qrypfx}_vs_${refpfx}
log=${outpfx}_log.txt

# create empty log file
cat /dev/null > ${outfolder}/${log}

thr=${optthr:-4}

###################
# minimap2 command
cmd="minimap2 -t ${thr} \
    -cx asm5 \
    ${reference} \
    ${query} \
    -o ${outfolder}/${outpfx}.paf"

echo "# ${cmd}" | tee -a ${outfolder}/${log}
eval ${cmd} 2>&1 | tee -a ${outfolder}/${log}

# check for failure
if [ $? -ne 0 ]; then
    echo "# the minimap2 command failed, please check your inputs"
    exit 0
fi

##################################
# pafCoordsDotPlotly_SP.R command
# the script has been modified to:
# * not remove the path before the output named
# * allow passing x_label and y_label through arguments

minqlen=${optmqlen:-1000}
minalilen=${optmalen:-1000}
plotsize=15

# reference limits
if [ -n "${optrlim}" ]; then
xlimit="-k ${optrlim}"
fi

# over-rules the previous choice
if [ -n "${optrlst}" ]; then
xlimit="-r ${optrlst}"
fi

xlabel=${optxlabel:-"$(basename ${reference%%.fa*})"}
ylabel=${optylabel:-"$(basename ${query%%.fa*})"}

cmd="pafCoordsDotPlotly_SP.R \
    -i ${outfolder}/${outpfx}.paf \
    -o ${outfolder}/${outpfx} \
    -q ${minqlen} \
    -m ${minalilen} \
    -p ${plotsize} \
    -a ${xlabel} \
    -b ${ylabel} \
    -x \
    ${xlimit} ${colored} ${ontarget} ${drawlines}"

echo "# ${cmd}" | tee -a ${outfolder}/${log}
eval ${cmd} 2>&1 | tee -a ${outfolder}/${log}
