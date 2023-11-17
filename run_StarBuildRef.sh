#!/bin/bash
# StarBuildRef.sh
# create reference genome index for STAR
# takes a fasta file and read length
# requires a lot of RAM for large genomes (>32GB for human)
#
# Stéphane Plaisance - VIB-BITS - Mar-22-2012 v1
# small edits, 2016-03-13 v1.1
# small edits, 2019-06-25 v1.2
# small edits, 2020-02-28 v1.3
# small edits, 2023-11-03 v1.4

version="1.4, 2023_11_03"

usage='# Usage: StarBuildRef.sh
# -i <reference.fasta>
# -g <matching gtf file>
# -l <read length (default: 75)>
# -o <STAR-indices root folder (if $STAR_INDEXES is not defined)>
# -T <title for this new STAR-indices folder (or defined from inputs)
# -m <memory (default: 64G)>
# -t <threads (default: 8)>
# -h <this help text>
# script version '${version}

# Default values
mram="64G"
thr=8

while getopts "o:i:g:l:T:m:t:h" opt; do
  case $opt in
    o)
      outdir=${OPTARG}
      ;;
    T)
      outfolder=${OPTARG}
      ;;
    i)
      reffasta=${OPTARG}
      ;;
    g)
      usrgtf=${OPTARG}
      ;;
    l)
      readlen=${OPTARG:-75}
      ;;
    m)
      mram=${OPTARG:-64G}
      ;;
    t)
      thr=${OPTARG:-8}
      ;;
    h)
      echo "${usage}" >&2
      exit 0
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      echo "${usage}" >&2
      exit 1
      ;;
    *)
      echo "this command requires arguments, try -h" >&2
      exit 1
      ;;
  esac
done

# variables
star_folder=${outdir:-"$STAR_INDEXES"}

if [ -z "$reffasta" ]
then
   echo "# Error: ${reffasta} fasta file not found!";
   echo "${usage}" >&2
   exit 1
fi

# check for GTF annotations
refpref=$(basename ${reffasta%%.fa*})
refgtf=${usrgtf:-${reffasta/.fa/.gtf}}
echo "# fasta : "${reffasta}
echo "# gtf   : "${refgtf}
echo "# prefix: "${refpref}

# test fasta presence
if [ ! -f ${reffasta} ]; then
    echo "# Error: ${reffasta} reference sequence not found!";
    echo "${usage}" >&2
    exit 1
fi

# test gtf presence
if [ ! -f ${refgtf} ]; then
    echo "# Error: ${refgtf} annotations not found!";
    echo "${usage}" >&2
    exit 1
fi

# test output folder exists
if [ ! -d ${star_folder} ]; then
    echo "# Error: ${star_folder} folder not found!";
    exit 1
fi

# idx length from provided read length
idxlen=$( echo ${readlen}-1 | bc)

# test valid value
re='^[0-9]+$'
if ! [[ ${idxlen} =~ $re ]] ; then
    echo "# Error: ${idxlen} is not a number!"; 
    exit 1
fi

# create folder for star files
destfolder=${outfolder:-"${refpref}-index_${readlen}"}
mkdir -p ${star_folder}/${destfolder}

# save all to log for ezrror tracing ...
starttime=$(date +%s)
logfile=${star_folder}/${destfolder}/starindex_build-log.txt
exec &> >(tee -a "${logfile}")

# build STAR index
cmd="STAR \
    --outFileNamePrefix ${star_folder}/${destfolder} \
    --runThreadN ${thr} \
    --runMode genomeGenerate \
    --genomeDir ${star_folder}/${destfolder} \
    --genomeFastaFiles ${reffasta} \
    --sjdbGTFfile ${refgtf} \
    --sjdbOverhang ${idxlen}"

echo "# ${cmd}"
eval ${cmd} && \
    mv ${star_folder}/${destfolder}Log.out ${star_folder}/${destfolder}/

# report run duration
endtime=$(date +%s)
dur=$(echo ${endtime}-${starttime} | bc)
echo "# ended at $date: $endtime"
echo "Done in ${dur} sec"

exit 0
