#!/bin/bash
# StarBuildRef.sh
# create reference genome index for STAR
# takes a fasta file and read length
# requires a lot of RAM for large genomes (>32GB for human)
#
# Stéphane Plaisance - VIB-BITS - Mar-22-2012 v1
# small edits, 2016-03-13 v1.1
# small edits, 2019-06-25 v1.2

version="1.2, 2019_06_25"

usage='# Usage: StarBuildRef.sh
# -i <reference.fasta>
# -g <matching gtf file>
# -l <read length (default: 75)>
# -o <STAR-indices root (if $STAR_INDEXES is not defined)>
# -h <this help text>
# script version '${version}

while getopts "o:i:g:l:h" opt; do
  case $opt in
    o)
      outdir=${OPTARG}
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

############# EDIT HERE ############
# hardware limits
thr=8
mram="64G"
##########STOP EDITING HERE ########

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
destfolder=${refpref}-index_${readlen}
mkdir -p ${star_folder}/${destfolder}

# save all to log for ezrror tracing ...
starttime=$(date +%s)
logfile=${star_folder}/star_${destfolder}-log.txt
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
eval ${cmd}

# report run duration
endtime=$(date +%s)
dur=$(echo ${endtime}-${starttime} | bc)
echo "# ended at $date: $endtime"
echo "Done in ${dur} sec"

exit 0
