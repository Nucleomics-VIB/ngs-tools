#!/bin/bash

# script: mapped2fastx.sh
# search a WGS readset for subset mapping to a chromosome(s) (eg mito)
#  create a new BWA index with a subset of chromosome(s)
#  align all reads to the new index and keep/filter mapped reads
#  pipe and convert back selected reads to fastx
#
# Stephane Plaisance (VIB-NC) 2017/11/30; v1.0

# required (tested versions)
# samtools (1.x, htslib)
# bwa (0.7.15+)

read -d '' usage <<- EOF
Usage: mapped2fastx.sh 
#   -i <reads.fastx (required)>
#   -r <reference (mutlifasta with target regions, required)>
#   -x <bwa mem read type (na, 'pacbio', 'ont2d'), opt)>
#   -o <output format ('fasta', 'fastq') (default to 'fasta', opt)>
#   -t <threads for processing (default to 4)>
#   -a <keep only mapped reads (default, otherwise filter them out)>
#   -h <show this help>
EOF

while getopts "i:r:x:o:t:ah" opt; do
  case $opt in
    i)
      reads=${OPTARG}
      ;;
    r)
      reference=${OPTARG}
      ;;
	x)
	  rtype=${OPTARG}
	  ;;
	o)
	  outformat=${OPTARG}
	  ;;
	t)
	  threads=${OPTARG}
	  ;;
	a)
	  keepmapped=1
	  ;;
    h)
      echo "${usage}"
      exit 0
      ;;      
    \?)
      echo "${usage}"
      exit 1
      ;;
    *)
      echo "${usage}" >&2
      exit 1
      ;;
  esac
done

# test if minimal arguments were provided
if [ -z "${reads}" ]; then
echo "# no read data provided!"
echo "${usage}"
exit 1
fi

# test if minimal arguments were provided
if [ -z "${reference}" ]; then
echo "# no reference provided!"
echo "${usage}"
exit 1
fi

# check if requirements are present
$( hash samtools 2>/dev/null ) || ( echo "# samtools not found in PATH"; exit 1 )
$( hash bwa 2>/dev/null ) || ( echo "# bwa not found in PATH"; exit 1 )

# define name for bwa index
idxname=$(basename ${reference%.*})

# reads are pacbio or ONT type
if [ -n "${rtype}" ]; then
	xval="-x ${rtype}"
else
	xval=""
fi

# check for both inputs
if [ ! -f "${reads}" ]; then
    echo "${reads} file not found!";
    exit 1
fi

if [ ! -f "${reference}" ]; then
    echo "${reference} file not found!";
    exit 1
fi

# output reads in format
outf=${outformat:-"fasta"}

# filter out or keep the spiked alignments
if [ -n "${keepmapped}" ]; then
	samfilter="-F 4"
	prefix="mapped"
else
	samfilter="-f 4"
	prefix="filtered"
fi

# what fastx output
if [ ${outf} == "fastq" ]; then
	fastx=${idxname}_${prefix}.fq.gz
else
	outf="fasta"
	fastx=${idxname}_${prefix}.fa.gz
fi

# threads
nthr=${threads:-4}

##################################################################
# aligning all reads to the index and saving mapped reads to fastx

outfolder="mapped2fastx"
mkdir -p ${outfolder}

# create new BWA index a ${reffasta} reference
cmd0="bwa index -p ${outfolder}/${idxname} ${reference}"
echo "# index build command: ${cmd0}"
eval ${cmd0}

# check for failure
if [ $? -ne 0 ]; then
	echo "# the bwa index could not be built, please check your inputs"
	exit 0
fi

# map all reads to the ${indexpath} bwa index
cmd1="bwa mem ${xval} -t ${threads} ${outfolder}/${idxname} ${reads}"
cmd2="samtools view -bS ${samfilter} - | samtools ${outf} - | bgzip -c > ${outfolder}/${fastx}"

##################################################################
# merging both commands to limit output to the de-spiked read file

mrgcmd="${cmd1} | ${cmd2}"
echo "# merged command: ${mrgcmd}"
eval ${mrgcmd}

exit 0
