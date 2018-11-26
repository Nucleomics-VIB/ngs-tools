#!/bin/bash

# usage: run_clumpify.sh -r <reads> -p <paired_reads>
# 
# find optical read duplicates and return counts
# opt: remove duplicates and create output read files
#
# Stephane Plaisance - VIB-Nucleomics Core - September-11-2017 v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# requirements
# BBTools clumpify.sh from https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/

version="1.0, 2018_11_26"

# <dedup distance
#     NextSeq      40  (and spany=t)
#     HiSeq 1T     40
#     HiSeq 2500   40
#     HiSeq 3k/4k  2500
#     Novaseq      12000
#
# command default Type presets
# -T N <Nextseq shorthand for 'dupedist=40 dedupe=t optical=t spany=t'>
# -T H <HiSeq shorthand for  'dupedist=2500 dedupe=t optical=t'>
# -T V <Novaseq shorthand for  'dupedist=12000 dedupe=t optical=t'>

usage='# Usage: run_clumpify.sh 
# -r <reads (fastq or gz archive)>
# -p <opt: paired-reads (fastq or gz archive)>
# -o <opt: <prefix>_summary.txt (default to basename of <reads>)>
# -T <one out of N|H|V (N=NextSeq, H=HiSeq3/4k, V=NovaSeq)>
# -w <opt: write deduplicated read set(s) (default to /dev/null)>
# -h <this help text>
# script version '${version}

while getopts "r:p:o:T:wh" opt; do
  case $opt in
    r) opt_reads=${OPTARG} ;;
    p) opt_paired=${OPTARG} ;;
    o) opt_prefix=${OPTARG} ;;
    T) opt_type=${OPTARG} ;;
    w) opt_write=1 ;;
    t) opt_tmp=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

#############################
# check executables present
declare -a arr=( "clumpify.sh" )
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || ( echo "# required ${prog} not found in PATH"; exit 1 )
done

clumpify=$(which clumpify.sh)

#############################
# reads
if [ -z "${opt_reads}" ]; then
   echo
   echo "# no read file provided!"
   echo "${usage}"
   exit 1
fi

if [ ! -f "${opt_reads}" ]; then
	echo
	echo "${reads} file not found!"
	exit 1
fi

if [ -n "${opt_paired}" ] && [ ! -f "${opt_paired}" ]; then
   echo
   echo "# paired reads provided but file not found!"
   echo "${usage}"
   exit 1
fi

##############################
# paired read command
if [ -n "${opt_paired}" ]; then
   cmd="${clumpify} in=${opt_reads} in2=${opt_paired}"
else
   # single read command
   cmd="${clumpify} in=${opt_reads}"
fi

#############################
# write out or not
if [ -n "${opt_paired}" ]; then
    if [ -n "${opt_write}" ]; then
        out2="dedup_$(basename ${opt_paired})"
    else
        # create symlink in tmp to /dev/null
        ln -s /dev/null ./out2lnk
        out2="./out2lnk"
    fi
    cmd="${cmd} out2=${out2}"
fi

if [ -n "${opt_write}" ]; then
    out="dedup_$(basename ${opt_reads})"
else
    # create symlink in tmp to /dev/null
    ln -s /dev/null ./outlnk
    out="./outlnk"
fi

cmd="${cmd} out=${out}"

#############################
# Platform specific presets
if [ -z "${opt_type}" ]; then
   echo "# no platform type provided!"
   echo "${usage}"
   exit 1
fi

opts=''
case ${opt_type} in
N)
  opts='dupedist=40 dedupe=t optical=t spany=t'
  ;;
H)
  opts='dupedist=2500 dedupe=t optical=t spany=f'
  ;;
V)
  opts='dupedist=12000 dedupe=t optical=t spany=f'
  ;;
*)
  Message="type must be of N|H|V"
  echo "${usage}"
  exit 1
  ;;
esac

# because out /tmp partition is so small
cmd="$cmd ${opts} tmpdir=. usetmpdir=t"

#############################
# defaults summary name
prefix=${opt_prefix:-$(basename $opt_reads)}
# remove extension(s)
prefix=${prefix%.gz}
prefix=${prefix%.fa*}
# 2 output files
log=${prefix}_clumpify-log.txt
summary=${prefix}_summary.txt

# save stderr to file too
exec 2> >(tee -a ${log})

echo "# ${cmd}"
eval ${cmd}

# happy ending
if [ $? -eq 0 ]; then
    echo
    echo "# results saved in ${summary}"
    echo
	(echo "# "${prefix}"; tail -6 ${log} | head -5) | tee -a ${summary}
    # cleanup links
    if [ -z "${opt_write}" ]; then
        unlink ${out}
        unlink ${out2}
    fi
fi
