#!/bin/env bash

# add this function to your .bashrc or alike
# split large compressed read files into chunks of fixed size
# suffix is a three digit counter starting with 000
# take compressed input and compress output with pigz
# keeps the read-in-pair suffix in outputs
# requires pigz installed or modification to use gzip
#
# Stephane Plaisance - VIB-NC Feb-25-2020 v1.0

usage="# splitreads <reads.fastq.gz> <reads per chunk; default 10000000>\n";
    if [ $# -lt 1 ]; then
        echo;
        echo ${usage};
        return;
    fi;

# threads for pigz (adapt to your needs)
thr=8

input=$1

# extract prefix and read number in pair
# this code is adapted to paired reads
base=$(basename ${input%.f*.gz})
pref=$(basename ${input%_?.f*.gz})
readn="${base#"${base%%_*}"}"

# 10M reads (4 lines each)
binsize=$((${2:-10000000}*4))

# split in bins of ${binsize}
echo "# splitting ${input} in chuncks of $((${binsize}/4)) reads"

cmd="zcat ${input} \
  | split \
    -a 3 \
    -d \
    -l ${binsize} \
    --numeric-suffixes \
    --additional-suffix ${readn} \
    --filter='pigz -p ${thr} > \$FILE.fq.gz' \
    - ${pref}_"

echo "# ${cmd}"
eval ${cmd}
