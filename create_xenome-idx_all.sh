#!/bin/bash

# script: create_xenome-idx_all.sh

workdir=/data1/Xenome
mkdir -p ${workdir} && cd ${workdir}

# get the references from UCSC
# requires razip for parallel compression
# requires samtools for indexing the fasta
#for pfx in mm39 hg38; do
#mkdir reference && cd reference
#wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/'${pfx}'/chromosomes/chr*'
#for f in $(ls *.fa.gz | sort -k 1V,1); do echo $f; gzip -cd $f >> ${pfx}.fa; done && rm -f *.fa.gz
#razip -c ${pfx}.fa > ${pfx}.fa.gz && samtools faidx ${pfx}.fa.gz
#done && cd ${workdir}

mouseref=${workdir}/references/mm39.fa
humref=${workdir}/references/hg38.fa

# kmers
kmers=25

# server settings (standard: 8 threads of 12GB each=96GB)
nthr=8
ram=36

# build index host=mouse graft=human
mkdir -p index && cd index
xenome index -v \
-M ${ram} \
-T ${nthr} \
--tmp-dir ${workdir} \
-l ./index_creation_log.txt \
-P idx \
-H ${mouseref} \
-G ${humref} \
-K ${kmers}
