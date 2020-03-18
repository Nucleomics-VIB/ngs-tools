#!/bin/env bash

# run_StarAlign2pass.sh
# takes all fastq_files from the local 'reads' folder
# align reads to a reference
# perform two passes:
# first pass to add existing junctions to the reference index
# second pass to align to the extended reference index
#
# Stéphane Plaisance - VIB-BITS - Feb-27-2020 v1.1
# including NC stringency options 
#   but keeping other defaults as on STAR forum

version="1.1, 2020_03_04"

workdir=${1:-"/data2/NC_projects/DNBSEQG400_validation/DNBSEQG400_eval3"}

cd ${workdir}

thr=84
readlen=100
idxlen=$((${readlen}-1))
refindex=$STAR_INDEXES/GRCh38.99
reffolder=$BIODATA/references
refgtf=${reffolder}/Homo_sapiens.GRCh38.99/Homo_sapiens.GRCh38.99.chr.gtf
reffasta=${reffolder}/Homo_sapiens.GRCh38.99/Homo_sapiens.GRCh38.dna.primary_assembly.fa

#####################################
# Nucleomics STAR optional settings #
#####################################

read -r -d '' STAR_OPTIONS <<'EOF'
    --outFilterMismatchNmax 10 \
    --outFilterMismatchNoverLmax 0.3 \
    --alignSJDBoverhangMin 3 \
    --alignSJoverhangMin 5 \
    --alignIntronMin 21 \
    --alignIntronMax 500000 \
    --outFilterMultimapNmax 10 \
    --outSJfilterOverhangMin 12 30 30 30 \
    --outWigType None \
    --outSAMprimaryFlag OneBestScore
EOF

######################################
# build STAR index if does not exist #
######################################

if [ -d "${refindex}" ]; then

echo "# STAR index already present, passing"

else

mkdir -p ${refindex}

cmd="STAR \
    --runMode genomeGenerate \
    --runThreadN ${thr} \
    --genomeDir ${refindex} \
    --genomeFastaFiles ${reffasta} \
    --sjdbGTFfile ${refgtf} \
    --sjdbOverhang ${idxlen}"

echo "# ${cmd}"
eval ${cmd}

fi


##############################################################
# first pass on all samples to collect all possible junctions
##############################################################

mkdir -p ${workdir}/STAR_mappings_PASS-1

for reads1 in ${workdir}/reads/*_1.fq.gz; do
    samplename=$(basename ${reads1/%1.f*.gz})
    reads2=${reads1/_1.fq.gz/_2.fq.gz}

    outpfx=${workdir}/STAR_mappings_PASS-1/${samplename}

    cmd="STAR \
        --runMode alignReads \
        --genomeDir ${refindex} \
        --genomeLoad LoadAndKeep \
        --readFilesCommand zcat \
        --readFilesIn ${reads1} ${reads2} \
        --outSAMtype None \
        --runThreadN ${thr} \
        --outFileNamePrefix ${outpfx}\
        ${STAR_OPTIONS}"

    echo "# first PASS alignment for $(basename ${reads1}) $(basename ${reads2})"
    echo "# ${cmd}"
    eval ${cmd}
done

# unload genome
cmd="STAR \
    --genomeDir ${refindex} \
    --genomeLoad Remove \
    --outSAMtype None \
    --outFileNamePrefix /dev/null/"

echo "# unloading the reference index"
echo "# ${cmd}"
eval ${cmd}

#############################
# merge and filter junctions
#############################

# 1. Filter out the junctions on chrM, those are most likely to be false.
# 2. Filter out non-canonical junctions (column5 == 0).
# 3. Filter out junctions supported by multimappers only (column7==0)
# 4. Filter out junctions supported by too few reads (e.g. column7<=2)

cmd="cat ${workdir}/STAR_mappings_PASS-1/*_SJ.out.tab \
    | awk 'BEGIN {OFS=\"\\t\"; strChar[0]=\".\"; strChar[1]=\"+\"; strChar[2]=\"-\";}
      {if((\$1!=\"M\") && (\$5>0) && (\$7>2)){print \$1,\$2,\$3,strChar[\$4]}}' \
    | sort -k 1V,1 -k 2n,2 -k 3n,3 \
    | uniq > ${workdir}/SJ.out.all.tab"

echo "# merging and filtering all junction files"
echo "# ${cmd}"
eval ${cmd}

##################################
# create new index with junctions
##################################

mkdir -p ${workdir}/SJ_index
outpfx=${workdir}/SJ_index/

cmd="STAR \
    --runMode genomeGenerate \
    --genomeDir ${workdir}/SJ_index \
    --genomeFastaFiles ${reffasta} \
    --sjdbGTFfile ${refgtf} \
    --sjdbFileChrStartEnd ${workdir}/SJ.out.all.tab \
    --runThreadN ${thr} \
    --sjdbOverhang ${idxlen} \
    --outFileNamePrefix ${outpfx}"

echo "# creating junction-aware reference index"
echo "# ${cmd}"
eval ${cmd}

################################################
# second pass on all samples and counting genes
################################################

mkdir -p ${workdir}/STAR_mappings_PASS-2

for reads1 in reads/*_1.fq.gz; do
    samplename=$(basename ${reads1/%1.f*.gz})
    reads2=${reads1/_1.fq.gz/_2.fq.gz}

    outpfx=${workdir}/STAR_mappings_PASS-2/${samplename}

    cmd="STAR \
        --runMode alignReads \
        --runThreadN ${thr} \
        --genomeDir ${workdir}/SJ_index \
        --readFilesCommand zcat \
        --readFilesIn ${reads1} ${reads2} \
        --outFileNamePrefix ${outpfx} \
        --sjdbFileChrStartEnd ${workdir}/SJ.out.all.tab \
        --outFilterType BySJout \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes Standard \
        --outSAMunmapped Within \
        --quantMode GeneCounts \
        ${STAR_OPTIONS}"

    echo "# second PASS alignment for $(basename ${reads1}) $(basename ${reads2})"
    echo "# ${cmd}"
    eval ${cmd}

done

# cleanup
rm -rf _STARtmp

