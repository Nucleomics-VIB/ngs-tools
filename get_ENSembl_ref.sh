#!/bin/bash

# script: get_ENSembl_ref.sh
# Aim: download a set of homo sapiens reference files
# create accessory files and indices
#
# Stéphane Plaisance - VIB-Nucleomics Core - 2020-02-28 v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

usage='# Usage: get_ENSembl_ref.sh -o <organism (homo_sapiens)> -p <build name (Homo_sapiens.GRCh38)> -b <ensembl build (99)>
# script version '${version}'
# [optional: -h <this help text>]'

# check parameters for your system
version="1.0, 2020-02-28"

while getopts "o:p:b:h" opt; do
  case $opt in
    o) org_opt=${OPTARG} ;;
    p) pfx_opt=${OPTARG} ;;
    b) build_opt=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# ref genome default settings
# org=${org_opt:-"mus_musculus"}
# pfx=${pfx_opt:-"Mus_musculus.GRCm38"}
org=${org_opt:-"homo_sapiens"}
pfx=${pfx_opt:-"Homo_sapiens.GRCh38"}
build=${build_opt:-"99"}

#########################################
# get the data from ensEMBL FTP         #
#########################################

baseurl=ftp.ensembl.org/pub/release-${build}

# genome assembly
asmlnk="ftp://${baseurl}/fasta/${org}/dna/${pfx}.dna.primary_assembly.fa.gz"

# transcripts
translnk="ftp://${baseurl}/fasta/${org}/cdna/${pfx}.cdna.all.fa.gz"

# gene models
ann1lnk="ftp://${baseurl}/gff3/${org}/${pfx}.${build}.chr.gff3.gz"
ann2lnk="ftp://${baseurl}/gtf/${org}/${pfx}.${build}.chr.gtf.gz"

outfolder="${pfx}.${build}"
mkdir -p ${outfolder}

# get the goods
for lnk in ${asmlnk} ${translnk} ${ann1lnk} ${ann2lnk}; do
wget -P ${outfolder} ${lnk}
done

#########################################
# decompress and create accessory files #
#########################################

for f in ${outfolder}/*.gz; do
gunzip ${f}
done

# create fasta fai
samtools faidx ${outfolder}/${pfx}.dna.primary_assembly.fa
samtools faidx ${outfolder}/${pfx}.cdna.all.fa

# create .genome file and sequence dictionary
awk 'BEGIN{FS="\t"; OFS="\t"}{print $1, $2}' ${outfolder}/${pfx}.dna.primary_assembly.fa.fai \
	> ${outfolder}/${pfx}.dna.primary_assembly.genome

java -jar $PICARD/picard.jar CreateSequenceDictionary \
	R=${outfolder}/${pfx}.dna.primary_assembly.fa \
	O=${outfolder}/${pfx}.dna.primary_assembly.dict


exit 0

# GenCode
# human (build v33)
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gff3.gz
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz
# mouse (build vM24)
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.annotation.gff3.gz
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.basic.annotation.gtf.gz
