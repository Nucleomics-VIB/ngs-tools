#!/bin/bash

# script: get_ENSembl_ref.sh
# Aim: download a set of eg 'homo sapiens' reference files
# NOTE: will not work for all organisms due to ensembl file tree inconsistency
# create accessory files and indices
# required Picard and other tools
#
# Stéphane Plaisance - VIB-Nucleomics Core - 2020-02-28 v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

usage='# Usage: get_ENSembl_ref.sh -o <organism (homo_sapiens)> -p <build name (Homo_sapiens.GRCh38)> -b <ensembl build (106)>
# script version '${version}'
# [optional: -h <this help text>]'

# check parameters for your system
version="1.01, 2021-05-11"

# check executables present (not checking Picard)
declare -a arr=( "grep" "sed" "wget" "samtools")
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || ( echo "# required ${prog} not found in PATH"; exit 1 )
done

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

# get latest build number live
latest=$(wget  -q -O - http://ftp.ensembl.org/pub/current_README | grep "Ensembl Release" | sed 's/^Ensembl Release \(.*\) Databases.$/\1/')
build=${build_opt:-${latest}}

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

# variants (may not all exist!)
varlnk1="ftp://${baseurl}/variation/vcf/${org}/${org}.vcf.gz"
varlnk2="ftp://${baseurl}/variation/vcf/${org}/${org}_somatic.vcf.gz"
varlnk3="ftp://${baseurl}/variation/vcf/${org}/${org}_structural_variations.vcf.gz"

outfolder="${pfx}.${build}"
mkdir -p ${outfolder}

# get the goods
for lnk in ${asmlnk} ${translnk} ${ann1lnk} ${ann2lnk} ${varlnk1} ${varlnk2} ${varlnk3}; do
wget -P ${outfolder} ${lnk}
done

#########################################
# decompress and create accessory files #
#########################################

for f in ${outfolder}/*.gz; do
gunzip -f ${f}
done

# recompress VCF and index
for v in ${outfolder}/*.vcf; do
vcf2index ${outfolder}/*.vcf
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
# human (build v38)
# http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
# http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gff3.gz

# human (build v33)
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gff3.gz
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz
# mouse (build vM24)
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.annotation.gff3.gz
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.basic.annotation.gtf.gz
