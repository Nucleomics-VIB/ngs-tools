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

# last edits
version="1.02, 2022-06-01"

usage='# Usage: get_ENSembl_ref.sh -o <organism (homo_sapiens)> 
#  -p <build name (Homo_sapiens.GRCh38)> 
#  -b <ensembl build (latest)>
#  -t <compression threads (1)>
#  -h <this help text>
# script version '${version}

# check executables present (not checking Picard)
declare -a arr=( "grep" "sed" "wget" "samtools")
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || ( echo "# required ${prog} not found in PATH"; exit 1 )
done

while getopts "o:p:b:t:h" opt; do
  case $opt in
    o) org_opt=${OPTARG} ;;
    p) pfx_opt=${OPTARG} ;;
    b) build_opt=${OPTARG} ;;
    t) nthr_opt=${OPTARG} ;;
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
nthr=${nthr_opt:-1}

# get current ENSembl build number
current=$(wget  -q -O - http://ftp.ensembl.org/pub/current_README | \
  grep "Ensembl Release" | \
  sed 's/^Ensembl Release \(.*\) Databases.$/\1/')
build=${build_opt:-${current}}

asmname=$(wget -q -O -  https://fungi.ensembl.org/Saccharomyces_cerevisiae/Info/Index?db=core | \
  grep "Genome assembly" | \
  sed 's/.*<h2>Genome assembly: <a href=\".*\">\(.*\)<\/a><\/h2>.*/\1/')

#########################################
# get the data from ensEMBL FTP         #
#########################################

baseurl=ftp.ensembl.org/pub/release-${build}

# genome assembly
asmlnk1="ftp://${baseurl}/fasta/${org}/dna/${pfx}.dna.primary_assembly.fa.gz"
asmlnk2="ftp://${baseurl}/fasta/${org}/dna/${pfx}.dna.toplevel.fa.gz"

# transcripts
translnk="ftp://${baseurl}/fasta/${org}/cdna/${pfx}.cdna.all.fa.gz"

# gene models
ann1lnk="ftp://${baseurl}/gff3/${org}/${pfx}.${build}.chr.gff3.gz"
ann2lnk="ftp://${baseurl}/gtf/${org}/${pfx}.${build}.chr.gtf.gz"

# variants (may not all exist!)
varlnk1="ftp://${baseurl}/variation/vcf/${org}/${org}.vcf.gz"
varlnk2="ftp://${baseurl}/variation/vcf/${org}/${org}_somatic.vcf.gz"
varlnk3="ftp://${baseurl}/variation/vcf/${org}/${org}_structural_variations.vcf.gz"
http://ftp.ensembl.org/pub/release-106/fasta/saccharomyces_cerevisiae/cdna/

outfolder="${pfx}.${build}"
mkdir -p ${outfolder}

# get the goods
for lnk in ${asmlnk1} ${asmlnk2} ${translnk} ${ann1lnk} ${ann2lnk} ${varlnk1} ${varlnk2} ${varlnk3}; do
wget -P ${outfolder} ${lnk}
done

#########################################
# decompress and create accessory files #
#########################################

for f in ${outfolder}/*.gz; do
bgzip -@${nthr} -d -f ${f}
done

# re-compress VCF and index
for v in ${outfolder}/*.vcf; do
bgzip -@${nthr} ${v} && \
  tabix -p vcf ${v}.gz
done

# set main assembly
if [ -f "${outfolder}/${pfx}.dna.primary_assembly.fa" ]; then
asmfile="${outfolder}/${pfx}.dna.primary_assembly.fa"
elif [ -f "${outfolder}/${pfx}.dna.toplevel.fa" ]; then
asmfile="${outfolder}/${pfx}.dna.toplevel.fa"
else
echo "assembly file not found, quitting!"
exit 1
fi

# create fasta fai
samtools faidx ${asmfile}
samtools faidx ${outfolder}/${pfx}.cdna.all.fa

# create .genome file and sequence dictionary
awk 'BEGIN{FS="\t"; OFS="\t"}{print $1, $2}' ${asmfile}.fai \
	> ${asmfile%.fa}.genome

java -jar $PICARD/picard.jar CreateSequenceDictionary \
	R=${asmfile} \
	O=${asmfile%.fa}.dict

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
