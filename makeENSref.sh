#!/bin/bash
# script name: makeENSref.sh
# download human reference files from the ensembl ftp
# add indexes
## Requirements:
# wget to download files from the internet
# samtools for indexing
# picard tools for producing the dict file
# a recent bash for capitalisation
#
# Stephane Plaisance (VIB-NC+BITS) 2017/03/03; v1.0
#
# visit our Git: https://github.com/BITS-VIB

version="1.0, 2017_04_05"

usage='# Usage: makehumanref.sh
# -o <organism (default to <homo_sapiens>)> 
# -b <build number (default to <GRCh38>)> 
# -r <release number (default to 88)>
# script version '${version}'
# [-h for this help]'

while getopts "o:b:r:h" opt; do
  case $opt in
    o) opto=${OPTARG} ;;
    b) optb=${OPTARG} ;;
    r) optr=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# check file presence before download
function validate_url(){
  if [[ $(wget -S --spider $1 2>&1 | grep "File.*exists\.") ]]; then return 0; else return 1; fi
}

organism=${opto:-"homo_sapiens"}
build=${optb:-"GRCh38"}
release=${optr:-88}

# FTP base at ENSembl
baseurl=ftp://ftp.ensembl.org/pub/release-${release}/

# create new folder
mkdir -p ${build}.${release} && cd ${build}.${release}

# get genome
oriname=${organism^}.${build}.dna.primary_assembly.fa.gz
url=${baseurl}/fasta/${organism}/dna/${oriname}

# test for presence of a primary assembly
if `validate_url $url >/dev/null` ; then
fullpfx=${organism^}.${build}.${release}.dna.primary_assembly
else
# echo "## no primary-assembly found, trying with toplevel assembly"
oriname=${organism^}.${build}.dna.toplevel.fa.gz
url=${baseurl}/fasta/${organism}/dna/${oriname}
fullpfx=${organism^}.${build}.${release}.dna.toplevel

# test once more for good measure
if ! $(validate_url ${url} >/dev/null); then
echo "# reference was not found using our euristic!"
exit 0
fi
fi

wget ${baseurl}/fasta/${organism}/dna/${oriname}
gunzip -cd ${oriname} > ${fullpfx}.fa && rm ${oriname}
samtools faidx ${fullpfx}.fa
java -jar /opt/biotools/picard/picard.jar CreateSequenceDictionary R=${fullpfx}.fa O=${fullpfx}.dict
grep "^@SQ" ${fullpfx}.dict | awk '{split($2,name,":"); split($3,len,":"); print name[2]"\t"len[2]}' > ${fullpfx}.genome

# get cDNA
oriname=${organism^}.${build}.cdna.all.fa.gz
fullpfx=${organism^}.${build}.${release}.cdna.all
wget ${baseurl}/fasta/${organism}/cdna/${oriname}
gunzip -cd ${oriname} > ${fullpfx}.fa && rm ${oriname}
samtools faidx ${fullpfx}.fa

# get annotations
gtfname=${organism^}.${build}.${release}.gtf.gz
gffname=${organism^}.${build}.${release}.gff3.gz
wget ftp://ftp.ensembl.org/pub/release-${release}/gtf/${organism}/${gtfname} && gunzip -d ${gtfname}
wget ftp://ftp.ensembl.org/pub/release-${release}/gff3/${organism}/${gffname} && gunzip -d ${gffname}

# list folder for control
ls -lah
