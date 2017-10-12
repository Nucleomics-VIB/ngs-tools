#!/bin/bash
# script name: bcftoolsParallel.sh
# call variants using bcftools on a samtools mpileup
# use parallel to speed up the process
# merge, sort, compress, and index the resulting VCF
#
## Requirements:
# reference fasta used for mapping
# BAM file from the BWA mapping
# samtools, bcftools, vcftools, sort, bgzip, tabix
# if not available, please instal these dependencies first
#
# Stephane Plaisance, VIB-NC 2017/10/12; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

version="1.0, 2017_09_20"

# check if executables are present
declare -a arr=("samtools" "bcftools" "vcf-concat" "bgzip" "tabix" "sort")
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || ( echo "# required ${prog} not found in PATH"; exit 1 )
done

usage='# Usage: bcftoolsParallel.sh
# -r <reference assembly fasta> 
# -b <bam file with mappings against the same reference>
# -o <prefix for output data (default to input file prefix)>
# -t <max number of threads (default to 1)>
## samtools mpileup arguments 
# -B | --BAQ <enable BAQ (per-Base Alignment Quality) [default OFF]>
# -Q | --min-BQ <skip bases with baseQ/BAQ smaller than (13)>
# -d | --max-depth <max per-file depth; avoids excessive memory usage (250)>
# script version '${version}'
# [-h for this help]'

while getopts "r:b:o:t:QBh" opt; do
  case $opt in
    r | --ref) ref=${OPTARG} ;;
    b | --bam) bam=${OPTARG} ;;
    o | --out) opt_o=${OPTARG} ;;
    t | --threads) opt_t=${OPTARG} ;;
    B | --BAQ) opt_baq=true ;;
    Q | --min-BQ) opt_minbq=${OPTARG} ;;
    d | --max-depth) opt_maxdepth=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# test if minimal arguments were provided
if [ -z "${ref}" ]
then
   echo "# no reference fasta provided!"
   echo "${usage}"
   exit 1
fi

# check if exists or die
[ -f ${ref} ] || ( echo "## ERROR! ${ref} input not found" ; exit 1 )

if [ -z "${bam}" ]
then
   echo "# no bam mappings provided!"
   echo "${usage}"
   exit 1
fi

# check if exists or die
[ -f ${bam} ] || ( echo "## ERROR! ${bam} input not found" ; exit 1 )

# other defaults
outpref=${opt_o:-"$(basename ${bam%.*})"}
maxthr=${opt_t:-1}

# default to 'no-BAQ'
if [ -n opt_baq ]; then
nobaq=""
else
nobaq="-B"
fi

ninbq=${opt_minbq:-13}
maxdepth=${opt_maxdepth:-250}

# create new folder
basedir=$(dirname ${bam})
outdir="${basedir}/${outpref}"
mkdir -p ${outdir}

# run with parallel threads
cmd="samtools view -H ${bam} \
	| grep \"\@SQ\" \
	| sed 's/^.*SN://g' \
	| cut -f 1 \
	| xargs -I {} -n 1 -P ${maxthr} \
		sh -c \"samtools mpileup \
		${nobaq} \
		-Q${ninbq} \
		-d ${maxdepth} \
		-uf ${ref} \
		-r \\\"{}\\\" ${bam} \
		| bcftools call -cv \
		> ${outdir}/\\\"{}\\\".vcf\""
			
echo "# ${cmd}"
eval ${cmd}

# concat all results
concat=${outdir}/merged_results.vcf

if [ $? -eq 0 ]; then
vcf-concat $(find ${outdir} -type f ! -size 0 -name "*.vcf" \
	-not -name "merged_results.vcf" \
	-exec echo -n \'"{}"\'\  \; | tr '\n' ' ') \
	> ${concat}
fi

# sort compress index
if [ -f ${concat} ]; then
( grep ^"#" ${concat}; grep -v ^"#" ${concat} | \
	sort "-k 1V,1 -k 2n,2" ) | \
	bgzip -c > ${outdir}/${outpref}".vcf.gz" && \
	tabix -f -p vcf ${outdir}/${outpref}".vcf.gz"

echo "# all done."
else
echo "# merging vcf results failed, please check!"
fi
