#!/bin/bash
# script name: varscan2mpileup2somatic.sh
# call variants using varscan2 on a samtools mpileup
# process in order normal tumour
# merge, sort, compress, and index the resulting VCF
#
## Requirements:
# reference fasta used for mapping
# BAM file from the BWA mapping
# samtools, varscan2, vcftools, sort, bgzip, tabix
# varscan.jar should be present in a defined $VARSCAN folder path,
# if not available, please instal these dependencies first
#
# Stephane Plaisance, VIB-NC 2017/11/06; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

version="1.0, 2017_11_06"

# check if executables are present
declare -a arr=("samtools" "$VARSCAN/varscan.jar" "vcf-concat" "bgzip" "tabix" "sort")
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || ( echo "# required ${prog} not found in PATH"; exit 1 )
done

usage='# Usage: varscan2mpileup2somatic.sh
# -r <reference assembly fasta> 
# -N <bam file with normal sample mappings against the same reference>
# -T <bam file with tumor sample mappings against the same reference>
# -o <prefix for output data (default to input file prefix)>
# -m <max ram for each varscan2 jobs (default to 4G)>
# -t <max number of threads (default to 1)>
# script version '${version}'
# [-h for this help]'

while getopts "r:N:T:o:h" opt; do
  case $opt in
    r) ref=${OPTARG} ;;
    N) bam_normal=${OPTARG} ;;
    T) bam_tumour=${OPTARG} ;;
    o) opt_o=${OPTARG} ;;
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

if [ -z "${bam_normal}" ] || [ -z "${bam_tumour}" ]
then
   echo "# no BAM files provided!"
   echo "${usage}"
   exit 1
fi

# check if exists or die
[ -f ${bam_normal} ] || ( echo "## ERROR! Father ${bam_normal} input not found" ; exit 1 )
[ -f ${bam_tumour} ] || ( echo "## ERROR! Mother ${bam_tumour} input not found" ; exit 1 )

# other defaults
outpref=${opt_o:-"$(basename ${bam_tumour%.*})"}

# create new folder based on child BAM path
basedir=$(dirname ${bam_child})
outdir="${basedir}/${outpref}"
mkdir -p ${outdir}

# run with parallel threads
cmd="samtools mpileup -B -f ${ref} -r \\\"{}\\\" ${bam_normal} ${bam_tumour} \
		| java -Xmx${maxmem} -jar $VARSCAN/varscan.jar somatic --mpileup ${outdir}/varscan_somatic"
			
echo "# ${cmd}"
eval ${cmd}

# post-process results
java -jar /opt/biotools/varscan/varscan.jar processSomatic ${outdir}/varscan_somatic.snp
java -jar /opt/biotools/varscan/varscan.jar processSomatic ${outdir}/varscan_somatic.indel

# concat all results
concat=${outdir}/varscan_somatic_merged.vcf

if [ $? -eq 0 ]; then
vcf-concat $(find ${outdir} -type f ! -size 0 -name "*.vcf" -not -name "merged_results.vcf" -exec echo -n \'"{}"\'\  \; | tr '\n' ' ') \
	> ${concat}
fi

# sort compress index
if [ -f ${concat} ]; then
( grep ^"#" ${concat}; grep -v ^"#" ${concat} | sort "-k 1V,1 -k 2n,2" ) | \
	bgzip -c > ${outdir}/${outpref}".vcf.gz" && \
	tabix -f -p vcf ${outdir}/${outpref}".vcf.gz"

echo "# all done."
else
echo "# merging vcf results failed, please check!"
fi

###### nothing living below this line #####
exit 0

# USAGE: VarScan somatic [normal_pileup] [tumor_pileup] [Opt: output] OPTIONS
#         normal_pileup - The SAMtools pileup file for Normal
#         tumor_pileup - The SAMtools pileup file for Tumor
#         output - Output base name for SNP and indel output
# 
# OPTIONS:
#++       --mpileup - input from a double mpileup (normal, tumour) 
#!!                   not compatible with some other parameters
#         --output-snp - Output file for SNP calls [output.snp]
#         --output-indel - Output file for indel calls [output.indel]
#         --min-coverage - Minimum coverage in normal and tumor to call variant [8]
#         --min-coverage-normal - Minimum coverage in normal to call somatic [8]
#         --min-coverage-tumor - Minimum coverage in tumor to call somatic [6]
#         --min-var-freq - Minimum variant frequency to call a heterozygote [0.10]
#         --min-freq-for-hom      Minimum frequency to call homozygote [0.75]
#         --normal-purity - Estimated purity (non-tumor content) of normal sample [1.00]
#         --tumor-purity - Estimated purity (tumor content) of tumor sample [1.00]
#         --p-value - P-value threshold to call a heterozygote [0.99]
#         --somatic-p-value - P-value threshold to call a somatic site [0.05]
#         --strand-filter - If set to 1, removes variants with >90% strand bias [0]
#         --validation - If set to 1, outputs all compared positions even if non-variant
#         --output-vcf - If set to 1, output VCF instead of VarScan native format