#!/bin/bash
# script name: varscan2mpileup2trio.sh
# call variants using varscan2 on a samtools mpileup
# process a trio in order father mother child
# use parallel to speed up the process
# merge, sort, compress, and index the resulting VCF
#
## Requirements:
# reference fasta used for mapping
# BAM file from the BWA mapping
# samtools, varscan2, vcftools, sort, bgzip, tabix
# varscan.jar should be present in a defined $VARSCAN folder path,
# if not available, please instal these dependencies first
#
# Stephane Plaisance, VIB-NC 2017/09/20; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

version="1.0, 2017_10_30"

# check if executables are present
declare -a arr=("samtools" "$VARSCAN/varscan.jar" "vcf-concat" "bgzip" "tabix" "sort")
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || ( echo "# required ${prog} not found in PATH"; exit 1 )
done

usage='# Usage: varscan2mpileup2trio.sh
# -r <reference assembly fasta> 
# -X <bam file with mappings against the same reference>
# -Y <bam file with mappings against the same reference>
# -Z <bam file with mappings against the same reference>
# -o <prefix for output data (default to input file prefix)>
# -m <max ram for each varscan2 jobs (default to 4G)>
# -t <max number of threads (default to 1)>
# -C | --min-coverage <Minimum read depth at a position to make a call (8)>
# -A | --min-reads2 <Minimum supporting reads at a position to call variants (2)>
# -Q | --min-avg-qual <Minimum base quality at a position to count a read (15)>
# -F | --min-var-freq <Minimum variant allele frequency threshold (0.01)>
# -H | --min-freq-for-hom <Minimum frequency to call homozygote (0.75)>
# -P | --p-value <Default p-value threshold for calling variants (99e-2)>
# -S | --strand-filter <Ignore variants with >90% support on one strand (1, unset with 0)>
# script version '${version}'
# [-h for this help]'

while getopts "r:X:Y:Z:o:m:t:C:A:Q:F:H:P:Sh" opt; do
  case $opt in
    r) ref=${OPTARG} ;;
    X) bam_father=${OPTARG} ;;
    Y) bam_mother=${OPTARG} ;;
    Z) bam_child=${OPTARG} ;;
    o) opt_o=${OPTARG} ;;
    m) opt_m=${OPTARG} ;;
    t) opt_t=${OPTARG} ;;
    C | --min-coverage) opt_mincov=${OPTARG} ;;
    A | --min-reads2) opt_minr2=${OPTARG} ;; 
    Q | --min-avg-qual) opt_avgqual=${OPTARG} ;;
    F | --min-var-freq) opt_minvarf=${OPTARG} ;;
    H | --min-freq-for-hom) opt_minfhom=${OPTARG} ;;
    P | --p-value) opt_pval=${OPTARG} ;;
    S | --strand-filter) opt_strand=${OPTARG} ;;
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

if [ -z "${bam_father}" || -z "${bam_mother}" || -z "${bam_child}"]
then
   echo "# Three BAM files are required in order: Father Mother Child!"
   echo "${usage}"
   exit 1
fi

# check if exists or die
[ -f ${bam_father} ] || ( echo "## ERROR! Father ${bam_father} input not found" ; exit 1 )
[ -f ${bam_mother} ] || ( echo "## ERROR! Mother ${bam_mother} input not found" ; exit 1 )
[ -f ${bam_child} ] || ( echo "## ERROR! Child ${bam_child} input not found" ; exit 1 )

# other defaults
outpref=${opt_o:-"$(basename ${bam_child%.*})"}
maxmem=${opt_m:-"4G"}
maxthr=${opt_t:-1}
mincvg=${opt_mincov:-8}
minread2=${opt_minr2:-2}
minavgq=${opt_avgqual:-15}
minvarf=${opt_minvarf:-0.01}
minfhom=${opt_minfhom:-0.75}
pval=${opt_pval:-"99e-03"}
strandfilter=${opt_strand:-1}

# create new folder based on child BAM path
basedir=$(dirname ${bam_child})
outdir="${basedir}/${outpref}"
mkdir -p ${outdir}

# run with parallel threads
cmd="samtools view -H ${bam_child} \
	| grep \"\@SQ\" \
	| sed 's/^.*SN://g' \
	| cut -f 1 \
	| xargs -I {} -n 1 -P ${maxthr} \
		sh -c \"samtools mpileup -f ${ref} -r \\\"{}\\\" ${bam_father} ${bam_mother} ${bam_child} \
		| java -Xmx${maxmem} -jar $VARSCAN/varscan.jar trio \
			--output-name ${outdir}/\\\"{}\\\".vcf \
			--min-coverage ${mincvg} \
			--min-reads2 ${minread2} \
			--min-avg-qual ${minavgq} \
			--min-var-freq ${minvarf} \
			--min-freq-for-hom ${minfhom}
			--p-value ${pval} \
			--strand-filter ${strandfilter}\""
			
echo "# ${cmd}"
eval ${cmd}

# concat all results
concat=${outdir}/merged_results.vcf

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
