#!/bin/bash

# compare VCF calls to same genome VCF goldstandard reference
# using Hap.py (docker)
#
# Requirements:
# your user is part of the docker group to avoid sudo
# run on a unix computer installed with
#   docker image: jmcdani20/hap.py:v0.3.12
#   # get with: sudo docker pull jmcdani20/hap.py:v0.3.12
#   VCF goldstandard for the same genome (GIAB or other caller VCF)
#   VCF calls (eg from Deepvariants)
#   reference FASTA used for mapping and calling
#
# this script will write results in the current folder
#
# Stephane Plaisance (VIB-NC) 2022/06/01; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# version number
version="1.0, 2023_06_15"

usage='# Usage: run_Hap_py.sh -v <VCF calls> -g <VCF reference (goldstandard)> -r <FASTA reference>
# script version '${version}'
# [optional: -n <outfolder name|prefix of VCF-calls>]
# [optional: -t <threads|4>]'

while getopts "v:g:r:n:t:h" opt; do
  case $opt in
    v) opt_calls=${OPTARG} ;;
    g) opt_gold=${OPTARG} ;;
    r) opt_ref=${OPTARG} ;;
    n) opt_name=${OPTARG} ;;
    t) opt_thr=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2;
       exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; 
       exit 1 ;;
  esac
done

# set defaults
nthr=${opt_thr:-4}

# check if all dependencies are present
declare -a arr=("samtools")
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || \
  ( echo "# required ${prog} not found in PATH"
exit 1 )
done

# Hap.py docker
IMAGE="jmcdani20/hap.py"
BIN_VERSION="v0.3.12"
img="${IMAGE}:${BIN_VERSION}"

# check for the docker image
res=$(docker images | grep ${img} | cut -d " " -f 1)
if [ ! ${res} == "${IMAGE}" ]; then
  echo "docker image not found"
  exit 1
fi

# check for VCF calls to be tested
# test if minimal arguments were provided
if [ -z "${opt_calls}" ]; then
   echo "# no VCF calls provided!"
   echo "${usage}"
   exit 1
fi

if [ ! -f "${opt_calls}" ]; then
   echo "# ${opt_calls} file not found!"
   echo "${usage}"
   exit 1
fi

callfile=$(basename ${opt_calls})

# check for VCF reference calls to be compared to (GIAB)
if [ -z "${opt_gold}" ]; then
	echo "# no VCF reference provided"
	echo "${usage}"
	exit 1
fi

if [ ! -f "${opt_gold}" ]; then
    echo "# ${opt_gold} file not found!";
    exit 1
fi

goldfile=$(basename ${opt_gold})

# check for FASTA reference used for mapping
if [ -z "${opt_ref}" ]; then
	echo "# no fasta reference provided"
	echo "${usage}"
	exit 1
fi

if [ ! -f "${opt_ref}" ]; then
    echo "# ${opt_ref} file not found!";
    exit 1
fi

reffile=$(basename ${opt_ref})

# get file prefix and sample_name
pfx=$(basename ${opt_calls%.vcf})

workdir=$PWD

outfolder=${opt_name:-"happy_${pfx}_results"}
input_dir="${workdir}/${outfolder}/input"
output_dir="${workdir}/${outfolder}"

# create local folders for inputs and results
# docker does not work with non-local paths
mkdir -p ${workdir}/${outfolder}/{input,logs}

#############################################
# copy reference and mappings               #
#############################################

echo "# REM: docker only works with local input files"

# copy if not yet there
if [ ! -f "${input_dir}/${callfile}" ]; then
echo "# copying test VCF file (and index if present)"
cp ${opt_calls}* ${input_dir}/
fi

# copy if not yet there
if [ ! -f "${input_dir}/${goldfile}" ]; then
echo "# copying VCF reference (goldstandard) file (and index if present)"
cp ${opt_gold}* ${input_dir}/
fi

# copy if not yet there
if [ ! -f "${input_dir}/${reffile}.fai" ]; then
echo "# copying reference FASTA file"
cp ${opt_ref} ${input_dir}/${reffile} && \
  samtools faidx ${input_dir}/${reffile}
fi

#################################
# run google deepvariant        #
#################################

# docker options
engine="vcfeval"
filter="--pass-only"

echo "# running Deepvariant with docker command:"

IFS='' read -r -d '' CMD <<EOF
time docker run \
  --rm \
  -it \
  -u "$(id -u):$(id -g)" \
  -v "${input_dir}":"/input" \
  -v "${output_dir}:/output" \
  ${img} \
  /opt/hap.py/bin/hap.py \
  /input/${goldfile} \
  /input/${callfile} \
  -r /input/${reffile} \
  -o /output/happy.${pfx}.output \
  --engine=${engine} \
  --threads ${nthr} \
  --logfile /output/logs/happy.${pfx}.output \
  ${filter} 
EOF

echo "# ${CMD}"
echo
eval ${CMD} > ${output_dir}/happy.${pfx}.output.txt 2>&1 

exit 0

######################################################
################### Command info #####################
######################################################

usage: Haplotype Comparison [-h] [-v] [-r REF] [-o REPORTS_PREFIX]
                            [--scratch-prefix SCRATCH_PREFIX] [--keep-scratch]
                            [-t {xcmp,ga4gh}] [-f FP_BEDFILE]
                            [--stratification STRAT_TSV]
                            [--stratification-region STRAT_REGIONS]
                            [--stratification-fixchr] [-V] [-X]
                            [--no-write-counts] [--output-vtc]
                            [--preserve-info] [--roc ROC] [--no-roc]
                            [--roc-regions ROC_REGIONS]
                            [--roc-filter ROC_FILTER] [--roc-delta ROC_DELTA]
                            [--ci-alpha CI_ALPHA] [--no-json]
                            [--location LOCATIONS] [--pass-only]
                            [--filters-only FILTERS_ONLY] [-R REGIONS_BEDFILE]
                            [-T TARGETS_BEDFILE] [-L] [--no-leftshift]
                            [--decompose] [-D] [--bcftools-norm] [--fixchr]
                            [--no-fixchr] [--bcf] [--somatic]
                            [--set-gt {half,hemi,het,hom,first}]
                            [--filter-nonref] [--convert-gvcf-truth]
                            [--convert-gvcf-query]
                            [--gender {male,female,auto,none}]
                            [--preprocess-truth] [--usefiltered-truth]
                            [--preprocessing-window-size PREPROCESS_WINDOW]
                            [--adjust-conf-regions] [--no-adjust-conf-regions]
                            [--unhappy] [-w WINDOW]
                            [--xcmp-enumeration-threshold MAX_ENUM]
                            [--xcmp-expand-hapblocks HB_EXPAND]
                            [--threads THREADS]
                            [--engine {xcmp,vcfeval,scmp-somatic,scmp-distance}]
                            [--engine-vcfeval-path ENGINE_VCFEVAL]
                            [--engine-vcfeval-template ENGINE_VCFEVAL_TEMPLATE]
                            [--scmp-distance ENGINE_SCMP_DISTANCE]
                            [--lose-match-distance ENGINE_SCMP_DISTANCE]
                            [--logfile LOGFILE] [--verbose | --quiet]
                            [_vcfs [_vcfs ...]]

positional arguments:
  _vcfs                 Two VCF files.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Show version number and exit.
  -r REF, --reference REF
                        Specify a reference file.
  -o REPORTS_PREFIX, --report-prefix REPORTS_PREFIX
                        Filename prefix for report output.
  --scratch-prefix SCRATCH_PREFIX
                        Directory for scratch files.
  --keep-scratch        Filename prefix for scratch report output.
  -t {xcmp,ga4gh}, --type {xcmp,ga4gh}
                        Annotation format in input VCF file.
  -f FP_BEDFILE, --false-positives FP_BEDFILE
                        False positive / confident call regions (.bed or
                        .bed.gz). Calls outside these regions will be labelled
                        as UNK.
  --stratification STRAT_TSV
                        Stratification file list (TSV format -- first column
                        is region name, second column is file name).
  --stratification-region STRAT_REGIONS
                        Add single stratification region, e.g.
                        --stratification-region TEST:test.bed
  --stratification-fixchr
                        Add chr prefix to stratification files if necessary
  -V, --write-vcf       Write an annotated VCF.
  -X, --write-counts    Write advanced counts and metrics.
  --no-write-counts     Do not write advanced counts and metrics.
  --output-vtc          Write VTC field in the final VCF which gives the
                        counts each position has contributed to.
  --preserve-info       When using XCMP, preserve and merge the INFO fields in
                        truth and query. Useful for ROC computation.
  --roc ROC             Select a feature to produce a ROC on (INFO feature,
                        QUAL, GQX, ...).
  --no-roc              Disable ROC computation and only output summary
                        statistics for more concise output.
  --roc-regions ROC_REGIONS
                        Select a list of regions to compute ROCs in. By
                        default, only the '*' region will produce ROC output
                        (aggregate variant counts).
  --roc-filter ROC_FILTER
                        Select a filter to ignore when making ROCs.
  --roc-delta ROC_DELTA
                        Minimum spacing between ROC QQ levels.
  --ci-alpha CI_ALPHA   Confidence level for Jeffrey's CI for recall,
                        precision and fraction of non-assessed calls.
  --no-json             Disable JSON file output.
  --location LOCATIONS, -l LOCATIONS
                        Comma-separated list of locations [use naming after
                        preprocessing], when not specified will use whole VCF.
  --pass-only           Keep only PASS variants.
  --filters-only FILTERS_ONLY
                        Specify a comma-separated list of filters to apply (by
                        default all filters are ignored / passed on.
  -R REGIONS_BEDFILE, --restrict-regions REGIONS_BEDFILE
                        Restrict analysis to given (sparse) regions (using -R
                        in bcftools).
  -T TARGETS_BEDFILE, --target-regions TARGETS_BEDFILE
                        Restrict analysis to given (dense) regions (using -T
                        in bcftools).
  -L, --leftshift       Left-shift variants safely.
  --no-leftshift        Do not left-shift variants safely.
  --decompose           Decompose variants into primitives. This results in
                        more granular counts.
  -D, --no-decompose    Do not decompose variants into primitives.
  --bcftools-norm       Enable preprocessing through bcftools norm -c x -D
                        (requires external preprocessing to be switched on).
  --fixchr              Add chr prefix to VCF records where necessary
                        (default: auto, attempt to match reference).
  --no-fixchr           Do not add chr prefix to VCF records (default: auto,
                        attempt to match reference).
  --bcf                 Use BCF internally. This is the default when the input
                        file is in BCF format already. Using BCF can speed up
                        temp file access, but may fail for VCF files that have
                        broken headers or records that don't comply with the
                        header.
  --somatic             Assume the input file is a somatic call file and
                        squash all columns into one, putting all FORMATs into
                        INFO + use half genotypes (see also --set-gt). This
                        will replace all sample columns and replace them with
                        a single one.
  --set-gt {half,hemi,het,hom,first}
                        This is used to treat Strelka somatic files Possible
                        values for this parameter: half / hemi / het / hom /
                        half to assign one of the following genotypes to the
                        resulting sample: 1 | 0/1 | 1/1 | ./1. This will
                        replace all sample columns and replace them with a
                        single one.
  --filter-nonref       Remove any variants genotyped as <NON_REF>.
  --convert-gvcf-truth  Convert the truth set from genome VCF format to a VCF
                        before processing.
  --convert-gvcf-query  Convert the query set from genome VCF format to a VCF
                        before processing.
  --gender {male,female,auto,none}
                        Specify sex. This determines how haploid calls on chrX
                        get treated: for male samples, all non-ref calls (in
                        the truthset only when running through hap.py) are
                        given a 1/1 genotype.
  --preprocess-truth    Preprocess truth file with same settings as query
                        (default is to accept truth in original format).
  --usefiltered-truth   Use filtered variant calls in truth file (by default,
                        only PASS calls in the truth file are used)
  --preprocessing-window-size PREPROCESS_WINDOW
                        Preprocessing window size (variants further apart than
                        that size are not expected to interfere).
  --adjust-conf-regions
                        Adjust confident regions to include variant locations.
                        Note this will only include variants that are included
                        in the CONF regions already when viewing with
                        bcftools; this option only makes sure insertions are
                        padded correctly in the CONF regions (to capture
                        these, both the base before and after must be
                        contained in the bed file).
  --no-adjust-conf-regions
                        Do not adjust confident regions for insertions.
  --unhappy, --no-haplotype-comparison
                        Disable haplotype comparison (only count direct GT
                        matches as TP).
  -w WINDOW, --window-size WINDOW
                        Minimum distance between variants such that they fall
                        into the same superlocus.
  --xcmp-enumeration-threshold MAX_ENUM
                        Enumeration threshold / maximum number of sequences to
                        enumerate per block.
  --xcmp-expand-hapblocks HB_EXPAND
                        Expand haplotype blocks by this many basepairs left
                        and right.
  --threads THREADS     Number of threads to use.
  --engine {xcmp,vcfeval,scmp-somatic,scmp-distance}
                        Comparison engine to use.
  --engine-vcfeval-path ENGINE_VCFEVAL
                        This parameter should give the path to the "rtg"
                        executable. The default is
                        /opt/hap.py/lib/python27/Haplo/../../../libexec/rtg-
                        tools-install/rtg
  --engine-vcfeval-template ENGINE_VCFEVAL_TEMPLATE
                        Vcfeval needs the reference sequence formatted in its
                        own file format (SDF -- run rtg format -o ref.SDF
                        ref.fa). You can specify this here to save time when
                        running hap.py with vcfeval. If no SDF folder is
                        specified, hap.py will create a temporary one.
  --scmp-distance ENGINE_SCMP_DISTANCE
                        For distance-based matching (vcfeval and scmp), this
                        is the distance between variants to use.
  --lose-match-distance ENGINE_SCMP_DISTANCE
                        For distance-based matching (vcfeval and scmp), this
                        is the distance between variants to use.
  --logfile LOGFILE     Write logging information into file rather than to
                        stderr
  --verbose             Raise logging level from warning to info.
  --quiet               Set logging level to output errors only.