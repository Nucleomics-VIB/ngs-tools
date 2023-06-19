#!/bin/bash

# call variants from long reads and a reference
# using PEPPER Margin Deepvariant
#
# Requirements:
# your user is part of the docker group to avoid sudo
# run on a unix computer installed with
#   docker image: google/deepvariant:latest (now 1.5.0)
#   # get with: sudo docker pull google/deepvariant:1.5.0
#   samtools, minimap2,
#   reference genome fasta present and indexed
#   a minimum of 4 threads
#
# this script will write results in the current folder
#
# Stephane Plaisance (VIB-NC) 2022/06/01; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# version number
version="1.0, 2023_06_15"

usage='# Usage: run_Deepvariant.sh -b <Aligned, sorted, indexed BAM file> -r <reference fasta>
# script version '${version}'
# [optional: -n <sample name (default BAM prefix)>]
# [optional: -m <model (default to WGS)>]
# [optional: -t <threads|4>]
# [optional: -X <extra parameter string between quotes for the docker command 
# (eg --dry_run to get a list of commands without running them)>]'

while getopts "b:r:n:m:t:X:h" opt; do
  case $opt in
    b) opt_mappings=${OPTARG} ;;
    r) opt_ref=${OPTARG} ;;
    n) opt_name=${OPTARG} ;;
    m) opt_model=${OPTARG} ;;
    t) opt_thr=${OPTARG} ;;
    X) opt_extra=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2;
       exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; 
       exit 1 ;;
  esac
done

# set defaults
nthr=${opt_thr:-4}
extra=${opt_extra:-""}

# check if all dependencies are present
declare -a arr=("samtools")
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || \
  ( echo "# required ${prog} not found in PATH"
exit 1 )
done

# check for the docker image
res=$(docker images | grep google/deepvariant:latest | cut -d " " -f 1)
if [ ! ${res} == "google/deepvariant" ]; then
  echo "docker image not found"
  exit 1
fi

# test if minimal arguments were provided
if [ -z "${opt_mappings}" ]; then
   echo "# no mappings provided!"
   echo "${usage}"
   exit 1
fi

if [ ! -f "${opt_mappings}" ]; then
   echo "# ${opt_mappings} file not found!"
   echo "${usage}"
   exit 1
fi

mappingfile=$(basename ${opt_mappings})

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
pfx=$(basename ${opt_name:-${mappingfile%.bam}})

# other parameters or defaults
workdir=$PWD

out_pfx="DeepVariant_${pfx}_results"
input_dir="${workdir}/${out_pfx}/input"
output_dir="${workdir}/${out_pfx}"

# create local folders for inputs and results
# docker does not work with non-local paths
mkdir -p ${workdir}/${out_pfx}/{input,intermediate_files,logs}

#############################################
# copy reference and mappings               #
#############################################

echo "# REM: docker only works with local input files"

# copy if not yet there
if [ ! -f "${input_dir}/${reffile}.fai" ]; then
echo "# copying reference file and indexing"
cp ${opt_ref} ${input_dir}/${reffile} && \
  samtools faidx ${input_dir}/${reffile}
fi

# copy if not yet there
if [ ! -f "${input_dir}/${mappingfile}.bai" ]; then
echo "# copying mapping file and indexing"
cp ${opt_mappings} ${input_dir}/${mappingfile} && \
  samtools index -@4 ${input_dir}/${mappingfile}
fi

#################################
# run google deepvariant        #
#################################

# DeepVariant docker
BIN_VERSION="latest"
img="google/deepvariant:${BIN_VERSION}"

# define model
# WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMINA
model=${opt_model:-"WGS"}

echo "# running Deepvariant with docker command:"
echo ${reffile}
IFS='' read -r -d '' CMD <<EOF
time docker run \
  --rm \
  -u "$(id -u):$(id -g)" \
  -v "${input_dir}":"/input" \
  -v "${output_dir}":"/output" \
  ${img} \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=${model} \
  --ref=/input/${reffile} \
  --reads=/input/${mappingfile} \
  --sample_name=${pfx} \
  --output_vcf=/output/${pfx}.vcf.gz \
  --output_gvcf=/output/${pfx}.g.vcf.gz \
  --vcf_stats_report \
  --intermediate_results_dir=/output/intermediate_results_dir \
  --logging_dir=/output/logs \
  --num_shards=${nthr} \
  ${extra}
EOF

echo "# ${CMD}"
echo
eval ${CMD}

exit 0

######################################################
################### Command info #####################
######################################################

For more details, see:
https://github.com/google/deepvariant/blob/r1.5/docs/deepvariant-quick-start.md

flags:

/opt/deepvariant/bin/run_deepvariant.py:
  --call_variants_extra_args: A comma-separated list of flag_name=flag_value.
    "flag_name" has to be valid flags for call_variants.py. If the flag_value is
    boolean, it has to be flag_name=true or flag_name=false.
  --customized_model: Optional. A path to a model checkpoint to load for the
    `call_variants` step. If not set, the default for each --model_type will be
    used
  --[no]dry_run: Optional. If True, only prints out commands without executing
    them.
    (default: 'false')
  --intermediate_results_dir: Optional. If specified, this should be an existing
    directory that is visible insider docker, and will be used to to store
    intermediate outputs.
  --logging_dir: Optional. Directory where we should write log files for each
    stage and optionally runtime reports.
  --make_examples_extra_args: A comma-separated list of flag_name=flag_value.
    "flag_name" has to be valid flags for make_examples.py. If the flag_value is
    boolean, it has to be flag_name=true or flag_name=false.
  --model_type: <WGS|WES|PACBIO|ONT_R104|HYBRID_PACBIO_ILLUMINA>: Required. Type
    of model to use for variant calling. Set this flag to use the default model
    associated with each type, and it will set necessary flags corresponding to
    each model. If you want to use a customized model, add --customized_model
    flag in addition to this flag.
  --num_shards: Optional. Number of shards for make_examples step.
    (default: '1')
    (an integer)
  --output_gvcf: Optional. Path where we should write gVCF file.
  --output_vcf: Required. Path where we should write VCF file.
  --postprocess_variants_extra_args: A comma-separated list of
    flag_name=flag_value. "flag_name" has to be valid flags for
    postprocess_variants.py. If the flag_value is boolean, it has to be
    flag_name=true or flag_name=false.
  --reads: Required. Aligned, sorted, indexed BAM file containing the reads we
    want to call. Should be aligned to a reference genome compatible with --ref.
  --ref: Required. Genome reference to use. Must have an associated FAI index as
    well. Supports text or gzipped references. Should match the reference used
    to align the BAM file provided to --reads.
  --regions: Optional. Space-separated list of regions we want to process.
    Elements can be region literals (e.g., chr20:10-20) or paths to BED/BEDPE
    files.
  --[no]runtime_report: Output make_examples runtime metrics and create a visual
    runtime report using runtime_by_region_vis. Only works with --logging_dir.
    (default: 'false')
  --sample_name: Sample name to use instead of the sample name from the input
    reads BAM (SM tag in the header). This flag is used for both make_examples
    and postprocess_variants.
  --[no]use_hp_information: (Deprecated in v1.4.0) Optional. If True,
    corresponding flags will be set to properly use the HP information present
    in the BAM input.
  --[no]vcf_stats_report: Optional. Output a visual report (HTML) of statistics
    about the output VCF.
    (default: 'true')
  --[no]version: Optional. If true, print out version number and exit.

Try --helpfull to get a list of all flags.