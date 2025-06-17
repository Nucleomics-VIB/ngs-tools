#!/usr/bin/env bash

# ============================================================================
# Unite_classify_ITS.sh - Taxonomic classification of ITS sequences using QIIME2
# 
# This script:
# 1. Downloads UNITE reference data (if needed)
# 2. Creates a classifier (if needed)
# 3. Classifies OTU sequences against the UNITE database
# 
# Usage: ./Unite_classify_ITS.sh -i INPUT.fa [options]
# SP@NC - email nucleomics@vib.be
# ============================================================================

# Script metadata
SCRIPT_VERSION="1.0.1"
SCRIPT_AUTHOR="SP@NC (+AI)"
REVISION_DATE="June 17, 2025"

# Exit on error
set -e

# Function to echo commands before executing them
echo_cmd() {
  echo "# $@"
  "$@"
}

# Default values
INPUT_FASTA=""
UNITE_VERSION="2025-02-19"
TAXON_GROUP="fungi"
CLUSTER_ID="99"
QIIME_ENV="qiime2-amplicon-2025.4"
UNITE_DIR="/data/biodata/qiime2_ITS"
OUTPUT_DIR="$(pwd)"
RESULTS_FOLDER="classification_results"
THREADS="4"  # Default number of threads
BATCH_SIZE="auto"  # Default batch size for classification

# Display usage information
usage() {
  echo "Unite_classify_ITS.sh v${SCRIPT_VERSION} - ${REVISION_DATE}"
  echo "Author: ${SCRIPT_AUTHOR}"
  echo
  echo "Usage: $0 -i INPUT.fa [options]"
  echo
  echo "Required arguments:"
  echo "  -i FILE     Input fasta file with OTU sequences (can be gzipped)"
  echo
  echo "Optional arguments:"
  echo "  -o DIR      Output directory (default: current directory)"
  echo "  -r NAME     Results folder name (default: ${RESULTS_FOLDER})"
  echo "  -v VERSION  UNITE version (default: ${UNITE_VERSION})"
  echo "  -t TAXON    Taxon group (default: ${TAXON_GROUP})"
  echo "  -c ID       Cluster ID (default: ${CLUSTER_ID})"
  echo "  -e ENV      Conda environment name (default: ${QIIME_ENV})"
  echo "  -d DIR      UNITE directory (default: ${UNITE_DIR})"
  echo "  -p NUM      Number of threads to use (default: ${THREADS})"
  echo "  -b SIZE     Batch size for classification (default: ${BATCH_SIZE}, use a number for manual sizing)"
  echo "  -h          Display this help message"
  echo
  exit 1
}

# Parse command-line options
while getopts "i:o:r:v:t:c:e:d:p:b:h" opt; do
  case $opt in
    i) INPUT_FASTA="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    r) RESULTS_FOLDER="$OPTARG" ;;
    v) UNITE_VERSION="$OPTARG" ;;
    t) TAXON_GROUP="$OPTARG" ;;
    c) CLUSTER_ID="$OPTARG" ;;
    e) QIIME_ENV="$OPTARG" ;;
    d) UNITE_DIR="$OPTARG" ;;
    p) THREADS="$OPTARG" ;;
    b) BATCH_SIZE="$OPTARG" ;;
    h) usage ;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
    :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# Check if input file was provided
if [ -z "$INPUT_FASTA" ]; then
  echo "Error: Input fasta file (-i) is required."
  usage
fi

# Check if input file exists
if [ ! -f "$INPUT_FASTA" ]; then
  echo "Error: Input file '$INPUT_FASTA' does not exist."
  exit 1
fi

# Function to prepare and get classifier
get_classifier() {
  local unite_dir="$1"
  local unite_version="$2"
  local taxon_group="$3"
  local cluster_id="$4"
  local threads="$5"
  
  # Create UNITE directory if it doesn't exist (silently)
  mkdir -p "${unite_dir}" >&2
  
  # Prepare classifier name with version information
  local classifier_name="unite-classifier-${taxon_group}_${cluster_id}pc_${unite_version}.qza"
  local taxonomy_name="unite-taxonomy-${taxon_group}_${cluster_id}pc_${unite_version}.qza"
  local sequences_name="unite-sequences-${taxon_group}_${cluster_id}pc_${unite_version}.qza"
  
  # Check if classifier already exists
  if [ ! -f "${unite_dir}/${classifier_name}" ]; then
    echo "Classifier not found. Building new classifier..." >&2
    
    # Download UNITE reference data
    echo "Downloading UNITE reference data (version: ${unite_version})..." >&2
    echo "# qiime rescript get-unite-data --p-version ${unite_version} --p-taxon-group ${taxon_group} --p-cluster-id ${cluster_id} --o-taxonomy ${unite_dir}/${taxonomy_name} --o-sequences ${unite_dir}/${sequences_name}" >&2
    qiime rescript get-unite-data \
      --p-version "${unite_version}" \
      --p-taxon-group "${taxon_group}" \
      --p-cluster-id "${cluster_id}" \
      --o-taxonomy "${unite_dir}/${taxonomy_name}" \
      --o-sequences "${unite_dir}/${sequences_name}" >&2
    
    # Build the classifier (with multithreading)
    echo "Building naive Bayes classifier (using ${threads} threads)..." >&2
    echo "# qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ${unite_dir}/${sequences_name} --i-reference-taxonomy ${unite_dir}/${taxonomy_name} --o-classifier ${unite_dir}/${classifier_name} --p-n-jobs ${threads}" >&2
    qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads "${unite_dir}/${sequences_name}" \
      --i-reference-taxonomy "${unite_dir}/${taxonomy_name}" \
      --o-classifier "${unite_dir}/${classifier_name}" \
      --p-n-jobs "${threads}" >&2
    
    echo "Classifier built and saved to ${unite_dir}/${classifier_name}" >&2
  else
    echo "Using existing classifier: ${unite_dir}/${classifier_name}" >&2
  fi
  
  # Return ONLY the full path to the classifier
  echo "${unite_dir}/${classifier_name}"
}

# Function to process a fasta file
process_fasta() {
  local input_file="$1"
  local classifier="$2"
  local results_dir="$3"
  local threads="$4"
  local batch_size="$5"
  
  local base_name=$(basename "$input_file")
  local input_name_no_ext="${base_name%.*}"
  
  # Check if input file is gzipped and decompress during copy if needed
  local copied_input=""
  if [[ "$input_file" == *.gz ]]; then
    echo "Detected gzipped input file"
    # Remove .gz extension to get the base name
    local uncompressed_name="${input_name_no_ext}"
    
    echo "Decompressing input file to results directory..."
    echo "# gunzip -c ${input_file} > ${results_dir}/${uncompressed_name}"
    gunzip -c "$input_file" > "${results_dir}/${uncompressed_name}"
    copied_input="${results_dir}/${uncompressed_name}"
    
    echo "Created uncompressed version: ${copied_input}"
  else
    # Copy the input file to the results directory for portability
    echo "Copying input file to results directory..."
    echo "# cp ${input_file} ${results_dir}/"
    cp "$input_file" "${results_dir}/"
    copied_input="${results_dir}/${base_name}"
  fi
  
  # Import sequences - use the copied/uncompressed file
  echo "Importing sequences from input file..."
  echo "# qiime tools import --type 'FeatureData[Sequence]' --input-path ${copied_input} --output-path ${results_dir}/sequences.qza"
  qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path "${copied_input}" \
    --output-path "${results_dir}/sequences.qza"
  
  # Classify sequences with generic output name (with multithreading)
  echo "Classifying sequences against UNITE database (using ${threads} threads)..."
  echo "# qiime feature-classifier classify-sklearn --i-classifier ${classifier} --i-reads ${results_dir}/sequences.qza --o-classification ${results_dir}/taxonomy.qza --p-n-jobs ${threads} --p-reads-per-batch ${batch_size} --p-pre-dispatch '2*n_jobs'"
  qiime feature-classifier classify-sklearn \
    --i-classifier "${classifier}" \
    --i-reads "${results_dir}/sequences.qza" \
    --o-classification "${results_dir}/taxonomy.qza" \
    --p-n-jobs "${threads}" \
    --p-reads-per-batch "${batch_size}" \
    --p-pre-dispatch "2*n_jobs"
  
  # Create a temporary directory for export
  local tmp_export_dir="/tmp/taxonomy-export-$$-${RANDOM}"
  
  # Export taxonomy results to temporary directory
  echo "Exporting taxonomy results..."
  echo "# qiime tools export --input-path ${results_dir}/taxonomy.qza --output-path ${tmp_export_dir}"
  qiime tools export \
    --input-path "${results_dir}/taxonomy.qza" \
    --output-path "${tmp_export_dir}"
  
  # Move only the taxonomy.tsv file to the results directory with the desired name
  echo "Moving taxonomy file to results directory..."
  echo "# mv ${tmp_export_dir}/taxonomy.tsv ${results_dir}/classified_${input_name_no_ext}.tsv"
  mv "${tmp_export_dir}/taxonomy.tsv" "${results_dir}/classified_${input_name_no_ext}.tsv"
  
  # Clean up the temporary directory
  echo "Cleaning up temporary files..."
  echo "# rm -rf ${tmp_export_dir}"
  rm -rf "${tmp_export_dir}"
  
  echo "Classification complete for ${base_name}"
  echo "Results saved as: ${results_dir}/classified_${input_name_no_ext}.tsv"
}

# Create output directory if it doesn't exist
echo "# mkdir -p ${OUTPUT_DIR}"
mkdir -p "$OUTPUT_DIR"

# Create results directory
RESULTS_DIR="${OUTPUT_DIR}/${RESULTS_FOLDER}"
echo "# mkdir -p ${RESULTS_DIR}"
mkdir -p "$RESULTS_DIR"
echo "All results will be saved to: ${RESULTS_DIR}"

# Create UNITE directory if it doesn't exist
echo "# mkdir -p ${UNITE_DIR}"
mkdir -p "$UNITE_DIR"

# Check if required conda environment exists
echo "# Checking for conda environment: ${QIIME_ENV}"
if ! conda env list | grep -q "${QIIME_ENV}"; then
  echo "Error: Required conda environment '${QIIME_ENV}' not found!"
  echo "Please create it with: conda create -n ${QIIME_ENV} -c qiime2 qiime2"
  exit 1
fi

# Activate QIIME2 environment
echo "Activating QIIME2 environment..."
echo "# conda activate ${QIIME_ENV}"
eval "$(conda shell.bash hook)"
conda activate "${QIIME_ENV}"

# Get or build the classifier
echo "# Getting classifier..."
CLASSIFIER_PATH=$(get_classifier "$UNITE_DIR" "$UNITE_VERSION" "$TAXON_GROUP" "$CLUSTER_ID" "$THREADS")

# Save a copy of the classifier information
echo "# Creating classification_info.txt"
echo "Classification parameters:" > "${RESULTS_DIR}/classification_info.txt"
echo "Date: $(date)" >> "${RESULTS_DIR}/classification_info.txt"
echo "Input file: ${INPUT_FASTA}" >> "${RESULTS_DIR}/classification_info.txt"
echo "UNITE version: ${UNITE_VERSION}" >> "${RESULTS_DIR}/classification_info.txt"
echo "Taxon group: ${TAXON_GROUP}" >> "${RESULTS_DIR}/classification_info.txt"
echo "Cluster ID: ${CLUSTER_ID}" >> "${RESULTS_DIR}/classification_info.txt"
echo "Threads: ${THREADS}" >> "${RESULTS_DIR}/classification_info.txt"
echo "Classifier: ${CLASSIFIER_PATH}" >> "${RESULTS_DIR}/classification_info.txt"

# Process the input fasta file with generic output names
echo "# Processing input fasta file..."
process_fasta "$INPUT_FASTA" "$CLASSIFIER_PATH" "$RESULTS_DIR" "$THREADS" "$BATCH_SIZE"

echo "ITS classification completed!"
echo "Results saved to: ${RESULTS_DIR}"