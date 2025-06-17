#!/bin/bash

###############################################################################
# Script: Eukaryome_classify_ITS.sh
# Description: Classify OTU representative sequences with the Eukaryome reference DB
#
# Usage: ./Eukaryome_classify_ITS.sh -i INPUT.fa [options]
# SP@NC - email nucleomics@vib.be
# ============================================================================

SCRIPT_VERSION="1.0"
SCRIPT_AUTHOR="SP@NC"
REVISION_DATE="2025-06-17"

# Default values
INPUT_FASTA=""
EUKARYOME_VERSION="v1.9.4"
TAXON_GROUP="fungi"  # Keeping this as fungi as requested
QIIME_ENV="qiime2-amplicon-2025.4"
EUKARYOME_DIR="/data/biodata/qiime2_ITS"
OUTPUT_DIR="$(pwd)"
RESULTS_FOLDER="classification_results_Eukaryome"
THREADS="4"  # Default number of threads
BATCH_SIZE="auto"  # Default batch size for classification

# Display usage information
usage() {
  echo "Eukaryome_classify_ITS.sh v${SCRIPT_VERSION} - ${REVISION_DATE}"
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
  echo "  -v VERSION  Eukaryome version (default: ${EUKARYOME_VERSION})"
  echo "  -t TAXON    Taxon group (default: ${TAXON_GROUP})"
  echo "  -e ENV      Conda environment name (default: ${QIIME_ENV})"
  echo "  -d DIR      Eukaryome directory (default: ${EUKARYOME_DIR})"
  echo "  -p NUM      Number of threads to use (default: ${THREADS})"
  echo "  -b SIZE     Batch size for classification (default: ${BATCH_SIZE}, use a number for manual sizing)"
  echo "  -h          Display this help message"
  echo
  exit 1
}

# Parse command line arguments
while getopts "i:o:r:v:t:e:d:p:b:h" opt; do
  case $opt in
    i) INPUT_FASTA="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    r) RESULTS_FOLDER="$OPTARG" ;;
    v) EUKARYOME_VERSION="$OPTARG" ;;
    t) TAXON_GROUP="$OPTARG" ;;
    e) QIIME_ENV="$OPTARG" ;;
    d) EUKARYOME_DIR="$OPTARG" ;;
    p) THREADS="$OPTARG" ;;
    b) BATCH_SIZE="$OPTARG" ;;
    h) usage ;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
    :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# Check for required input file
if [ -z "$INPUT_FASTA" ]; then
  echo "Error: Input fasta file is required" >&2
  usage
fi

# Check if input file exists
if [ ! -f "$INPUT_FASTA" ]; then
  echo "Error: Input file not found: $INPUT_FASTA" >&2
  exit 1
fi

# Create full path to results directory
RESULTS_DIR="${OUTPUT_DIR}/${RESULTS_FOLDER}"
mkdir -p "$RESULTS_DIR" || { echo "Error: Could not create results directory"; exit 1; }

# Eukaryome classifier filename
EUKARYOME_CLASSIFIER="General_EUK_longread_${EUKARYOME_VERSION}_classifier.qza"

# Full path to the classifier
CLASSIFIER_PATH="${EUKARYOME_DIR}/${EUKARYOME_CLASSIFIER}"

# Check if classifier exists
if [ ! -f "$CLASSIFIER_PATH" ]; then
  echo "Error: Eukaryome classifier not found at: $CLASSIFIER_PATH" >&2
  exit 1
fi

# Log file
LOG_FILE="${RESULTS_DIR}/eukaryome_classification.log"

# Start logging
echo "Starting Eukaryome classification at $(date)" | tee -a "$LOG_FILE"
echo "Input file: $INPUT_FASTA" | tee -a "$LOG_FILE"
echo "Output directory: $RESULTS_DIR" | tee -a "$LOG_FILE"
echo "Eukaryome version: $EUKARYOME_VERSION" | tee -a "$LOG_FILE"
echo "Classifier: $CLASSIFIER_PATH" | tee -a "$LOG_FILE"

# Function to check if input is in FASTA format
check_fasta_format() {
  local file="$1"
  
  # Handle gzipped files
  if [[ "$file" == *.gz ]]; then
    if ! zcat "$file" | head -n 1 | grep -q '^>'; then
      echo "Error: Input file is not in FASTA format" | tee -a "$LOG_FILE" >&2
      exit 1
    fi
  else
    if ! head -n 1 "$file" | grep -q '^>'; then
      echo "Error: Input file is not in FASTA format" | tee -a "$LOG_FILE" >&2
      exit 1
    fi
  fi
}

# Save classifier information
echo "Saving Eukaryome classifier information..." | tee -a "$LOG_FILE"
INFO_FILE="${RESULTS_DIR}/classification_info.txt"
cat > "$INFO_FILE" << EOT
# Eukaryome Classifier Information
# Generated on $(date)
# by script: $(basename "$0") version $SCRIPT_VERSION

Classifier: $EUKARYOME_CLASSIFIER
Version: $EUKARYOME_VERSION
Location: $CLASSIFIER_PATH
Input file: $INPUT_FASTA
Output directory: $RESULTS_DIR
Threads used: $THREADS
QIIME environment: $QIIME_ENV
Batch size: $BATCH_SIZE
Taxon group: $TAXON_GROUP
EOT

# Activate conda environment
echo "Activating QIIME2 conda environment: $QIIME_ENV" | tee -a "$LOG_FILE"
eval "$(conda shell.bash hook)"
if ! conda activate "$QIIME_ENV"; then
  echo "Error: Failed to activate conda environment: $QIIME_ENV" | tee -a "$LOG_FILE" >&2
  exit 1
fi

# Check and process input file
check_fasta_format "$INPUT_FASTA"

# Determine if input is gzipped
if [[ "$INPUT_FASTA" == *.gz ]]; then
  echo "Decompressing gzipped input file..." | tee -a "$LOG_FILE"
  TEMP_FASTA="${RESULTS_DIR}/temp_input.fa"
  zcat "$INPUT_FASTA" > "$TEMP_FASTA"
  FASTA_TO_IMPORT="$TEMP_FASTA"
else
  FASTA_TO_IMPORT="$INPUT_FASTA"
fi

# Import FASTA to QIIME2 format
echo "Importing FASTA sequences to QIIME2 format..." | tee -a "$LOG_FILE"
SEQUENCES_QZA="${RESULTS_DIR}/sequences.qza"

qiime tools import \
  --input-path "$FASTA_TO_IMPORT" \
  --output-path "$SEQUENCES_QZA" \
  --type 'FeatureData[Sequence]'

# Clean up temporary files if created
if [[ -n "$TEMP_FASTA" && -f "$TEMP_FASTA" ]]; then
  rm "$TEMP_FASTA"
fi

# Set batch size parameter if not auto
BATCH_PARAM=""
if [[ "$BATCH_SIZE" != "auto" ]]; then
  BATCH_PARAM="--p-batch-size $BATCH_SIZE"
fi

# Run classification
echo "Running taxonomic classification with Eukaryome classifier..." | tee -a "$LOG_FILE"
echo "Started at $(date)" | tee -a "$LOG_FILE"

qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER_PATH" \
  --i-reads "$SEQUENCES_QZA" \
  --o-classification "${RESULTS_DIR}/taxonomy.qza" \
  --p-n-jobs "$THREADS" \
  $BATCH_PARAM

echo "Classification completed at $(date)" | tee -a "$LOG_FILE"

# Create visualization of the taxonomy
echo "Creating taxonomy visualization..." | tee -a "$LOG_FILE"
qiime metadata tabulate \
  --m-input-file "${RESULTS_DIR}/taxonomy.qza" \
  --o-visualization "${RESULTS_DIR}/taxonomy.qzv"

# Export taxonomy to TSV format
echo "Exporting taxonomy to TSV format..." | tee -a "$LOG_FILE"
qiime tools export \
  --input-path "${RESULTS_DIR}/taxonomy.qza" \
  --output-path "${RESULTS_DIR}/taxonomy_export"

# Move and rename the exported taxonomy file
if [ -f "${RESULTS_DIR}/taxonomy_export/taxonomy.tsv" ]; then
  mv "${RESULTS_DIR}/taxonomy_export/taxonomy.tsv" "${RESULTS_DIR}/classified_OTUs_LULU.tsv"
  # Clean up the temporary export directory
  rm -rf "${RESULTS_DIR}/taxonomy_export"
  echo "Exported taxonomy to: ${RESULTS_DIR}/classified_OTUs_LULU.tsv" | tee -a "$LOG_FILE"
else
  echo "Warning: Could not find exported taxonomy file" | tee -a "$LOG_FILE" >&2
fi

# Deactivate conda environment
conda deactivate

echo "All classification tasks completed successfully at $(date)" | tee -a "$LOG_FILE"
echo "Results saved to: $RESULTS_DIR" | tee -a "$LOG_FILE"