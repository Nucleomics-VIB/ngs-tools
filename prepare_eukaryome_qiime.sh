#!/bin/bash

###############################################################################
# Script: prepare_eukaryome_qiime.sh
#
# Description:
#   This script automates the process of downloading a zipped EUKARYOME reference
#   database (or a similar resource), extracting the FASTA file, generating a
#   QIIME 2-compatible taxonomy mapping file, and creating a new FASTA file with
#   simplified headers (unique IDs only). It then imports both files into QIIME 2
#   artifacts using the qiime2-amplicon-2025.4 conda environment and
#   trains a Naive Bayes classifier.
#
# Usage:
#   ./prepare_eukaryome_qiime.sh [DOWNLOAD_LINK]
#
#   - If DOWNLOAD_LINK is not provided, the script defaults to the official
#     EUKARYOME v1.9.4 longread reference.
#
# Output:
#   - QIIME 2 sequence artifact (.qza)
#   - QIIME 2 taxonomy artifact (.qza)
#   - QIIME 2 trained classifier (.qza)
#
# Requirements:
#   - wget, unzip, awk, sed, grep, conda, QIIME 2 (qiime2-amplicon-2025.4)
###############################################################################

# Script metadata
SCRIPT_VERSION="1.0.0"
SCRIPT_AUTHOR="SP@NC (+AI)"
REVISION_DATE="June 17, 2025"

# Exit on error
set -e

# Default configuration
DEFAULT_ZIP_URL="https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_longread_v1.9.4.zip"
QIIME_ENV="qiime2-amplicon-2025.4"

# Global variables
ZIP_URL=""
ZIP_FILE=""
ORIG_FASTA=""
QIIME_FASTA=""
TAX_TSV=""
SEQ_QZA=""
TAX_QZA=""
CLASSIFIER_QZA=""

# Function to display script information
display_info() {
    echo "============================================================================"
    echo "prepare_eukaryome_qiime.sh v${SCRIPT_VERSION}"
    echo "Author: ${SCRIPT_AUTHOR}"
    echo "Last revised: ${REVISION_DATE}"
    echo "Contact: nucleomics@vib.be"
    echo "============================================================================"
}

# Function to initialize variables
init_variables() {
    # Use first argument as URL if provided, else default
    ZIP_URL="${1:-$DEFAULT_ZIP_URL}"
    
    # Extract the base filename from the URL for naming
    ZIP_FILE=$(basename "$ZIP_URL")
    
    # Set filenames based on the zip file name
    ORIG_FASTA="General_EUK_longread_v1.9.4.fasta"
    QIIME_FASTA="General_EUK_longread_v1.9.4_qiime.fasta"
    TAX_TSV="General_EUK_longread_v1.9.4_taxonomy.tsv"
    SEQ_QZA="General_EUK_longread_v1.9.4_sequences.qza"
    TAX_QZA="General_EUK_longread_v1.9.4_taxonomy.qza"
    CLASSIFIER_QZA="General_EUK_longread_v1.9.4_classifier.qza"
    
    echo "Download link: $ZIP_URL"
    echo "Zip file: $ZIP_FILE"
    echo "Assumed FASTA file inside zip: $ORIG_FASTA"
}

# Function to check prerequisites
check_prerequisites() {
    local required_commands=("wget" "unzip" "awk" "sed" "grep" "conda")
    local missing_commands=()
    
    echo "Checking prerequisites..."
    
    for cmd in "${required_commands[@]}"; do
        if ! command -v "$cmd" &> /dev/null; then
            missing_commands+=("$cmd")
        fi
    done
    
    if [ ${#missing_commands[@]} -gt 0 ]; then
        echo "ERROR: The following required commands are missing:"
        for cmd in "${missing_commands[@]}"; do
            echo "  - $cmd"
        done
        echo "Please install them before running this script."
        return 1
    fi
    
    # Check if QIIME 2 environment exists
    if ! conda env list | grep -q "$QIIME_ENV"; then
        echo "ERROR: QIIME 2 environment '$QIIME_ENV' not found."
        echo "Please create it first with: conda create -n $QIIME_ENV -c qiime2 qiime2"
        return 1
    fi
    
    echo "All prerequisites satisfied."
    return 0
}

# Function to download the zip file
download_zip() {
    if [ -f "$ZIP_FILE" ]; then
        echo "$ZIP_FILE already exists. Skipping download."
        return 0
    fi
    
    echo "Downloading reference database..."
    if wget -N "$ZIP_URL"; then
        echo "Download complete."
        return 0
    else
        echo "ERROR: Failed to download $ZIP_URL"
        return 1
    fi
}

# Function to extract the zip file
extract_zip() {
    if [ -f "$ORIG_FASTA" ]; then
        echo "$ORIG_FASTA already exists. Skipping unzip."
        return 0
    fi
    
    echo "Unzipping the downloaded file..."
    if unzip -o "$ZIP_FILE"; then
        echo "Unzip complete."
        if [ ! -f "$ORIG_FASTA" ]; then
            echo "ERROR: Expected FASTA file $ORIG_FASTA not found in zip."
            echo "Please check the zip contents and update the script if needed."
            return 1
        fi
        return 0
    else
        echo "ERROR: Failed to unzip $ZIP_FILE"
        return 1
    fi
}

# Function to create simplified FASTA for QIIME 2
create_qiime_fasta() {
    if [ -f "$QIIME_FASTA" ]; then
        echo "$QIIME_FASTA already exists. Skipping FASTA header simplification."
        return 0
    fi
    
    echo "Creating QIIME-compatible FASTA file with simplified headers..."
    awk '/^>/ {split($0, a, ";"); print a[1]} !/^>/' "$ORIG_FASTA" > "$QIIME_FASTA"
    
    if [ -f "$QIIME_FASTA" ]; then
        echo "QIIME FASTA creation complete."
        return 0
    else
        echo "ERROR: Failed to create $QIIME_FASTA"
        return 1
    fi
}

# Function to create taxonomy TSV file
create_taxonomy_tsv() {
    if [ -f "$TAX_TSV" ]; then
        echo "$TAX_TSV already exists. Skipping taxonomy TSV generation."
        return 0
    fi
    
    echo "Generating taxonomy TSV file..."
    grep '^>' "$ORIG_FASTA" | \
      sed 's/^>//' | \
      awk -F';' '{printf "%s\t", $1; for(i=2;i<=NF;i++) {printf "%s%s", $i, (i<NF?";":"")}; print ""}' \
      > "$TAX_TSV"
    
    if [ ! -f "$TAX_TSV" ]; then
        echo "ERROR: Failed to create taxonomy TSV file."
        return 1
    fi
    
    echo "Adding header to taxonomy TSV file..."
    sed -i '1iFeature ID\tTaxon' "$TAX_TSV"
    
    echo "Taxonomy TSV creation complete."
    return 0
}

# Function to activate QIIME 2 environment
activate_qiime() {
    echo "Activating QIIME 2 environment ($QIIME_ENV)..."
    source "$(conda info --base)/etc/profile.d/conda.sh"
    if conda activate "$QIIME_ENV"; then
        echo "QIIME 2 environment activated."
        return 0
    else
        echo "ERROR: Failed to activate QIIME 2 environment."
        return 1
    fi
}

# Function to import sequences into QIIME 2
import_sequences() {
    if [ -f "$SEQ_QZA" ]; then
        echo "$SEQ_QZA already exists. Skipping QIIME 2 sequence import."
        return 0
    fi
    
    echo "Importing sequences into QIIME 2 artifact..."
    if qiime tools import \
         --type 'FeatureData[Sequence]' \
         --input-path "$QIIME_FASTA" \
         --output-path "$SEQ_QZA"; then
        echo "Sequence import complete."
        return 0
    else
        echo "ERROR: Failed to import sequences into QIIME 2."
        return 1
    fi
}

# Function to import taxonomy into QIIME 2
import_taxonomy() {
    if [ -f "$TAX_QZA" ]; then
        echo "$TAX_QZA already exists. Skipping QIIME 2 taxonomy import."
        return 0
    fi
    
    echo "Importing taxonomy into QIIME 2 artifact..."
    if qiime tools import \
         --type 'FeatureData[Taxonomy]' \
         --input-format HeaderlessTSVTaxonomyFormat \
         --input-path "$TAX_TSV" \
         --output-path "$TAX_QZA"; then
        echo "Taxonomy import complete."
        return 0
    else
        echo "ERROR: Failed to import taxonomy into QIIME 2."
        return 1
    fi
}

# Function to train Naive Bayes classifier
train_classifier() {
    if [ -f "$CLASSIFIER_QZA" ]; then
        echo "$CLASSIFIER_QZA already exists. Skipping classifier training."
        return 0
    fi
    
    echo "Training Naive Bayes classifier (this step can take a long time)..."
    if qiime feature-classifier fit-classifier-naive-bayes \
         --i-reference-reads "$SEQ_QZA" \
         --i-reference-taxonomy "$TAX_QZA" \
         --o-classifier "$CLASSIFIER_QZA"; then
        echo "Classifier training complete."
        return 0
    else
        echo "ERROR: Failed to train Naive Bayes classifier."
        return 1
    fi
}

# Function to create documentation file
create_documentation() {
    local doc_file="General_EUK_longread_v1.9.4_processing_info.txt"
    local current_date=$(date "+%B %d, %Y")
    
    echo "Creating documentation file..."
    
    cat > "$doc_file" << EOF
# EUKARYOME QIIME 2 Processing Documentation

## Processing Information
- Date: $current_date
- Script Used: $(basename "$0")
- Processing Purpose: Preparation of EUKARYOME reference database for QIIME 2 analysis

## Input Source
- Database: EUKARYOME v1.9.4 longread reference
- Download URL: $ZIP_URL
- Original Format: Zipped FASTA file with taxonomic information embedded in headers

## Processing Steps
The script performed the following operations:
1. Downloaded the reference database zip file
2. Extracted the FASTA file
3. Created a simplified FASTA file with QIIME 2-compatible headers
4. Generated a taxonomy mapping file in TSV format
5. Imported sequences and taxonomy into QIIME 2 artifacts
6. Trained a Naive Bayes classifier

## Output Files
- QIIME 2 Sequence Artifact: $SEQ_QZA
- QIIME 2 Taxonomy Artifact: $TAX_QZA
- QIIME 2 Trained Classifier: $CLASSIFIER_QZA

## Notes
These files can be used with QIIME 2 for taxonomic classification of eukaryotic amplicon sequences.
The classifier is particularly suited for long-read sequencing data.
EOF

    if [ -f "$doc_file" ]; then
        echo "Documentation file created: $doc_file"
        return 0
    else
        echo "ERROR: Failed to create documentation file."
        return 1
    fi
}

# Main function to orchestrate the workflow
main() {
    local download_url="$1"
    
    # Initialize variables
    init_variables "$download_url"
    
    # Check prerequisites
    if ! check_prerequisites; then
        echo "Prerequisite check failed. Exiting."
        exit 1
    fi
    
    # Execute each step in sequence, stopping if any fails
    if ! download_zip; then
        echo "Download failed. Exiting."
        exit 1
    fi
    
    if ! extract_zip; then
        echo "Extraction failed. Exiting."
        exit 1
    fi
    
    if ! create_qiime_fasta; then
        echo "FASTA preparation failed. Exiting."
        exit 1
    fi
    
    if ! create_taxonomy_tsv; then
        echo "Taxonomy preparation failed. Exiting."
        exit 1
    fi
    
    if ! activate_qiime; then
        echo "QIIME 2 activation failed. Exiting."
        exit 1
    fi
    
    if ! import_sequences; then
        echo "Sequence import failed. Exiting."
        exit 1
    fi
    
    if ! import_taxonomy; then
        echo "Taxonomy import failed. Exiting."
        exit 1
    fi
    
    if ! train_classifier; then
        echo "Classifier training failed. Exiting."
        exit 1
    fi
    
    if ! create_documentation; then
        echo "Documentation creation failed. Exiting."
        exit 1
    fi
    
    echo "All steps complete."
    echo "QIIME 2 artifacts created (if not already present):"
    echo "  Sequences:  $SEQ_QZA"
    echo "  Taxonomy:   $TAX_QZA"
    echo "  Classifier: $CLASSIFIER_QZA"
    
    return 0
}

# Execute the main function with all arguments
main "$@"
