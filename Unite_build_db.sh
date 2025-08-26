#!/usr/bin/env bash

# ============================================================================
# Unite_build_db.sh - build Unite database files for QIIME2
# 
# This script:
# 1. Downloads UNITE reference data using qiime2 rescript
# 2. Builds raw versions of the data
# 3. Build classifiers using qiime2
# Usage: ./Unite_build_db.sh
# SP@NC - email nucleomics@vib.be
# ============================================================================

# Script metadata
SCRIPT_VERSION="1.0.0"
SCRIPT_AUTHOR="SP@NC (+AI)"
REVISION_DATE="August 26, 2025"

# requires the following conda env for qiime2
QIIME_ENV="qiime2-amplicon-2025.7"
echo "Activating conda environment: $QIIME_ENV"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$QIIME_ENV"

CLUSTER_IDS=(97 99)
UNITE_DIR="$(pwd)"
UNITE_VERSION="2025-02-19"
TAXON_GROUP="fungi"
TMP_DIR="$PWD/tmp"
mkdir -p "${TMP_DIR}"

{
for CLUSTER_ID in "${CLUSTER_IDS[@]}"; do
  TAXONOMY_QZA="${UNITE_DIR}/unite-taxonomy-${TAXON_GROUP}_${CLUSTER_ID}pc_${UNITE_VERSION}.qza"
  SEQUENCES_QZA="${UNITE_DIR}/unite-sequences-${TAXON_GROUP}_${CLUSTER_ID}pc_${UNITE_VERSION}.qza"
  CLASSIFIER_QZA="${UNITE_DIR}/unite-classifier-${TAXON_GROUP}_${CLUSTER_ID}pc_${UNITE_VERSION}.qza"
  TAXONOMY_TSV="${UNITE_DIR}/unite-taxonomy-${TAXON_GROUP}_${CLUSTER_ID}pc_${UNITE_VERSION}.tsv"
  TRAINSET_FA="${UNITE_DIR}/unite-trainset-${TAXON_GROUP}_${CLUSTER_ID}pc_${UNITE_VERSION}.fa.gz"

  # Download data if not present
  if [[ ! -f "$TAXONOMY_QZA" || ! -f "$SEQUENCES_QZA" ]]; then
    echo "Downloading UNITE data for cluster ID $CLUSTER_ID..."
    qiime rescript get-unite-data \
      --p-version "$UNITE_VERSION" \
      --p-taxon-group "$TAXON_GROUP" \
      --p-cluster-id "$CLUSTER_ID" \
      --o-taxonomy "$TAXONOMY_QZA" \
      --o-sequences "$SEQUENCES_QZA"
  else
    echo "UNITE data for cluster ID $CLUSTER_ID already exists. Skipping download."
  fi

  # Export sequences (FASTA) if not present
  if [[ ! -f "$TRAINSET_FA" ]]; then
    echo "Exporting sequences for cluster ID $CLUSTER_ID..."
    qiime tools export \
      --input-path "$SEQUENCES_QZA" \
      --output-path "${TMP_DIR}" && \
      bgzip -c ${TMP_DIR}/dna-sequences.fasta >"$TRAINSET_FA"
  else
    echo "Trainset FASTA for cluster ID $CLUSTER_ID already exists. Skipping export."
  fi

  # Export taxonomy (TSV) if not present
  if [[ ! -f "$TAXONOMY_TSV" ]]; then
    echo "Exporting taxonomy for cluster ID $CLUSTER_ID..."
    qiime tools export \
      --input-path "$TAXONOMY_QZA" \
      --output-path "${TMP_DIR}" && \
      mv ${TMP_DIR}/taxonomy.tsv "$TAXONOMY_TSV"
  else
    echo "Taxonomy TSV for cluster ID $CLUSTER_ID already exists. Skipping export."
  fi

  # Train classifier if not present
  if [[ ! -f "$CLASSIFIER_QZA" ]]; then
    echo "Training classifier for cluster ID $CLUSTER_ID..."
    qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads "$SEQUENCES_QZA" \
      --i-reference-taxonomy "$TAXONOMY_QZA" \
      --o-classifier "$CLASSIFIER_QZA"
  else
    echo "Classifier for cluster ID $CLUSTER_ID already exists. Skipping training."
  fi
done
}

if [[ $? -eq 0 ]]; then
  rm -rf "${TMP_DIR}"
else
  echo "Error occurred during processing. Temporary files retained in ${TMP_DIR}." >&2
fi

exit 0

# qiime rescript get-unite-data --help
# Parameters:
#   --p-version TEXT Choices('2025-02-19', '2024-04-04', '2023-07-18',
#     '2022-10-16', '2021-05-10', '2020-02-20')
#                           UNITE version to download.   [default: '2025-02-19']
#   --p-taxon-group TEXT Choices('fungi', 'eukaryotes')
#                           Download a database with only 'fungi' or including
#                           all 'eukaryotes'.            [default: 'eukaryotes']
#   --p-cluster-id TEXT Choices('99', '97', 'dynamic')
#                           Percent similarity at which sequences in the of
#                           database were clustered.             [default: '99']
#   --p-singletons / --p-no-singletons
#                           Include singleton clusters in the database.
#                                                               [default: False]