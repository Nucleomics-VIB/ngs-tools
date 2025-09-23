#!/bin/bash

# download_Refseq_NR.sh
# Script to download RefSeq non-redundant protein files from NCBI
# using aria2c with multithreading, from a configurable base URL.
# SP@NC 2025-09-23; v1.0 (+Perplexity)
#
# Script to download all RefSeq non-redundant protein files
# (<clade>.nonredundant_protein.<part#>.protein.faa.gz)
# from user-selectable RefSeq taxonomy subsets (plus 'complete') 
# using aria2c multithreading.
#
# Usage:
#   ./download_Refseq_NR.sh [-t TYPE]
#
# Allowed taxonomy TYPE values (no 'release*' folders):
# - archaea
# - bacteria
# - complete
# - fungi
# - invertebrate
# - microbial
# - mitochondrion
# - other
# - plant
# - plasmid
# - plastid
# - protozoa
# - vertebrate_mammalian
# - vertebrate_other
# - viral
#
# Requires: curl, grep, sed, awk, aria2c

set -e

# Default subset
TYPE="complete"

# Parse optional arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -t|--type)
      TYPE="$2"
      shift; shift
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: $0 [-t|--type TYPE]"
      exit 1
      ;;
  esac
done

# Allowed taxonomy subsets including 'complete'
ALLOWED_TYPES=(
  "archaea" "bacteria" "complete" "fungi" "invertebrate" "microbial" "mitochondrion"
  "other" "plant" "plasmid" "plastid" "protozoa" "vertebrate_mammalian" "vertebrate_other" "viral"
)

# Validate input - exclude anything starting with "release"
if [[ "$TYPE" == release* ]]; then
  echo "Error: subsets starting with 'release' are not allowed."
  exit 1
fi

if [[ ! " ${ALLOWED_TYPES[@]} " =~ " ${TYPE} " ]]; then
  echo "Error: type must be one of: ${ALLOWED_TYPES[*]}"
  exit 1
fi

# Mirror base URL
MIRROR_BASE="https://ftp.funet.fi/pub/mirrors/ftp.ncbi.nlm.nih.gov/refseq/release"

# Full base URL
BASE_URL="${MIRROR_BASE}/${TYPE}/"

echo "Downloading RefSeq nonredundant proteins from subset: $TYPE"
echo "Using base URL: $BASE_URL"

# Download directory listing HTML silently
curl -s "$BASE_URL" > index.html

# Extract matching files, clean, prepend base URL, output to urls.txt
grep -oP 'href="[^"]*nonredundant_protein\.\d+\.protein\.faa\.gz"' index.html | \
sed 's/href="//; s/"//' | \
awk -v base="$BASE_URL" '{print base $0}' > urls.txt

num_files=$(wc -l < urls.txt)
echo "Found $num_files protein files to download."

if [ "$num_files" -eq 0 ]; then
  echo "No matching protein files found in $BASE_URL. Exiting."
  exit 1
fi

# Perform download
aria2c -x 16 -s 16 -j 10 -i urls.txt
