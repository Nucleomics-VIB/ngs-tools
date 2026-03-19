#!/bin/bash

# script: filter_pe_fastp_sequential.sh
# author SP@NC (+AI)
# 2026-03-09; version 1.4 (Fixed path handling)
#
###############################################################################
# filter_pe_fastp_sequential.sh - Sequential Quality + Length (+Optional UMI)
#
# Purpose: TWO/THREE-STEP filtering on paired-end FASTQ reads (PE300):
#   [0.] UMI extraction (first N nt R1+R2 → read names) IF -u specified
#   1. Quality filtering (average Phred > Q30 for BOTH reads)
#   2. Length filtering (BOTH reads >= MIN_LEN post-trimming) on QC-passed only
#
# Key Features:
#   - OPTIONAL UMI: -u [N] (extract first N nt from read starts, DEFAULT N=8)
#   - SEQUENTIAL: UMI→QC→length (or QC→length if no UMI)
#   - Paired-end: If ONE read fails ANY criterion, ENTIRE PAIR rejected
#   - ROBUST PATHS: Works with files in folders (/path/to/*.fq.gz)
#   - Outputs: good/, failed_quality/, failed_length/ directories
#   - NO short reads (< MIN_LEN) in ANY final output set
#   - Multi-threaded fastp (THREADS=8 default)
#   - HTML/JSON reports for each step
#   - Auto cleanup of intermediate files
#
# Workflow (with -u):
# input ──(UMI:Nnt)───(QC:avg>Q30)──> good/_qc.fq.gz ──(Len>=MIN_LEN)──> good/_final.fq.gz
#        ────────────(QC fail)──────> failed_quality/          failed_length/
#
# Usage: $0 -i input_R1.fq.gz -I input_R2.fq.gz [-u [UMI_LEN]] [-l MIN_LEN] [-t THREADS]
# Example: $0 -i sample_R1.fq.gz -I sample_R2.fq.gz -l 250
# Example: $0 -i /data/sample_R1.fq.gz -I /data/sample_R2.fq.gz -u 8 -l 250
###############################################################################

# Parse command line arguments with getopts
while getopts "i:I:u::l:t:" opt; do
  case $opt in
    i) R1_IN="$OPTARG" ;;
    I) R2_IN="$OPTARG" ;;
    u) DO_UMI=1; UMI_LEN="${2:-8}" ; shift 2 ;;  # -u [N], default 8
    l) MIN_LEN="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
  esac
done

# Validate required arguments
if [[ -z "$R1_IN" || -z "$R2_IN" ]]; then
    echo "Usage: $0 -i input_R1.fq.gz -I input_R2.fq.gz [-u [UMI_LEN]] [-l MIN_LEN] [-t THREADS]"
    echo "  -u [N] : Extract UMI from first N nt (default N=8 if -u alone)"
    echo "Example: $0 -i sample_R1.fq.gz -I sample_R2.fq.gz -l 250"
    echo "Example: $0 -i sample_R1.fq.gz -I sample_R2.fq.gz -u -l 250"
    echo "Example: $0 -i /data/sample_R1.fq.gz -I /data/sample_R2.fq.gz -u 12 -l 250"
    exit 1
fi

# Validate input files exist
if [[ ! -f "$R1_IN" || ! -f "$R2_IN" ]]; then
    echo "Error: Input files not found: $R1_IN, $R2_IN" >&2
    exit 1
fi

# Set defaults
MIN_LEN=${MIN_LEN:-250}
THREADS=${THREADS:-8}
DO_UMI=${DO_UMI:-0}
UMI_LEN=${UMI_LEN:-8}
BASE=$(basename "${R1_IN%%_*}")  # sample_R1.fq.gz → "sample_R1" → strip _R1 later if needed
BASE=${BASE%%_R*}               # sample_R1 → "sample" (robust for folder paths)
mkdir -p good failed_quality failed_length

# Step 0: UMI extraction (OPTIONAL)
if [[ $DO_UMI -eq 1 ]]; then
    echo "Step 0: UMI extraction (first $UMI_LEN nt R1+R2 → read names)..."
    fastp -i "$R1_IN" -I "$R2_IN" \
          -o "good/${BASE}_R1_umi.fq.gz" -O "good/${BASE}_R2_umi.fq.gz" \
          --umi --umi_loc=per_read --umi_len=$UMI_LEN --umi_skip=0 \
          --thread $THREADS -h "${BASE}_umi.html" -j "${BASE}_umi.json"
    UMI_R1="good/${BASE}_R1_umi.fq.gz"
    UMI_R2="good/${BASE}_R2_umi.fq.gz"
else
    UMI_R1="$R1_IN"
    UMI_R2="$R2_IN"
fi

# Step 1: Quality filtering
echo "Step 1: Quality filtering (avg Q30)..."
fastp -i "$UMI_R1" -I "$UMI_R2" \
      -o "good/${BASE}_R1_qc.fq.gz" -O "good/${BASE}_R2_qc.fq.gz" \
      --failed_out "failed_quality/${BASE}_R1_failed.fq.gz" \
      --average_qual 30 -q 30 -u 10 \
      --thread $THREADS -h "${BASE}_qc.html" -j "${BASE}_qc.json"

# Step 2: Length filtering
echo "Step 2: Length filtering (>= $MIN_LEN bp)..."
fastp -i "good/${BASE}_R1_qc.fq.gz" -I "good/${BASE}_R2_qc.fq.gz" \
      -o "good/${BASE}_R1_final.fq.gz" -O "good/${BASE}_R2_final.fq.gz" \
      --failed_out "failed_length/${BASE}_R1_short.fq.gz" \
      -l $MIN_LEN \
      --thread $THREADS -h "${BASE}_length.html" -j "${BASE}_length.json"

# Cleanup intermediates
rm -f "good/${BASE}_R*_qc.fq.gz"
if [[ $DO_UMI -eq 1 ]]; then
    rm -f "good/${BASE}_R*_umi.fq.gz"
fi

echo "Done. Final good pairs: good/${BASE}_R*_final.fq.gz $([[ $DO_UMI -eq 1 ]] && echo "(UMI:$UMI_LEN nt extracted)")"
