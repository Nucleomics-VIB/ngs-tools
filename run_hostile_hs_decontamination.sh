#!/bin/bash
#
# Script: run_hostile_hs_decontamination.sh
# Usage: ./run_hostile_hs_decontamination.sh -1 <read1.fastq.gz> -2 <read2.fastq.gz> -o <output_prefix> [OPTIONS]
# Description: Decontaminates paired-end FASTQ files using hostile by removing reads aligning to the human genome.
# Hostile Repository: https://github.com/czbiohub/hostile
# Author: SP@NC (+AI)
# Date: 2025-03-17
# Version: 1.0
#

# --- Configuration ---
# Required: Path to your human genome index for hostile
HOSTILE_INDEX="/data/biodata/hostile_db/human-t2t-hla"
THREADS=8           # Adjust as needed
CONDA_ENV="hostile" #The name of the conda environment
OUTPUT="./decontaminated"  # Default output directory

# --- Check Conda Environment ---
if ! conda info --envs | grep -q "$CONDA_ENV"; then
  echo "Error: Conda environment '$CONDA_ENV' not found." >&2
  echo "Please create and activate the environment before running this script." >&2
  echo "  Example: conda create -n $CONDA_ENV -c bioconda hostile" >&2
  exit 1
fi

# --- Activate Conda Environment ---
source $(conda info --base)/etc/profile.d/conda.sh # Needed for newer versions of conda
conda activate "$CONDA_ENV" || { echo "Error activating conda environment '$CONDA_ENV'." >&2; exit 1; }

# --- Function Usage ---
usage() {
  echo "Usage: $0 -1 <read1.fastq.gz> -2 <read2.fastq.gz> -o <output_prefix> [OPTIONS]"
  echo ""
  echo "Description: Decontaminates paired-end FASTQ files using hostile clean."
  echo ""
  echo "Required Arguments:"
  echo "  -1, --read1 <file>      Path to read 1 FASTQ file (gzipped)."
  echo "  -2, --read2 <file>      Path to read 2 FASTQ file (gzipped)."
  echo "  -o, --output <prefix>   Output file prefix. Output files will be placed in $OUTPUT and named:"
  echo "                          $OUTPUT/<prefix>.R1.clean.fastq.gz"
  echo "                          $OUTPUT/<prefix>.R2.clean.fastq.gz"
  echo ""
  echo "Optional Arguments:"
  echo "  -a, --aligner {bowtie2,minimap2,auto} Alignment algorithm (auto by default)"
  echo "  -i, --index <path>      Path to the hostile index.  Defaults to '$HOSTILE_INDEX'"
  echo "  -t, --threads <int>     Number of threads to use. Defaults to '$THREADS'"
  echo "  -h, --help              Show this help message and exit."
  echo ""
  echo "Example:"
  echo "  $0 -1 read1.fastq.gz -2 read2.fastq.gz -o output -t 8 -a bowtie2 -i /path/to/index"
  exit 1
}

# --- Argument Parsing ---
# Initialize variables
READ1=""
READ2=""
OUTPUT="decontaminated"
INDEX="$HOSTILE_INDEX"  # Use default
THREADS="$THREADS"
ALIGNER="auto"

# Use getopt for robust argument parsing
while getopts "1:2:o:i:t:a:h" opt; do
  case "$opt" in
    1) READ1="$OPTARG" ;;
    2) READ2="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
    i) INDEX="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    a) ALIGNER="$OPTARG" ;;
    h) usage ;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
    :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# Check for required arguments
if [ -z "$READ1" ] || [ -z "$READ2" ] || [ -z "$OUTPUT" ]; then
  echo "Error: Missing required arguments." >&2
  usage
fi

# Check if input files exist
if [ ! -f "$READ1" ]; then
  echo "Error: Read 1 file not found: $READ1" >&2
  exit 1
fi

if [ ! -f "$READ2" ]; then
  echo "Error: Read 2 file not found: $READ2" >&2
  exit 1
fi

# Create output directory
mkdir -p "$OUTPUT"

# --- Hostile Execution ---
echo "Running hostile clean with the following parameters:"
echo "  Read 1: $READ1"
echo "  Read 2: $READ2"
echo "  Aligner: $ALIGNER"
echo "  Index: $INDEX"
echo "  Threads: $THREADS"
echo "  Output Directory: $OUTPUT"

# Construct the hostile clean command
HOSTILE_CMD="hostile clean --fastq1 \"$READ1\" --fastq2 \"$READ2\" --index \"$INDEX\" -t \"$THREADS\" -o \"$OUTPUT\" --aligner \"$ALIGNER\""

# Execute the hostile clean command
echo "# command: $HOSTILE_CMD"
eval "$HOSTILE_CMD"

# Check if hostile clean ran successfully
if [ $? -ne 0 ]; then
  echo "Error: hostile clean failed.  Check the hostile output for details." >&2
  exit 1
fi

echo "Hostile decontamination complete."
exit 0
