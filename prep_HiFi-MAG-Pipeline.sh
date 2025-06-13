#!/bin/bash

# Title       : prep_HiFi-MAG-Pipeline.sh
# Description : Prepare hifiasm_meta assembly and HiFi reads for the HiFi-MAG-Pipeline.
# link:       : https://github.com/PacificBiosciences/pb-metagenomics-tools/blob/master/docs/Tutorial-HiFi-MAG-Pipeline.md
# Author      : SP@NC
# Version     : 1.1; 2023-06-13

# Check for required commands
for cmd in awk sed zcat grep; do
    if ! command -v $cmd &> /dev/null; then
        echo "Error: Required command '$cmd' not found."
        exit 1
    fi
done

# Usage function
usage() {
    echo "Usage: $0 -c <path_to_contig_gfa> -r <path_to_reads_fastq.gz> -s <sample_name>"
    echo "Options:"
    echo "  -c   Path to the contig GFA file from hifiasm_meta."
    echo "  -r   Path to the HiFi reads FASTQ.GZ file."
    echo "  -s   Sample name to use for renaming files (alphanumeric and underscore only)."
    exit 1
}

# Parse command-line options
while getopts "c:r:s:" opt; do
    case $opt in
        c) CONTIG_GFA="$OPTARG" ;;
        r) READS_FASTQ_GZ="$OPTARG" ;;
        s) SAMPLE_NAME="$OPTARG" ;;
        *) usage ;;
    esac
done

# Check if all required arguments are provided
if [ -z "$CONTIG_GFA" ] || [ -z "$READS_FASTQ_GZ" ] || [ -z "$SAMPLE_NAME" ]; then
    echo "Error: Missing required arguments."
    usage
fi

# Check if input files exist
if [ ! -f "$CONTIG_GFA" ]; then
    echo "Error: Contig GFA file '$CONTIG_GFA' not found."
    exit 1
fi

if [ ! -f "$READS_FASTQ_GZ" ]; then
    echo "Error: Reads FASTQ.GZ file '$READS_FASTQ_GZ' not found."
    exit 1
fi

# Validate input file formats (basic checks)
if ! head -n 1 "$CONTIG_GFA" | grep -q "^S"; then
    echo "Warning: Contig GFA file may not be in the expected format. First line doesn't start with 'S'."
fi

if ! gzip -t "$READS_FASTQ_GZ" 2>/dev/null; then
    echo "Error: Reads file is not a valid gzip file."
    exit 1
fi

# Create the inputs directory if it doesn't exist
mkdir -p inputs

# Convert the contig GFA file to FASTA and rename it to SAMPLE_contigs.fasta
echo "Converting $CONTIG_GFA to FASTA and saving as inputs/${SAMPLE_NAME}_contigs.fasta..."
awk '/^S/{print ">"$2"\n"$3}' "$CONTIG_GFA" > "inputs/${SAMPLE_NAME}_contigs.fasta"
if [ $? -ne 0 ]; then
    echo "Error converting contig GFA to FASTA."
    exit 1
fi

# Convert the reads FASTQ.GZ file to FASTA and rename it to SAMPLE.fasta
echo "Converting $READS_FASTQ_GZ to FASTA and saving as inputs/${SAMPLE_NAME}.fasta..."
zcat "$READS_FASTQ_GZ" | sed -n '1~4s/^@/>/p;2~4p' > "inputs/${SAMPLE_NAME}.fasta"
if [ $? -ne 0 ]; then
    echo "Error converting reads FASTQ.GZ to FASTA."
    exit 1
fi

# Update Sample-Config.yaml
CONFIG_FILE="configs/Sample-Config.yaml"
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Creating new Sample-Config.yaml file..."
    mkdir -p configs
    echo "samplenames:" > "$CONFIG_FILE"
fi

# Check if the sample is already in the config file
if ! grep -q "  - $SAMPLE_NAME" "$CONFIG_FILE"; then
    echo "Adding $SAMPLE_NAME to Sample-Config.yaml..."
    echo "  - $SAMPLE_NAME" >> "$CONFIG_FILE"
else
    echo "Sample $SAMPLE_NAME already exists in Sample-Config.yaml."
fi

echo "Files successfully prepared:"
echo "  Contigs: inputs/${SAMPLE_NAME}_contigs.fasta"
echo "  Reads:   inputs/${SAMPLE_NAME}.fasta"
echo "Sample-Config.yaml updated in configs/ directory."
