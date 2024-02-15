#!/bin/bash

# script: fastqsplit.sh
# SP@NC 2024-02-15, v1.0
# depends on fastqsplitter (https://github.com/LUMC/fastqsplitter)

myenv=fastqsplitter
source /etc/profile.d/conda.sh
conda activate ${myenv} || \
  ( echo "# the conda environment ${myenv} was not found on this machine" ;
    echo "# please read the top part of the script!" \
    && exit 1 )

# Check if fastqsplitter command is available
if ! command -v fastqsplitter &> /dev/null; then
    echo "Error: fastqsplitter command not found. Please install fastqsplitter (https://github.com/LUMC/fastqsplitter)" >&2
    exit 1
fi

usage="Usage: $(basename "$0") -i <input_fastq> [-c <chunks (8)>] [-t <num_threads (2)>]"

while getopts ":i:c:t:" opt; do
    case $opt in
        i) infastq=$OPTARG ;;
        c) opt_chunks=$OPTARG ;;
        t) opt_nthr=$OPTARG ;;
        *) echo ${usage} >&2
           exit 1 ;;
    esac
done

# Default values
chunks=${opt_chunks:-8}
nthr=${opt_nthr:-2}

# Check if input fastq is provided
if [ -z "${infastq}" ]; then
    echo "Error: Input FASTQ file not provided." >&2
    echo ${usage} >&2
    exit 1
fi

# Extract prefix from input fastq
pfx=$(basename "${infastq%.f*}")
outputs=$(for chunk in $(seq -f "${pfx}_chunk_%g.fq.gz" 1 $chunks); do echo -n "-o $chunk "; done)

# Generate command for fastqsplitter
cmd="fastqsplitter -i ${infastq} ${outputs} -t ${nthr}"

# Execute command
echo "# ${cmd}"
eval ${cmd}
