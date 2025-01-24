#!/bin/bash

# script: fastp_benchmark.sh
# author: SP@NC (+AI)
# 2024-01-24; v1.0
# Description: This script benchmarks fastp performance using 1 to 8 threads
# for paired-end reads. It runs fastp on the input files, times each execution,
# and saves the results to a benchmarking_result folder.

# Check if required arguments are provided
# 4920_01_A1_20487881F4509f06A1_S5_L001_R1_001.fastq.gz
# 4920_01_A1_20487881F4509f06A1_S5_L001_R2_001.fastq.gz

if [ $# -ne 1 ]; then
    echo "Usage: $0 <input>_R1_001.fastq.gz"
    exit 1
fi

# deduce R2 from R1
r1=$1
r2=${r1/R1/R2}

# number of reads to process 0=all
limit=10000000

pfx=$(basename $r1 | sed 's/_R1.*\.fastq\.gz$//')
echo $pfx

# Activate conda environment
myenv=fastp
source /etc/profile.d/conda.sh
conda activate ${myenv} || \
  ( echo "# the conda environment ${myenv} was not found on this machine" ;
    echo "# please read the top part of the script!" \
    && exit 1 )

# Create output directories
mkdir -p benchmarking_result

# Function to run fastp with specified number of threads
run_fastp() {
    nth=$1
    output_prefix="benchmarking_result/${pfx}_${nth}threads"

    start_time=$(date +%s.%N)
    fastp -i ${r1} -o ${output_prefix}_R1.fq.gz \
          -I ${r2} -O ${output_prefix}_R2.fq.gz \
          -w ${nth} \
          --reads_to_process ${limit} \
          -R ${pfx}_${lane}_${nth}threads_report \
          -j ${output_prefix}_report.json \
          -h ${output_prefix}_report.html > /dev/null 2>&1
    end_time=$(date +%s.%N)

    runtime=$(echo "${end_time} - ${start_time}" | bc)
    echo "${nth},${runtime}" >> benchmarking_result/fastp_benchmark_results.csv
}

# Initialize results file
echo "Threads,Runtime(s)" > benchmarking_result/fastp_benchmark_results.csv

# Run benchmarks
for nth in {1..8}; do
    echo "Running fastp with ${nth} thread(s)..."
    run_fastp ${nth}
done

# Generate results table
echo "Benchmark results:"
column -t -s',' benchmarking_result/fastp_benchmark_results.csv | sed 's/^/| /' | sed 's/$/ |/'
echo "Results saved in benchmarking_result/fastp_benchmark_results.csv"

