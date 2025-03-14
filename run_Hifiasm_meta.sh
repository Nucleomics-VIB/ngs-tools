#!/bin/bash

# Script: run_Hifiasm_meta.sh
# Description: Process PacBio HiFi data through conversion (if needed) and assembly
#              Two modes are available:
#              1. BAM mode (default): Converts BAM to FASTQ, then assembles
#              2. FASTQ mode: Directly assembles from existing FASTQ files
# Author: SP@NC (+AI)
# Version: 1.0; 2025-03-14

version="1.0; 2025-03-14"

usage() {
    echo "Usage: $0 [-q] [-b INPUT_BAM_DIR] [-f FASTQ_DIR] [-o OUTPUT_ASM_DIR] [-n THREADS] [-j MAX_JOBS] [-t ASM_THREADS]"
    echo "Options:"
    echo "  -q          start directly from FASTQ files (default: process BAM files)"
    echo "  -b DIR      input BAM directory (default: bam_data)"
    echo "  -f DIR      FASTQ input/output directory (default: fastq_data)"
    echo "  -o DIR      assembly output directory (default: asm_results)"
    echo "  -n INT      threads for bam2fastq (default: 4)"
    echo "  -j INT      max concurrent jobs (default: 4)"
    echo "  -t INT      threads per hifiasm_meta (default: 20)"
    echo " version ${version}"
    exit 1
}

# Initialize default values
start_from_fastq=false
inbam="bam_data"
infastq="fastq_data"
outasm="asm_results"
n=4
j=4
t=20

# Check if required tools are installed
check_tools() {
    for tool in bam2fastq hifiasm_meta parallel; do
        if ! command -v $tool &> /dev/null; then
            echo "ERROR: $tool is not installed. Please install it and try again."
            exit 1
        fi
    done
}

check_tools

# Parse command-line options
while getopts ":qb:f:o:n:j:t:h" opt; do
    case $opt in
        q) start_from_fastq=true ;;
        b) inbam="$OPTARG" ;;
        f) infastq="$OPTARG" ;;
        o) outasm="$OPTARG" ;;
        n) n="$OPTARG" ;;
        j) j="$OPTARG" ;;
        t) t="$OPTARG" ;;
        h) usage ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# Create log directory
log_dir="hifiasm_meta_logs_$(date +%Y%m%d_%H%M%S)"
mkdir -p "${log_dir}"

# Check directories based on the starting point
if [[ "$start_from_fastq" == false ]]; then
    # Check BAM directory
    if [[ ! -d "$inbam" ]] || [[ -z "$(find "$inbam" -type f -name "*.bam" -print -quit)" ]]; then
        echo "ERROR: Input BAM directory '$inbam' does not exist or is empty."
        exit 1
    fi
else
    # Check FASTQ directory for .fastq.gz files
    if [[ ! -d "$infastq" ]] || [[ -z "$(find "$infastq" -type f -name "*.fastq.gz" -print -quit)" ]]; then
        echo "ERROR: Input FASTQ directory '$infastq' does not exist or contains no .fastq.gz files."
        exit 1
    fi
fi

echo "Starting pipeline. Logs will be saved in ${log_dir}."

process_bam() {
    local bam="$1"
    local pfx=$(basename "${bam}" .bam)
    local fastq_flag="${infastq}/${pfx}.fastq.gz.ok"
    local log_file="${log_dir}/${pfx}_processing.log"

    mkdir -p "${infastq}"

    echo "Processing ${pfx} - Log: ${log_file}"
    echo "Starting processing of ${pfx} at $(date)" > "${log_file}"

    # BAM to FASTQ conversion
    if [[ -f "${fastq_flag}" ]]; then
        echo "Skipping FASTQ conversion for ${pfx} (already completed)" | tee -a "${log_file}"
    else
        echo "Converting ${pfx} to FASTQ..." | tee -a "${log_file}"
        local conv_start=$(date +%s)
        
        local cmd="bam2fastq -j ${n} ${bam} -o ${infastq}/${pfx}"
        echo "Executing command: ${cmd}" | tee -a "${log_file}"
        eval ${cmd} 2>> "${log_file}"
        local exit_status=$?
        local conv_end=$(date +%s)
        local conv_duration=$((conv_end - conv_start))

        if [[ ${exit_status} -eq 0 ]] && [[ -f "${infastq}/${pfx}.fastq.gz" ]]; then
            touch "${fastq_flag}"
            echo "FASTQ conversion successful - Duration: ${conv_duration}s" | tee -a "${log_file}"
        else
            echo "WARNING: Conversion issues detected for ${pfx}" | tee -a "${log_file}"
            echo "FASTQ conversion failed - Duration: ${conv_duration}s" | tee -a "${log_file}"
            return 1
        fi
    fi

    # Call the assembly function for this prefix
    assemble_fastq "${pfx}" "${log_file}"
}

process_fastq() {
    local fastq="$1"
    local pfx=$(basename "$fastq" .fastq.gz)
    local log_file="${log_dir}/${pfx}_processing.log"

    echo "Processing ${pfx} - Log: ${log_file}"
    echo "Starting processing of ${pfx} at $(date)" > "${log_file}"

    assemble_fastq "$pfx" "$log_file"
}

assemble_fastq() {
    local pfx="$1"
    local log_file="$2"
    local asm_dir="${outasm}/${pfx}_asm"
    local asm_flag="${asm_dir}/assembly.ok"

    if [[ -f "${asm_flag}" ]]; then
        echo "Skipping assembly for ${pfx} (already completed)" | tee -a "${log_file}"
    else
        if [[ -f "${infastq}/${pfx}.fastq.gz" ]]; then
            echo "Assembling ${pfx}..." | tee -a "${log_file}"
            mkdir -p "${asm_dir}"
            local asm_start=$(date +%s)
            
            local cmd="hifiasm_meta -t ${t} -o ${asm_dir}/${pfx} ${infastq}/${pfx}.fastq.gz"
            echo "Executing command: ${cmd}" | tee -a "${log_file}"
            eval ${cmd} 2>> "${log_file}"
            local exit_status=$?
            local asm_end=$(date +%s)
            local asm_duration=$((asm_end - asm_start))

            if [[ ${exit_status} -eq 0 ]]; then
                touch "${asm_flag}"
                echo "Assembly successful - Duration: ${asm_duration}s" | tee -a "${log_file}"
            else
                echo "WARNING: Assembly issues detected for ${pfx}" | tee -a "${log_file}"
                echo "Assembly failed - Duration: ${asm_duration}s" | tee -a "${log_file}"
                return 1
            fi
        else
            echo "Skipping assembly for ${pfx} (missing FASTQ)" | tee -a "${log_file}"
            return 1
        fi
    fi

    echo "Completed processing of ${pfx} at $(date)" | tee -a "${log_file}"
}

# ==========================================
# Main Execution Logic Based on Starting Point
# ==========================================

if [[ "$start_from_fastq" == true ]]; then
    export -f process_fastq assemble_fastq
    export infastq outasm t log_dir

    echo "Starting assembly from FASTQ files with ${j} concurrent jobs..."
    cmd="find ${infastq} -name '*.fastq.gz' -print0 | parallel -0 -j ${j} process_fastq {}"
    echo "Executing command: ${cmd}" | tee -a "${log_dir}/main_execution.log"
    eval ${cmd}
else
    export -f process_bam assemble_fastq
    export infastq outasm n t log_dir

    echo "Starting BAM processing with ${j} concurrent jobs..."
    cmd="find ${inbam} -name '*.bam' -print0 | parallel -0 -j ${j} process_bam {}"
    echo "Executing command: ${cmd}" | tee -a "${log_dir}/main_execution.log"
    eval ${cmd}
fi

echo "Pipeline completed. Logs saved in: ${log_dir}"

exit 0
