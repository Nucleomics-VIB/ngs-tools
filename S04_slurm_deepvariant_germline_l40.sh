#!/bin/bash

#SBATCH --job-name="S04_DV_germline"
#SBATCH --output="%x_%j.out"
#SBATCH --error="%x_%j.err"
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=4
#SBATCH --time=0-03:00:00
#SBATCH --mem=256G
#SBATCH --account=s04
#SBATCH --gres=gpu:4
#SBATCH --partition=gpu_l40s_64C_128T_1TB

# Aim: run the NVidia clara parabricks deepvariant_germline (v4.4.0) on 1 sample read pair
# script: slurm_deepvariant_germline.sh
#   run with 'sbatch scripts/slurm_deepvariant_germline.sh reads/SRR29676022_1.fq.gz 4'
# also logs usage during run and ends logging before quitting
# author: SP@NC; 2025-01-21 v1.1

# Check if an argument is provided
if [ $# -eq 0 ]; then
    echo "Error: Please provide the path to the reads_1.fq.gz file. \
          Usage: sbatch scripts/slurm_deepvariant_germline.sh reads/<read1.fq.gz>"
    exit 1
fi

# Function to stop logging and wait for the process to finish
stop_logging() {
    if [ -n "$LOGGING_PID" ]; then
        kill $LOGGING_PID
        wait $LOGGING_PID 2>/dev/null
    fi
}

# Trap to ensure logging is stopped even if the script exits unexpectedly
trap stop_logging EXIT

# set the $1 argument to the path to read_1.fq.gz provided by the user
fq=${1}

# use 4 GPU by default unless stated with $2
numgpu=${2:-4}

# set WORKDIR
WORKDIR="/data/projects/s04/wgs_variant_analysis"

# set pfx
pfx=$(basename "${fq%_R1.fq.gz}")_${numgpu}

# output folder for mappings
outbam="mappings_${numgpu}gpu_${pfx}"
mkdir -p "${outbam}"

# output folder for variants
outvcf="variants_XYhap_${numgpu}gpu_${pfx}"
mkdir -p "${outvcf}"

# output folder for logs
outlogs="logs_XYhap_${numgpu}gpu_${pfx}"
mkdir -p "${outlogs}"

# start logging usage from custom scripts
interval=5
./scripts/logUsageGPU.sh -i ${interval} -o ${outlogs} &
LOGGING_PID=$!

#######################
# additional arguments
#######################

# locate the singularity SIF file for the latest "NVidia clara parabricks" workflow
singimg="${WORKDIR}/bin/clara-parabricks_4.4.0-1.sif"

# run settings (platform-dependent, this should migrate into a config file)
refidx="mRatBN7.2.fa"
knownsites="rattus_norvegicus.vcf.gz"
optdist=100
platform="Aviti"

# optional parameter for male rat genomes
# list of known haploid chromosomes in this assembly (X, Y, Y_unloc4)
hapchr="X,Y,MU150192.1"

# deduce fq1
fq1=$(basename "${fq}")

# deduce fq2
fq2=${fq1/_R1.fq.gz/_R2.fq.gz}

# build and run singularity command
echo "Starting slurm_deepvariant_germline.sh job execution..."

singularity run --nv \
    --bind "${WORKDIR}:/workdir" \
    --bind "${TMPDIR}:/tmp" \
    --pwd "/workdir" \
    ${singimg} \
    pbrun deepvariant_germline \
    --gpusort \
    --gpuwrite \
    --num-gpus ${numgpu} \
    --bwa-options '-M -K 10000000' \
    --optical-duplicate-pixel-distance "${optdist}" \
    --read-group-sm "${pfx}" \
    --read-group-lb "lib_${pfx}" \
    --read-group-pl "${platform}" \
    --read-group-id-prefix "${pfx}" \
    --ref "/workdir/bwaidx/${refidx}" \
    --in-fq "/workdir/reads/${fq1}" "/workdir/reads/${fq2}" \
    --knownSites "/workdir/reference/${knownsites}" \
    --haploid-contigs "${hapchr}" \
    --out-bam "/workdir/${outbam}/${pfx}_mrkdup_recal.bam" \
    --out-recal-file "/workdir/${outbam}/${pfx}_recal.txt" \
    --out-duplicate-metrics "/workdir/${outbam}/${pfx}_duplicate_metrics" \
    --out-variants "/workdir/${outvcf}/${pfx}_XYhap.g.vcf.gz" \
    --gvcf \
    --tmp-dir "/tmp"

# cleanup

# rename output and error files to include ${pfx} and move them
mv ${WORKDIR}/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out ${outlogs}/${pfx}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out
mv ${WORKDIR}/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err ${outlogs}/${pfx}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err

echo "Cleanup complete. Exiting."

# terminate slurm queued session now to run next in queue
exit 0
