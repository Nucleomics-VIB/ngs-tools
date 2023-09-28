#!/bin/env bash

# script: run_whatshap_parallel.sh
# run whatshap on a multi-VCF file and multiple BAM files
# split the unphased input multi-VCF to individual VCF files
# split the reference into single chromosome querries
# parallelise by phasing one chromosome and one sample per thread
#
# author:Stephane Plaisance (VIB-NC), 2023-08-31

# requirement:
# a conda env with:
# whatshap
# samtools, bcftools, bgzip, tabix
# a multi-VCF unphased input
# a fast reference file with fai index
# a sorted and indexed BAM file for each sample mapped to the reference
#
# commands from https://whatshap.readthedocs.io/en/latest/guide.html

# conda with whatshap v2.0
myenv=whatshap
source /etc/profile.d/conda.sh
conda activate ${myenv} || \
  ( echo "# the conda environment ${myenv} was not found on this machine" ;
    echo "# please read the top part of the script!" \
    && exit 1 )

logpfx="parallel_whatshap_$(date +%s).log"

#################
# USER VARIABLES
#################

# edit below to point to your own inputs

# multi-VCF unphased input database
mergedvcf="deepvariant_results/unphased_merged.vcf.gz"

# reference genome used for the mapping (with .fai index)
reffasta="reference/genome.fa"

# mapping folder in current path containing
# for each sample a BAM files post-processed with picard 
# - FixMateInformation
# - MarkDuplicates
# the files are named ${smpl}_FMIMD.bam
mappings="bwa__mappings"

# define chromosome and sample arrays
chrlist=$(cut -f 1 ${reffasta}.fai)

smpllist=($(bcftools query -l ${mergedvcf}))

# define output file name suffix
outsfx="phased.vcf.gz"

# create output folder
outfolder="parallel_whatshap_results"
mkdir -p ${outfolder}/tmpfiles

# export variables to have them in the custom function ran in parallel sub-shells
export "outfolder" "reffasta" "mergedvcf" "logpfx" "outsfx" "chrlist"

# run in parallel with N jobs
# Note: > 20GB RAM can be used at peak for human sized genome
pjob=10

################################################################################
# function to create chr subsets of the VCF for each sample (when not existing)
################################################################################

subset_vcf() {
local chr=$1;
local smpl=$2;

# create VCF chr subset of the full VCF
local tmp_vcf="${outfolder}/tmpfiles/tmp_${smpl}_${chr}_merged.vcf.gz";

# run only once, skip if already existing
if [ ! -f "${tmp_vcf}" ] ; then
echo "Creating subset for sample ${smpl} and chromosome ${chr}" >&2
bcftools view --threads 2 \
  -Oz \
  -s "${smpl}" \
  -r "${chr}" \
  -o "${tmp_vcf}" "${mergedvcf}" && \
    bcftools index -t "${tmp_vcf}"
fi
}

#######################################################################################
# function to run parallel jobs for each chromosome and each sample using GNU Parallel
#######################################################################################

run_whatshap() {
local chr=$1;
local smpl=$2;

# set VCF chr subset
local tmp_vcf="${outfolder}/tmpfiles/tmp_${smpl}_${chr}_merged.vcf.gz";

local bamfile=${mappings}/${smpl}_FMIMD.bam
local logfile="${outfolder}/${smpl}_${chr}_${logpfx}";

# run command with a single chr and a single sample
cmd="time whatshap phase \
  --reference ${reffasta} \
  --chromosome ${chr} \
  --sample ${smpl} \
  ${tmp_vcf} \
  ${bamfile} \
  2> >(tee -a "${logfile}" >&2)";

echo -e "# ${cmd}\n\n" > ${logfile};

eval ${cmd} | bgzip -c > ${outfolder}/${smpl}_${chr}_${outsfx} && \
    tabix -p vcf ${outfolder}/${smpl}_${chr}_${outsfx}
}

###########
# PIPELINE
###########

# export funtions to have them in parallel subshells
export -f subset_vcf
export -f run_whatshap

##############################
# split VCF into small subsets
parallel -j "${pjob}" subset_vcf ::: "${chrlist[@]}" ::: "${smpllist[@]}"

##############################
# run parallel whatshap jobs
parallel -j "${pjob}" run_whatshap ::: "${chrlist[@]}" ::: "${smpllist[@]}"

##########################################################
# run parallel whatshap jobs from a text file job_list.txt
# with a parameter pair on each line: sample,chromosome
#echo -e "sample1,chr1\nsample2,chr2" > job_list.txt
#parallel -j "${pjob}" --link --colsep ',' run_whatshap {2} {1} :::: "job_list.txt"

conda deactivate
