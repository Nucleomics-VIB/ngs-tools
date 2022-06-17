#!/bin/bash

# call variants from long reads and a reference
# using PEPPER Margin Deepvariant
#
# Requirements:
# you will be asked your sudo password (no sudo no run!)
# run on a unix computer installed with
#   docker image: kishwars/pepper_deepvariant:r0.8
#   # get with: sudo docker pull kishwars/pepper_deepvariant:r0.8
#   samtools, minimap2,
#   reference genome fasta present and indexed
#   a minimum of 4 threads
#
# this script will write results in the current folder
#
# Stephane Plaisance (VIB-NC) 2022/06/01; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# version number
version="1.0, 2022_06_01"

usage='# Usage: run_PMDeepvariant.sh -f <fastq reads> -r <fasta reference>
# script version '${version}'
# [optional: -n <sample name|sample>]
# [optional: -t <threads|4>]
# [optional: -X <extra parameter string between quotes for the docker command>]'

while getopts "f:r:n:t:X:h" opt; do
  case $opt in
    f) reads=${OPTARG} ;;
    r) ref=${OPTARG} ;;
    n) opt_name=${OPTARG} ;;
    t) opt_thr=${OPTARG} ;;
    X) opt_extra=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2;
       exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; 
       exit 1 ;;
  esac
done

# set defaults
nthr=${opt_thr:-4}
stthr=4
extra=${opt_extra:-""}
pfx=${opt_name:-"sample"}

# check if all dependencies are present
declare -a arr=("minimap2" "samtools")
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || \
  ( echo "# required ${prog} not found in PATH"
exit 1 )
done

# check for the docker image
res=$(sudo docker images | grep kishwars/pepper_deepvariant | cut -d " " -f 1)
if [ ! ${res} == "kishwars/pepper_deepvariant" ]; then
  (echo "docker image not found"; exit 1)
fi

# test if minimal arguments were provided
if [ -z "${reads}" ]; then
   echo "# no long reads provided!"
   echo "${usage}"
   exit 1
fi

if [ ! -f "${reads}" ]; then
    echo "${reads} file not found!";
    exit 1
fi

if [ -z "${ref}" ]
then
	echo "# no fasta reference provided"
	echo "${usage}"
	exit 1
fi

if [ ! -f "${ref}" ]; then
    echo "${ref} file not found!";
    exit 1
fi

# other parameters or defaults
workdir=$PWD

# create folders for results
mkdir -p input output/intermediate_files output/logs

#############################################
# map reads to the reference using minimap2 #
#############################################

reffile=$(basename ${ref})
cp ${ref} input/${reffile} && \
  samtools faidx input/${reffile}

# build the command
cmd="minimap2 \
  -a \
  -t ${nthr} \
  -x map-hifi \
  -R '@RG\tID:'${pfx}'\tSM:'${pfx} \
  --MD \
  input/${reffile} \
  ${reads} \
  > output/${pfx}_mappings.sam && \
    samtools sort -@${stthr} \
      output/${pfx}_mappings.sam \
      -O BAM \
      -o output/${pfx}_sorted.bam && \
        samtools index output/${pfx}_sorted.bam && \
          rm output/${pfx}_mappings.sam"

# show and execute
echo "# mapping reads to reference"
echo "# ${cmd}"
time eval ${cmd}

#################################
# run pepper margin deepvariant #
#################################

input_dir="${workdir}/input"
output_dir="${workdir}/output"
out_pfx="${pfx}_PEPPER_Margin_DeepVariant"

param_haplotag="allParams.haplotag.pb-hifi.hapDup.json"
param_phase="allParams.phase_vcf.pb-hifi.json"

# Run PEPPER-Margin-DeepVariant
echo "# running pepper_margin_deepvariant docker"
time sudo docker run \
  -u "$(id -u):$(id -g)" \
  -v "${input_dir}":"/input" \
  -v "${output_dir}":"/output" \
  kishwars/pepper_deepvariant:r0.8 \
    run_pepper_margin_deepvariant call_variant \
    -b "/output/${pfx}_sorted.bam" \
    -f "/input/${reffile}" \
    -o "/output" \
    -p "${out_pfx}"_q \
    -t "${nthr}" \
    -s ${pfx} \
    --keep_intermediate_bam_files \
    --sample_name "${pfx}" \
    --gvcf \
    --phased_output \
    --margin_haplotag_model "/opt/margin_dir/params/phase/${param_haplotag}" \
    --margin_phase_model "/opt/margin_dir/params/phase/${param_phase}" \
    --hifi \
    ${extra}

exit 0

######################################################
################### Command info #####################
######################################################

usage: run_pepper_margin_deepvariant call_variant -b BAM -f FASTA -o
                                                  OUTPUT_DIR -t THREADS
                                                  [-r REGION] [-g] [--gvcf]
                                                  [-s SAMPLE_NAME]
                                                  [--only_pepper]
                                                  [--only_haplotag]
                                                  [--phased_output]
                                                  [--skip_final_phased_bam]
                                                  [-p OUTPUT_PREFIX] [-k]
                                                  [--dry]
                                                  [--pepper_model PEPPER_MODEL]
                                                  [--pepper_quantized]
                                                  [--no_pepper_quantized]
                                                  [--pepper_downsample_rate PEPPER_DOWNSAMPLE_RATE]
                                                  [--pepper_region_size PEPPER_REGION_SIZE]
                                                  [--pepper_include_supplementary]
                                                  [--pepper_min_mapq PEPPER_MIN_MAPQ]
                                                  [--pepper_min_snp_baseq PEPPER_MIN_SNP_BASEQ]
                                                  [--pepper_min_indel_baseq PEPPER_MIN_INDEL_BASEQ]
                                                  [--pepper_snp_frequency PEPPER_SNP_FREQUENCY]
                                                  [--pepper_insert_frequency PEPPER_INSERT_FREQUENCY]
                                                  [--pepper_delete_frequency PEPPER_DELETE_FREQUENCY]
                                                  [--pepper_min_coverage_threshold PEPPER_MIN_COVERAGE_THRESHOLD]
                                                  [--pepper_candidate_support_threshold PEPPER_CANDIDATE_SUPPORT_THRESHOLD]
                                                  [--pepper_snp_candidate_frequency_threshold PEPPER_SNP_CANDIDATE_FREQUENCY_THRESHOLD]
                                                  [--pepper_indel_candidate_frequency_threshold PEPPER_INDEL_CANDIDATE_FREQUENCY_THRESHOLD]
                                                  [--pepper_skip_indels]
                                                  [--pepper_allowed_multiallelics PEPPER_ALLOWED_MULTIALLELICS]
                                                  [--pepper_snp_p_value PEPPER_SNP_P_VALUE]
                                                  [--pepper_insert_p_value PEPPER_INSERT_P_VALUE]
                                                  [--pepper_delete_p_value PEPPER_DELETE_P_VALUE]
                                                  [--pepper_snp_p_value_in_lc PEPPER_SNP_P_VALUE_IN_LC]
                                                  [--pepper_insert_p_value_in_lc PEPPER_INSERT_P_VALUE_IN_LC]
                                                  [--pepper_delete_p_value_in_lc PEPPER_DELETE_P_VALUE_IN_LC]
                                                  [--pepper_snp_q_cutoff PEPPER_SNP_Q_CUTOFF]
                                                  [--pepper_indel_q_cutoff PEPPER_INDEL_Q_CUTOFF]
                                                  [--pepper_snp_q_cutoff_in_lc PEPPER_SNP_Q_CUTOFF_IN_LC]
                                                  [--pepper_indel_q_cutoff_in_lc PEPPER_INDEL_Q_CUTOFF_IN_LC]
                                                  [--pepper_report_snp_above_freq PEPPER_REPORT_SNP_ABOVE_FREQ]
                                                  [--pepper_report_indel_above_freq PEPPER_REPORT_INDEL_ABOVE_FREQ]
                                                  [--margin_haplotag_model MARGIN_HAPLOTAG_MODEL]
                                                  [--margin_phase_model MARGIN_PHASE_MODEL]
                                                  [--dv_model DV_MODEL]
                                                  [--dv_model_snp DV_MODEL_SNP]
                                                  [--dv_model_indel DV_MODEL_INDEL]
                                                  [--dv_alt_aligned_pileup DV_ALT_ALIGNED_PILEUP]
                                                  [--dv_alt_aligned_pileup_snp DV_ALT_ALIGNED_PILEUP_SNP]
                                                  [--dv_alt_aligned_pileup_indel DV_ALT_ALIGNED_PILEUP_INDEL]
                                                  [--dv_pileup_image_width DV_PILEUP_IMAGE_WIDTH]
                                                  [--dv_realign_reads DV_REALIGN_READS]
                                                  [--dv_partition_size DV_PARTITION_SIZE]
                                                  [--dv_min_mapping_quality DV_MIN_MAPPING_QUALITY]
                                                  [--dv_min_base_quality DV_MIN_BASE_QUALITY]
                                                  [--dv_vsc_min_fraction_indels DV_VSC_MIN_FRACTION_INDELS]
                                                  [--dv_vsc_min_fraction_snps DV_VSC_MIN_FRACTION_SNPS]
                                                  [--dv_sort_by_haplotypes DV_SORT_BY_HAPLOTYPES]
                                                  [--dv_parse_sam_aux_fields DV_PARSE_SAM_AUX_FIELDS]
                                                  [--dv_add_hp_channel DV_ADD_HP_CHANNEL]
                                                  [--dv_use_hp_information DV_USE_HP_INFORMATION]
                                                  [--dv_use_multiallelic_mode DV_USE_MULTIALLELIC_MODE]
                                                  (--ont_r9_guppy5_sup | --ont_r10_q20 | --hifi)
                                                  [-h]

Run PEPPER-Margin-DeepVariant for variant calling.
Example run: run_pepper_margin_deepvariant -b <BAM> -f <FASTA> -o <OUTPUT_DIR> -t <THREADS> <--ont_r9_guppy5_sup/--ont_r10_q20/--hifi>
Output file description:
1) OUTPUT_DIR/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz: Final output file. Filename can be changed with -p option.
2) OUTPUT_DIR/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.visual_report.html: Basic variant statistics of the output file.
Intermediate output files:
1) OUTPUT_DIR/intermediate_files/PEPPER_VARIANT_FULL.vcf.gz: Full set of variants from PEPPER.
2) OUTPUT_DIR/intermediate_files/PEPPER_VARIANT_OUTPUT_PEPPER.vcf.gz: Variants from PEPPER that will NOT be re-genotyped with DeepVariant.
3) OUTPUT_DIR/intermediate_files/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING.vcf.gz: All variants from PEPPER that will be re-genotyped with DeepVariant.
4) OUTPUT_DIR/intermediate_files/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING_SNPs.vcf.gz: SNP variants from PEPPER that will be re-genotyped with DeepVariant.
5) OUTPUT_DIR/intermediate_files/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING_INDEL.vcf.gz: INDEL variants from PEPPER that will be re-genotyped with DeepVariant.
6) OUTPUT_DIR/intermediate_files/DEEPVARIANT_OUTPUT_SNP.vcf.gz: SNP variants reported by PEPPER_VARIANT_OUTPUT_VARIANT_CALLING_SNPs re-genotyped with DeepVariant
7) OUTPUT_DIR/intermediate_files/DEEPVARIANT_OUTPUT_INDEL.vcf.gz: INDEL variants reported by PEPPER_VARIANT_OUTPUT_VARIANT_CALLING_SNPs re-genotyped with DeepVariant
8) OUTPUT_DIR/intermediate_files/DEEPVARIANT_OUTPUT.vcf.gz: Variants reported by PEPPER_VARIANT_OUTPUT_VARIANT_CALLING re-genotyped with DeepVariant
Note 1: DEEPVARIANT_OUTPUT_SNP.vcf.gz and DEEPVARIANT_OUTPUT_INDEL.vcf.gz is subset of PEPPER_VARIANT_OUTPUT_VARIANT_CALLING.vcf.gz.
Note 2: DEEPVARIANT_OUTPUT.vcf.gz only exists if we use PEPPER_VARIANT_OUTPUT_VARIANT_CALLING.vcf.gz as candidates and do not do SNP and INDEL calling separately.
Note 3: The final PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz file is generated by merging PEPPER_VARIANT_OUTPUT_PEPPER.vcf.gz and DEEPVARIANT_OUTPUT* files.

Required Arguments:
  -b BAM, --bam BAM     Alignment containing mapping between reads and a reference.
  -f FASTA, --fasta FASTA
                        A reference file in FASTA format.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Path to output directory.
  -t THREADS, --threads THREADS
                        Number of threads to use.
  --ont_r9_guppy5_sup   Set to call variants on R9.4.1 Guppy 5+/6+ sup/hac Oxford Nanopore reads.
  --ont_r10_q20         Set to call variants on R10.4 Q20 Oxford Nanopore reads.
  --hifi                Set to call variants on PacBio HiFi reads.

Optional arguments:
  -r REGION, --region REGION
                        Region in [contig_name:start-end] format. Default is None.
  -g, --gpu             If set then will use GPUs for inference. CUDA required.
  --gvcf                If set then a gVCF output will be generated.
  -s SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Name of the sample to put in the output VCF.
  --only_pepper         If set then pipeline will stop after PEPPER.
  --only_haplotag       If set then pipeline will stop after Margin haplotag stage.
  --phased_output       If set then Margin phase will generate a phased VCF and BAM at the end when --phased_output is set.
  --skip_final_phased_bam
                        If true with phased output then the final output will not have a bam file when --phased_output is set.
  -p OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Prefix for output filename. Do not include extension in prefix.
  -k, --keep_intermediate_bam_files
                        If true then intermediate haplotagged bam files will be kept in the intermediate files directory. Default: False
  --dry                 If true then only the commands will be printed. [Debugging]
  -h, --help            Show this text and exit.

Arguments for PEPPER:
  --pepper_model PEPPER_MODEL
                        Use a custom PEPPER model.
  --pepper_quantized    PEPPER: Use quantization for inference while on CPU inference mode. Speeds up inference. Default is True.
  --no_pepper_quantized
                        Do not use quantization for inference while on CPU inference mode. Speeds up inference.
  --pepper_downsample_rate PEPPER_DOWNSAMPLE_RATE
                        Downsample rate of reads while generating images. Default is 1.0
  --pepper_region_size PEPPER_REGION_SIZE
                        Region size in bp used to chunk the genome. Default is 100000.
  --pepper_include_supplementary
                        If true then supplementary reads will be used. Default is False.
  --pepper_min_mapq PEPPER_MIN_MAPQ
                        Minimum mapping quality for read to be considered valid. Default is 5
  --pepper_min_snp_baseq PEPPER_MIN_SNP_BASEQ
                        Minimum base quality for base to be considered valid for SNP. Default is 1
  --pepper_min_indel_baseq PEPPER_MIN_INDEL_BASEQ
                        Minimum base quality for base to be considered valid for INDELs.
  --pepper_snp_frequency PEPPER_SNP_FREQUENCY
                        Minimum SNP frequency for a site to be considered to have a variant.
  --pepper_insert_frequency PEPPER_INSERT_FREQUENCY
                        Minimum insert frequency for a site to be considered to have a variant.
  --pepper_delete_frequency PEPPER_DELETE_FREQUENCY
                        Minimum delete frequency for a site to be considered to have a variant.
  --pepper_min_coverage_threshold PEPPER_MIN_COVERAGE_THRESHOLD
                        Minimum delete frequency for a site to be considered to have a variant.
  --pepper_candidate_support_threshold PEPPER_CANDIDATE_SUPPORT_THRESHOLD
                        Minimum number of reads supporting a variant to be considered as a candidate.
  --pepper_snp_candidate_frequency_threshold PEPPER_SNP_CANDIDATE_FREQUENCY_THRESHOLD
                        Minimum frequency for a SNP candidate to be considered to be a variant.
  --pepper_indel_candidate_frequency_threshold PEPPER_INDEL_CANDIDATE_FREQUENCY_THRESHOLD
                        Minimum frequency for an INDEL candidate to be considered to be a variant.
  --pepper_skip_indels  If set then INDEL calling will be skipped.
  --pepper_allowed_multiallelics PEPPER_ALLOWED_MULTIALLELICS
                        Max allowed multiallelic variants per site.
  --pepper_snp_p_value PEPPER_SNP_P_VALUE
                        Predicted value used for a SNP to be considered a candidate.
  --pepper_insert_p_value PEPPER_INSERT_P_VALUE
                        Predicted value used for a insert to be considered a candidate.
  --pepper_delete_p_value PEPPER_DELETE_P_VALUE
                        Predicted value used for a delete to be considered a candidate.
  --pepper_snp_p_value_in_lc PEPPER_SNP_P_VALUE_IN_LC
                        Predicted value used for a SNP to be considered a candidate in low complexity regions.
  --pepper_insert_p_value_in_lc PEPPER_INSERT_P_VALUE_IN_LC
                        Predicted value used for an insert to be considered a candidate in low complexity regions.
  --pepper_delete_p_value_in_lc PEPPER_DELETE_P_VALUE_IN_LC
                        Predicted value used for a delete to be considered a candidate in low complexity regions.
  --pepper_snp_q_cutoff PEPPER_SNP_Q_CUTOFF
                        GQ cutoff for a SNP variant to be re-genotyped with DeepVariant. Variants with GQ below this will be re-genotyped.
  --pepper_indel_q_cutoff PEPPER_INDEL_Q_CUTOFF
                        GQ cutoff for an INDEL variant to be re-genotyped with DeepVariant. Variants with GQ below this will be re-genotyped.
  --pepper_snp_q_cutoff_in_lc PEPPER_SNP_Q_CUTOFF_IN_LC
                        GQ cutoff for a SNP variant in low complexity region to be re-genotyped with DeepVariant. Variants with GQ below this will be re-genotyped.
  --pepper_indel_q_cutoff_in_lc PEPPER_INDEL_Q_CUTOFF_IN_LC
                        GQ cutoff for an INDEL variant in low complexity region to be re-genotyped with DeepVariant. Variants with GQ below this will be re-genotyped.
  --pepper_report_snp_above_freq PEPPER_REPORT_SNP_ABOVE_FREQ
                        Report all SNPs above frequency for re-genotyping even if the predicted value is low. Set 0 to disable this.
  --pepper_report_indel_above_freq PEPPER_REPORT_INDEL_ABOVE_FREQ
                        Report all INDELs above frequency for re-genotyping even if the predicted value is low. Set 0 to disable this.

Arguments for Margin:
  --margin_haplotag_model MARGIN_HAPLOTAG_MODEL
                        Custom margin model.
  --margin_phase_model MARGIN_PHASE_MODEL
                        Custom margin model.

Arguments for DeepVariant:
  --dv_model DV_MODEL   Custom DeepVariant model. If this parameter is set, then one model will be used for both SNP and INDEL calling.
  --dv_model_snp DV_MODEL_SNP
                        Custom DeepVariant model for calling SNP calling.
  --dv_model_indel DV_MODEL_INDEL
                        Custom DeepVariant model for calling INDEL calling.
  --dv_alt_aligned_pileup DV_ALT_ALIGNED_PILEUP
                        DeepVariant alt_align_pileup used for make_examples associated model: dv_model.
  --dv_alt_aligned_pileup_snp DV_ALT_ALIGNED_PILEUP_SNP
                        DeepVariant alt_align_pileup used for make_examples for snp calling associated model: dv_model_snp.
  --dv_alt_aligned_pileup_indel DV_ALT_ALIGNED_PILEUP_INDEL
                        DeepVariant alt_align_pileup used for make_examples for indel calling associated model: dv_model_indel.
  --dv_pileup_image_width DV_PILEUP_IMAGE_WIDTH
                        DeepVariant image width.
  --dv_realign_reads DV_REALIGN_READS
                        If true then local read alingment will be performed. [set: true/false]
  --dv_partition_size DV_PARTITION_SIZE
                        DeepVariant partition_size used for make_examples.
  --dv_min_mapping_quality DV_MIN_MAPPING_QUALITY
                        DeepVariant minimum mapping quality.
  --dv_min_base_quality DV_MIN_BASE_QUALITY
                        DeepVariant minimum base quality.
  --dv_vsc_min_fraction_indels DV_VSC_MIN_FRACTION_INDELS
                        DeepVariant minimum vsc fraction for indels.
  --dv_vsc_min_fraction_snps DV_VSC_MIN_FRACTION_SNPS
                        DeepVariant minimum vsc fraction for snps.
  --dv_sort_by_haplotypes DV_SORT_BY_HAPLOTYPES
                        If true then haplotype sorting will be used. [set: true/false]
  --dv_parse_sam_aux_fields DV_PARSE_SAM_AUX_FIELDS
                        If true then auxiliary field parsing is enabled. [set: true/false]
  --dv_add_hp_channel DV_ADD_HP_CHANNEL
                        If true then hp channel will be added. [set: true/false]
  --dv_use_hp_information DV_USE_HP_INFORMATION
                        If true then hp information will be properly used. [set: true/false]
  --dv_use_multiallelic_mode DV_USE_MULTIALLELIC_MODE
                        If true multiallelic model will be used during post-processing. [set: true/false]
