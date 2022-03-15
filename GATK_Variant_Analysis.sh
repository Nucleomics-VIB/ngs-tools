#!/bin/env bash
# author:Stephane Plaisance (VIB-NC), 2019-12-12
# from: https://github.com/BITS-VIB/NGS-Variant-Analysis-training-2020

## Where are we?
workdir="/data/analyses/SRR17149558_Gallus"

## remap reads from SRR17149558 to gallus_gallus GRCg6a.105
cd ${workdir}

# create tmp folder
mkdir -p tmpfiles

# create date-tagged logs 
mkdir -p logs
log=logs/runlog.txt
cat /dev/null > ${log}

# record all script actions
#set -x
#exec > ${log} 2>&1
exec &> >(tee -i ${log})

# multithreading, adapt to your cpu count
bwathr=84
samtoolsthr=24

# how much RAM can java use? - set these to less than your available RAM!
javaopts="-Xms24g -Xmx24g"

# my samtools is here
samtools=$BIOTOOLS/samtools/bin/samtools

function mytest {
	"$@"
	local status=$?
	if [ ${status} -ne 0 ]; then
		echo "error with $1" >&2
		exit
	fi
	return ${status}
}

#############################################
# Get GATK Bundle files
#############################################

# Please get all necessary reference files before running that code


#############################################
# BWA index & MEM mapping
#############################################

## create a bwa index when absent
bwaidxfolder=bwa_index
mkdir -p ${bwaidxfolder}

reference_fa=reference/Gallus_gallus.GRCg6a.dna.toplevel.fa
bwaidx="${bwaidxfolder}/GRCg6a.105"

if [ ! -f ${bwaidxfolder}/index_created ]; then
	echo "# creating BWA index"
	bwa index ${reference_fa} -p ${bwaidx} \
  	&& touch ${bwaidxfolder}/index_created
else
	echo "# BWA index already exists"
fi

## map reads to reference
reads_1="reads/SRR17149558_1.fq.gz"
reads_2="reads/SRR17149558_2.fq.gz"

# mapping settings
thr=84

# edit in the command below
outpfx="SRR17149558"
samplename="SRR17149558"

outfolder=bwa_mappings
mkdir -p ${workdir}/${outfolder}

# map using BWA mem (only once)
if [ ! -f bwa_mappings/mapping_done ]; then
	echo "# mapping reads with BWA mem"
	cmd="bwa mem -t ${bwathr} \
		-M \
		-R '@RG\tID:SRR17149558\tLB:SRR17149558\tPU:HiSeq1500\tPL:Illumina\tSM:SRR17149558' \
		${bwaidx} \
		${reads_1} \
		${reads_2} | ${samtools} view -b - -o ${outfolder}/${outpfx}_rawmappings.bam"
	echo "# ${cmd}"
	eval ${cmd}
	# get flagstats
	# then flag record that mapping was done for next run
	${samtools} flagstat -@ ${samtoolsthr} \
		${outfolder}/${outpfx}_rawmappings.bam \
		> ${outfolder}/${outpfx}_rawmappings_flagstats.txt && \
		touch bwa_mappings/mapping_done
else
	echo "# BWA mapping already done"
fi


#############################################
# PICARD Cleanup & MarkDuplicates
#############################################

# more records in RAM speeds up when enough RAM is present
recinram=10000000

if [ ! -f bwa_mappings/MarkDuplicates_done ]; then
# sort by queryname
java ${javaopts} -jar $PICARD/picard.jar \
	SortSam \
	I=${outfolder}/${outpfx}_rawmappings.bam \
	O=${outfolder}/${outpfx}_rawmappings_qrysrt.bam \
	SO=queryname \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/

# mark duplicates
pixdist=100
java ${javaopts} -jar $PICARD/picard.jar \
	MarkDuplicates \
	I=${outfolder}/${outpfx}_rawmappings_qrysrt.bam \
	O=${outfolder}/${outpfx}_mrkdup.bam \
	M=${outfolder}/${outpfx}_MarkDuplicates.txt \
	ASO=queryname \
	REMOVE_DUPLICATES=false \
	OPTICAL_DUPLICATE_PIXEL_DISTANCE=${pixdist} \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/

# sort by coordinate and index
java ${javaopts} -jar $PICARD/picard.jar \
	SortSam \
	I=${outfolder}/${outpfx}_mrkdup.bam \
	O=${outfolder}/${outpfx}_mrkdup_srt.bam \
	SO=coordinate \
	CREATE_INDEX=true \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/

# fix tags
java ${javaopts} -jar $PICARD/picard.jar \
	SetNmMdAndUqTags \
	I=${outfolder}/${outpfx}_mrkdup_srt.bam \
	O=${outfolder}/${outpfx}_mrkdup_srt-tags.bam \
	R=${reference_fa} \
	CREATE_INDEX=true \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/

# validate final BAM content
java ${javaopts} -jar $PICARD/picard.jar \
	ValidateSamFile \
	I=${outfolder}/${outpfx}_mrkdup_srt-tags.bam \
	O=${outfolder}/${outpfx}_mrkdup_srt-tags.bam_ValidateSamFile.txt \
	R=${reference_fa} \
	M=SUMMARY \
	MO=100 \
	IGNORE_WARNINGS=FALSE \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/

else
	echo "# GATK MarkDuplicates already done"
fi

#############################################
# GATK4 BAM RECALIBRATION
#############################################

outfolder=gatk_preprocessing
mkdir -p ${outfolder}

bamfile=bwa_mappings/${outpfx}_mrkdup_srt-tags.bam
recalbamfile=${outpfx}_mrkdup_srt_recal.bam

# base quality score recalibration
knownsites="reference/gallus_gallus_snv_cor_dedup.vcf.gz"
knownindels="reference/gallus_gallus_indels_cor_dedup.vcf.gz"

if [ ! -f gatk_preprocessing/recalibration_done ]; then
# compute table before
java ${javaopts} -jar $GATK/gatk.jar \
	BaseRecalibrator \
	-I ${bamfile} \
	-R ${reference_fa} \
	--known-sites ${knownsites} \
	--known-sites ${knownindels} \
	-O ${outfolder}/recal_data.table \
	--tmp-dir tmpfiles/

# apply recalibration table
java ${javaopts} -jar $GATK/gatk.jar \
	ApplyBQSR \
	-I ${bamfile} \
	-R ${reference_fa} \
	-bqsr ${outfolder}/recal_data.table \
	-O ${outfolder}/${recalbamfile} \
	--interval-padding 100 \
	--add-output-sam-program-record \
	--use-original-qualities \
	--static-quantized-quals 10 \
	--static-quantized-quals 20 \
	--static-quantized-quals 30 \
	--tmp-dir tmpfiles/

# compute table after
java ${javaopts} -jar $GATK/gatk.jar \
	BaseRecalibrator \
	-I ${outfolder}/${recalbamfile} \
	-R ${reference_fa} \
	--known-sites ${knownsites} \
	--known-sites ${knownindels} \
	-O ${outfolder}/recal_data_after.table \
	--tmp-dir tmpfiles/

# create plots from both tables
java ${javaopts} -jar $GATK/gatk.jar \
	AnalyzeCovariates \
	-before ${outfolder}/recal_data.table \
	-after ${outfolder}/recal_data_after.table \
	-plots ${outfolder}/BQSR_report.pdf \
	-csv ${outfolder}/BQSR-report.csv \
	--tmp-dir tmpfiles/

# Picard CollectMultipleMetrics on final BAM
java ${javaopts} -jar $PICARD/picard.jar \
	CollectMultipleMetrics \
	I=${outfolder}/${recalbamfile} \
	O=${outfolder}/${outpfx}_multiple_metrics \
	R=${reference_fa} \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/ && \
	touch gatk_preprocessing/recalibration_done
else
	echo "# GATK BAM recalibration already done"
fi


#############################################
# GATK4 VARIANT CALLING
#############################################

outfolder=gatk_variantcalling
mkdir -p ${outfolder}

# dbSNP for ID-field annotation
dbsnp=reference/gallus_gallus_snv.vcf.gz

if [ ! -f gatk_variantcalling/calling_done ]; then
# call short variants to gvcf format and save supporting reads

# parallelize the pair hidden Markov models (pair HMM) process
hmmt=16

java ${javaopts} -jar $GATK/gatk.jar \
	HaplotypeCaller  \
	--input gatk_preprocessing/${recalbamfile} \
	--output ${outfolder}/${samplename}.g.vcf.gz \
	--reference ${reference_fa} \
	--emit-ref-confidence GVCF \
	--sample-ploidy 2 \
	--native-pair-hmm-threads ${hmmt} \
	--bam-output ${outfolder}/${samplename}_HC_aligned_reads.bam \
	--tmp-dir tmpfiles/ 
	
# convert to multi-VCF (merged samples)
# https://software.broadinstitute.org/gatk/documentation/article?id=11813
java ${javaopts} -jar $GATK/gatk.jar \
	GenotypeGVCFs \
	--reference ${reference_fa} \
	--variant ${outfolder}/${samplename}.g.vcf.gz \
	--output ${outfolder}/${samplename}.vcf.gz \
	--dbsnp ${dbsnp} \
	--use-new-qual-calculator \
	--tmp-dir tmpfiles/
	
# mark ExcessHet
threshold=54.69
java ${javaopts} -jar $GATK/gatk.jar \
	VariantFiltration \
	--variant ${outfolder}/${samplename}.vcf.gz \
	--filter-expression "ExcessHet > ${threshold}" \
	--filter-name "ExcessHet" \
	-O ${outfolder}/${samplename}_excesshet_filtered.vcf.gz \
	--tmp-dir tmpfiles/ && \
	touch gatk_variantcalling/calling_done
else
	echo "# GATK calling already done"
fi


#############################################
# GATK4 VARIANT RECALIBRATION
#############################################

outfolder=gatk_variantrecalibration
mkdir -p ${outfolder}

# recalibration sources
## variants-sets

# Known sites resource, not used in training: dbSNP
snv=reference/gallus_gallus_snv.vcf.gz
indels=reference/gallus_gallus_indels.vcf.gz

# maxgaussians default to 8 is too high for chr22-only calls
maxSNPgaussians=6
maxINDELgaussians=4

if [ ! -f gatk_variantrecalibration/variantrecalibation_done ]; then
# copy last files from previous step
cp gatk_variantcalling/${samplename}_excesshet_filtered.vcf.gz* ${outfolder}/

# extract a 6-column version of the data for recalibration
java ${javaopts} -jar $GATK/gatk.jar \
	MakeSitesOnlyVcf \
	-I ${outfolder}/${samplename}_excesshet_filtered.vcf.gz \
	-O ${outfolder}/${samplename}_excesshet_sitesonly.vcf.gz

# Build the SNP recalibration model
java ${javaopts} -jar $GATK/gatk.jar \
	VariantRecalibrator \
	-R ${reference_fa} \
	-V ${outfolder}/${samplename}_excesshet_sitesonly.vcf.gz \
	-O ${outfolder}/${samplename}_recalibrate_SNP.recal.vcf.gz \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${snv} \
	--trust-all-polymorphic \
	--use-annotation DP \
	--use-annotation QD \
	--use-annotation FS \
	--use-annotation SOR \
	--use-annotation MQ \
	--use-annotation MQRankSum \
	--use-annotation ReadPosRankSum \
	--mode SNP \
	--max-gaussians ${maxSNPgaussians} \
	--tranches-file ${outfolder}/${samplename}_recalibrate_snp.tranches \
	--tranche 100.0 \
	--tranche 99.95 \
	--tranche 99.9 \
	--tranche 99.8 \
	--tranche 99.6 \
	--tranche 99.5 \
	--tranche 99.4 \
	--tranche 99.3 \
	--tranche 99.0 \
	--tranche 98.0 \
	--tranche 97.0 \
	--tranche 90.0 \
	--rscript-file ${outfolder}/${samplename}_recalibrate_snp_plots.R \
	--tmp-dir tmpfiles/


# Build the Indel recalibration model
java ${javaopts} -jar $GATK/gatk.jar \
	VariantRecalibrator \
	-R ${reference_fa} \
	-V ${outfolder}/${samplename}_excesshet_sitesonly.vcf.gz \
	-O ${outfolder}/${samplename}_recalibrate_indels.recal.vcf.gz \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2 ${indels} \
	--trust-all-polymorphic \
	--use-annotation QD \
	--use-annotation DP \
	--use-annotation FS \
	--use-annotation SOR \
	--use-annotation MQRankSum \
	--use-annotation ReadPosRankSum \
	--mode INDEL \
	--max-gaussians ${maxINDELgaussians} \
	--tranches-file ${outfolder}/${samplename}_recalibrate_indels.tranches \
	--tranche 100.0 \
	--tranche 99.95 \
	--tranche 99.9 \
	--tranche 99.5 \
	--tranche 99.0 \
	--tranche 97.0 \
	--tranche 96.0 \
	--tranche 95.0 \
	--tranche 94.0 \
	--tranche 93.5 \
	--tranche 93.0 \
	--tranche 92.0 \
	--tranche 91.0 \
	--tranche 90.0 \
	--rscript-file ${outfolder}/${samplename}_recalibrate_indels_plots.R \
	--tmp-dir tmpfiles/


# Apply the desired level of recalibration to the SNPs in the call set
java ${javaopts} -jar $GATK/gatk.jar \
	ApplyVQSR \
	-R ${reference_fa} \
	-V ${outfolder}/${samplename}_excesshet_filtered.vcf.gz \
	-O ${outfolder}/${samplename}_recalibrated_snps_raw_indels.vcf.gz \
	--mode SNP \
	--recal-file ${outfolder}/${samplename}_recalibrate_SNP.recal.vcf.gz \
	--tranches-file ${outfolder}/${samplename}_recalibrate_snp.tranches \
	--truth-sensitivity-filter-level 99.0 \
	--tmp-dir tmpfiles/

# Apply the desired level of recalibration to the Indels in the call set
java ${javaopts} -jar $GATK/gatk.jar \
	ApplyVQSR \
	-R ${reference_fa} \
	-V ${outfolder}/${samplename}_recalibrated_snps_raw_indels.vcf.gz \
	-O ${outfolder}/${samplename}_VQSR.vcf.gz \
	-mode INDEL \
	--recal-file ${outfolder}/${samplename}_recalibrate_indels.recal.vcf.gz \
	--tranches-file ${outfolder}/${samplename}_recalibrate_indels.tranches \
	--truth-sensitivity-filter-level 99.7 \
	--create-output-variant-index true \
	--tmp-dir tmpfiles/ && \
	touch gatk_variantrecalibration/variantrecalibation_done
else
	echo "# GATK variant recalibration already done"
fi


#############################################
# SnpEff on VQSR
#############################################

outfolder="snpeff"
mkdir -p ${outfolder}

build="GRCg6a.105"

if [ ! -f snpeff/VQSR_annotation_done ]; then
java ${javaopts} -jar $SNPEFF/snpEff.jar \
	-htmlStats ${outfolder}/gatk_recal_snpEff_summary.html \
	${build} \
	gatk_variantrecalibration/${samplename}_VQSR.vcf.gz | \
	bgzip -c > ${outfolder}/${samplename}_VQSR_snpeff.vcf.gz && \
	tabix -p vcf ${outfolder}/${samplename}_VQSR_snpeff.vcf.gz && \
	touch snpeff/VQSR_annotation_done
else
	echo "# SNPEff annotation for VQSR variant already done"
fi


#############################################
# GATK4 VARIANT HARD FILTERING
#############################################

# instructions and filters from:
# https://gatkforums.broadinstitute.org/gatk/discussion/23216/how-to-filter-variants-either-with-vqsr-or-by-hard-filtering
# Genepattern: QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0 || QUAL < 30

outfolder=gatk_varianthardfiltering
mkdir -p ${outfolder}

if [ ! -f gatk_varianthardfiltering/hardFiltering_recalibation_done ]; then

# copy the GenotypeGVCFs vcf output with marked excesshet output from above
# create local copy of previous file
cp gatk_variantcalling/${samplename}_excesshet_filtered.vcf.gz* ${outfolder}/


###########################################
# 2) hard-Filter SNPs on multiple metrics
###########################################

# produces a VCF with records with SNP-type variants only.
java ${javaopts} -jar $GATK/gatk.jar \
	SelectVariants \
	-V ${outfolder}/${samplename}_excesshet_filtered.vcf.gz \
	--select-type-to-include SNP \
	-O ${outfolder}/${samplename}_snp.vcf.gz

# This produces a VCF with the same variant records now annotated with filter status. Specifically, if a record passes all the filters, it receives a PASS label in the FILTER column.
# A record that fails a filter #receives the filter name in the FILTER column, e.g. SOR3.
# If a record fails multiple filters, then each failing filter name appears in the FILTER column separated by semi-colons ; e.g. "MQRankSum-12.5;ReadPosRankSum-8".

java ${javaopts} -jar $GATK/gatk.jar \
	VariantFiltration \
	-V ${outfolder}/${samplename}_snp.vcf.gz \
	-O ${outfolder}/${samplename}_snp_filtered.vcf.gz \
	--filter "QD < 2.0" --filter-name "QD2" \
	--filter "QUAL < 30.0" --filter-name "QUAL30" \
	--filter "SOR > 3.0" --filter-name "SOR3" \
	--filter "FS > 60.0" --filter-name "FS60" \
	--filter "MQ < 40.0" --filter-name "MQ40" \
	--filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
	--filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"


######################################################
# 3) hard-Filter INDELs and MIXED on multiple metrics
######################################################

# produces a VCF with records with INDEL-type variants only.
java ${javaopts} -jar $GATK/gatk.jar \
	SelectVariants \
	-V ${outfolder}/${samplename}_excesshet_filtered.vcf.gz \
	-O ${outfolder}/${samplename}_mixed_indels.vcf.gz \
	--select-type-to-include INDEL \
	--select-type-to-include MIXED

java ${javaopts} -jar $GATK/gatk.jar \
	VariantFiltration \
	-V ${outfolder}/${samplename}_mixed_indels.vcf.gz \
	-O ${outfolder}/${samplename}_mixed_indels_filtered.vcf.gz \
	--filter "QD < 2.0" --filter-name "QD2" \
	--filter "QUAL < 30.0" --filter-name "QUAL30" \
	--filter "FS > 200.0" --filter-name "FS200" \
	--filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"


###########################################
## 4) merge SNP and Indel filtered calls	
###########################################

# combine filtered into one file using Picard
java ${javaopts} -jar $GATK/gatk.jar \
	MergeVcfs \
	-I ${outfolder}/${samplename}_snp_filtered.vcf.gz \
	-I ${outfolder}/${samplename}_mixed_indels_filtered.vcf.gz \
	-R ${reference_fa} \
	-O ${outfolder}/${samplename}_snp_indel_filtered.vcf.gz && \
	touch ${outfolder}/hardFiltering_recalibation_done
else
	echo "# GATK hardFiltering recalibration already done"
fi


#############################################
# SnpEff on HardFiltering
#############################################

outfolder="snpeff"
mkdir -p ${outfolder}

build="GRCg6a.105"

if [ ! -f snpeff/hardFiltering_annotation_done ]; then
java ${javaopts} -jar $SNPEFF/snpEff.jar \
	-htmlStats ${outfolder}/gatk_hardfiltering_snpEff_summary.html \
	${build} \
	gatk_varianthardfiltering/${samplename}_snp_indel_filtered.vcf.gz | \
	bgzip -c > ${outfolder}/${samplename}_hardfiltering_snpeff.vcf.gz && \
	tabix -p vcf ${outfolder}/${samplename}_hardfiltering_snpeff.vcf.gz && \
	touch snpeff/hardFiltering_annotation_done
else
	echo "# SNPEff annotation for hard-filtered variant already done"
fi

# cleanup leftovers
rm tmpfiles/*
