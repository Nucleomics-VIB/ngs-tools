#!/bin/bash

# script: getCosmic.sh
# Aim: download Cosmic files except the large database dump
#
# requires Sanger registration
# requires jq to decode jason
#
# Stéphane Plaisance - VIB-Nucleomics Core - 2021-01-29 v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

read -p 'Sanger login email [stephane.plaisance@vib.be]: ' login
login=${login:-"stephane.plaisance@vib.be"}

read -sp 'Sanger Password: ' password

read -p 'Genome build: [GRCh38]: ' build
build=${build:-"GRCh38"}

token=$(echo "${login}:${password}" | base64)

# base URL for Cosmic download
baseurl="https://cancer.sanger.ac.uk/cosmic/file_download"
url="${baseurl}/${build}/cosmic/latest"

while read file; do 

filename=$(basename ${file})

# do not get the large database dump
if [ "${filename}" == "COSMIC_ORACLE_EXPORT.dmp.gz.tar" ]; then
echo "# not downloading the DB dump: ${filename}"
continue
fi

echo
echo "# downloading: ${filename}"

# get path
res=$(curl -s -H "Authorization: Basic ${token}" ${baseurl}/${file})
dlurl=$(echo ${res} | jq '.url')
# echo ${dlurl}

cmd="curl -o ${filename} ${dlurl}"
# echo "# ${cmd}"
eval ${cmd}

done < <(curl -s ${url} | jq -r '.[]')

exit 0

## "GRCh38/cosmic/v92/All_COSMIC_Genes.fasta.gz"
# NOT "GRCh38/cosmic/v92/COSMIC_ORACLE_EXPORT.dmp.gz.tar"
## "GRCh38/cosmic/v92/Cancer_Gene_Census_Hallmarks_Of_Cancer.tsv.gz"
## "GRCh38/cosmic/v92/CosmicBreakpointsExport.tsv.gz"
## "GRCh38/cosmic/v92/CosmicCompleteCNA.tsv.gz"
## "GRCh38/cosmic/v92/CosmicCompleteDifferentialMethylation.tsv.gz"
## "GRCh38/cosmic/v92/CosmicCompleteGeneExpression.tsv.gz"
## "GRCh38/cosmic/v92/CosmicCompleteTargetedScreensMutantExport.tsv.gz"
## "GRCh38/cosmic/v92/CosmicFusionExport.tsv.gz"
## "GRCh38/cosmic/v92/CosmicGenomeScreensMutantExport.tsv.gz"
## "GRCh38/cosmic/v92/CosmicHGNC.tsv.gz"
## "GRCh38/cosmic/v92/CosmicMutantExport.tsv.gz"
## "GRCh38/cosmic/v92/CosmicMutantExportCensus.tsv.gz"
## "GRCh38/cosmic/v92/CosmicMutationTracking.tsv.gz"
## "GRCh38/cosmic/v92/CosmicNCV.tsv.gz"
## "GRCh38/cosmic/v92/CosmicResistanceMutations.tsv.gz"
## "GRCh38/cosmic/v92/CosmicSample.tsv.gz"
## "GRCh38/cosmic/v92/CosmicStructExport.tsv.gz"
## "GRCh38/cosmic/v92/CosmicTranscripts.tsv.gz"
## "GRCh38/cosmic/v92/OracleSchemaDocumentation.pdf"
## "GRCh38/cosmic/v92/README-cosmic.txt"
## "GRCh38/cosmic/v92/VCF/CosmicCodingMuts.normal.vcf.gz"
## "GRCh38/cosmic/v92/VCF/CosmicCodingMuts.vcf.gz"
## "GRCh38/cosmic/v92/VCF/CosmicCodingMuts.vcf.gz.tbi"
## "GRCh38/cosmic/v92/VCF/CosmicNonCodingVariants.normal.vcf.gz"
## "GRCh38/cosmic/v92/VCF/CosmicNonCodingVariants.vcf.gz"
## "GRCh38/cosmic/v92/VCF/CosmicNonCodingVariants.vcf.gz.tbi"
## "GRCh38/cosmic/v92/ascat_acf_ploidy.tsv"
## "GRCh38/cosmic/v92/cancer_gene_census.csv"
## "GRCh38/cosmic/v92/classification.csv"
