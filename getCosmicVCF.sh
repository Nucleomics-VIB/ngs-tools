#!/bin/bash

# script: getCosmicVCF.sh
# Aim: download Cosmic VCF and info files only
#
# requires Sanger registration to access server
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
version=""

while read file; do 

filename=$(basename ${file})
version=$(echo "${file}" | sed "s|.*/v\(.*\)/.*|\1|")

# do not get the large database dump and other large files
if [ ${filename: -7} != ".vcf.gz" ] && \
   [ "${filename}" != "README-cosmic.txt" ] && \
   [ "${filename}" != "ascat_acf_ploidy.tsv" ] && \
   [ "${filename}" != "cancer_gene_census.csv" ] && \
   [ "${filename}" != "CosmicHGNC.tsv.gz" ] && \
   [ "${filename}" != "classification.csv" ]; then
echo "# file not downloaded: ${filename}"
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

# add flag to folder
touch "Cosmic-build_${build}-v${version}-$(date +%F)"

exit 0

# only files with '#+' are downloaded
## "GRCh38/cosmic/v92/All_COSMIC_Genes.fasta.gz"
## "GRCh38/cosmic/v92/COSMIC_ORACLE_EXPORT.dmp.gz.tar"
## "GRCh38/cosmic/v92/Cancer_Gene_Census_Hallmarks_Of_Cancer.tsv.gz"
## "GRCh38/cosmic/v92/CosmicBreakpointsExport.tsv.gz"
## "GRCh38/cosmic/v92/CosmicCompleteCNA.tsv.gz"
## "GRCh38/cosmic/v92/CosmicCompleteDifferentialMethylation.tsv.gz"
## "GRCh38/cosmic/v92/CosmicCompleteGeneExpression.tsv.gz"
## "GRCh38/cosmic/v92/CosmicCompleteTargetedScreensMutantExport.tsv.gz"
## "GRCh38/cosmic/v92/CosmicFusionExport.tsv.gz"
## "GRCh38/cosmic/v92/CosmicGenomeScreensMutantExport.tsv.gz"
#+ "GRCh38/cosmic/v92/CosmicHGNC.tsv.gz"
## "GRCh38/cosmic/v92/CosmicMutantExport.tsv.gz"
## "GRCh38/cosmic/v92/CosmicMutantExportCensus.tsv.gz"
## "GRCh38/cosmic/v92/CosmicMutationTracking.tsv.gz"
## "GRCh38/cosmic/v92/CosmicNCV.tsv.gz"
## "GRCh38/cosmic/v92/CosmicResistanceMutations.tsv.gz"
## "GRCh38/cosmic/v92/CosmicSample.tsv.gz"
## "GRCh38/cosmic/v92/CosmicStructExport.tsv.gz"
## "GRCh38/cosmic/v92/CosmicTranscripts.tsv.gz"
## "GRCh38/cosmic/v92/OracleSchemaDocumentation.pdf"
#+ "GRCh38/cosmic/v92/README-cosmic.txt"
#+ "GRCh38/cosmic/v92/VCF/CosmicCodingMuts.normal.vcf.gz"
#+ "GRCh38/cosmic/v92/VCF/CosmicCodingMuts.vcf.gz"
## "GRCh38/cosmic/v92/VCF/CosmicCodingMuts.vcf.gz.tbi"
#+ "GRCh38/cosmic/v92/VCF/CosmicNonCodingVariants.normal.vcf.gz"
#+ "GRCh38/cosmic/v92/VCF/CosmicNonCodingVariants.vcf.gz"
## "GRCh38/cosmic/v92/VCF/CosmicNonCodingVariants.vcf.gz.tbi"
## "GRCh38/cosmic/v92/ascat_acf_ploidy.tsv"
#+ "GRCh38/cosmic/v92/cancer_gene_census.csv"
#+ "GRCh38/cosmic/v92/classification.csv"
