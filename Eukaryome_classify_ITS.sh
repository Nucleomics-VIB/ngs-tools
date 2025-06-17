#!/bin/bash

# script: Eukaryome_classify_ITS.sh
# Classify OTU representative sequences with the Eukaryome reference DB
# modifies the script Unite_classify_ITS.sh to use the Eukaryome classifier
#
# Stephane Plaisance (VIB-NC) 2020/10/14; v1.0
# modified by SP@NC: 2025/06/17; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

version="1.0, 2025_06_17"

usage='# Usage: Eukaryome_classify_ITS.sh -i <rep-seqs.qza> 
# script version '${version}'
# [optional: -o <result_prefix|sample_name_Eukaryome_v1.9.4>]
# [optional: -t <threads|8>]
# [optional: -d <classifier_directory|/data/biodata/qiime2_ITS>]
# [optional: -c <classifier_filename|General_EUK_longread_v1.9.4_classifier.qza>]
# [optional: -v <classifier_version|v1.9.4>]
# [optional: -h <this help text>]'

while getopts "i:o:t:d:c:v:h" opt; do
  case $opt in
    i) infile=${OPTARG};;
    o) outpref=${OPTARG};;
    t) threads=${OPTARG};;
    d) EUK_DIR=${OPTARG};;
    c) eukaryomeClassifier=${OPTARG};;
    v) Version=${OPTARG};;
    h) echo "${usage}" >&2; exit 0;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1;;
    *) echo "Option -${OPTARG} requires an argument." >&2; exit 1;;
  esac
done

# Defaults
threads=${threads:-8}
EUK_DIR=${EUK_DIR:-"/data/biodata/qiime2_ITS"}
Version=${Version:-"v1.9.4"}
eukaryomeClassifier=${eukaryomeClassifier:-"General_EUK_longread_v1.9.4_classifier.qza"}

# Check if input file was provided
if [ -z "${infile}" ]; then
  echo "${usage}" >&2
  echo "# provide a rep-seqs.qza file as input" >&2
  exit 1
fi

# Create output prefix if not provided
if [ -z "${outpref}" ]; then
  # remove .qza suffix if present in input name
  outpref=${infile%.qza}
  # add the db type to the output name
  outpref=${outpref}_Eukaryome_${Version}
fi

# where the Eukaryome classifier is stored
dbpath="${EUK_DIR}/${eukaryomeClassifier}"

# check if classifier exists at specified path
if [ ! -f "${dbpath}" ]; then
  echo "# The classifier was not found at: ${dbpath}"
  echo "# check the path and version!"
  exit 1
fi

# run the classifier
echo "# classifying with Eukaryome ${Version}"
date

qiime feature-classifier classify-sklearn \
  --i-classifier ${dbpath} \
  --i-reads ${infile} \
  --p-n-jobs ${threads} \
  --o-classification ${outpref}-taxonomy.qza

# create visualization of the taxonomy
qiime metadata tabulate \
  --m-input-file ${outpref}-taxonomy.qza \
  --o-visualization ${outpref}-taxonomy.qzv

date
echo "# all done"