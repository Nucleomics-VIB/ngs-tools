#!/bin/bash

# run funAnnotate on a denovo assembly
# https://funannotate.readthedocs.io/en/latest/docker.html
#
# Requirements:
# you will be asked your sudo password (no sudo no run!)
# run on a unix computer installed with
#
# docker image: nextgenusfs/funannotate
# download/pull the image from docker hub
# $ sudo docker pull nextgenusfs/funannotate
# download bash wrapper script (optional)
# $ wget -O funannotate-docker \
#  https://raw.githubusercontent.com/nextgenusfs/funannotate/master/funannotate-docker
#
# might need to make this executable on your system
# $ chmod +x /path/to/funannotate-docker
#
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
version="1.0, 2022_06_03"

usage='# Usage: run_funAnnotate.sh -i <input assembly>
# -b busco augustus training set (get the list from a dry run with -l)
# [optional: -l get list of available busco traininng sets and stop]
# [optional: -s <Sample Species name|Species>]
# [optional: -S <Sample strain name|strain>]
# [optional: -n <output folder|fun>]
# [optional: -t <threads|4>]
# script version '${version}

# edit to your own path
wrapper=$(which funannotate-docker)

while getopts "i:b:s:S:n:t:lh" opt; do
  case $opt in
    i) opt_asm=${OPTARG} ;;
    b) opt_busco=${OPTARG} ;;
    s) opt_species=${OPTARG} ;;
    S) opt_strain=${OPTARG} ;;
    n) opt_name=${OPTARG} ;;
    l) echo "# list of augustus species";
       sudo ${wrapper} species;
       exit 0 ;;
    t) opt_thr=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2;
       exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; 
       exit 1 ;;
  esac
done

# check if all dependencies are present
declare -a arr=("interproscan.sh" "${wrapper}")
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || \
  ( echo "# required ${prog} not found in PATH"
exit 1 )
done

# check for the docker image
dckr=$(sudo docker images | grep nextgenusfs/funannotate | cut -d " " -f 1)
if [ ! ${dckr} == "nextgenusfs/funannotate" ]; then
  (echo "docker image not found"; exit 1)
fi

# check for the eggnog conda environment
eggn=$(conda env list | grep eggnog | cut -d " " -f 1)
if [ ! ${eggn} == "eggnog-mapper" ]; then
  (echo "eggnog env not found"; exit 1)
fi

# test if minimal arguments were provided
if [ -z "${opt_asm}" ]; then
   echo "# no assembly input provided!"
   echo "${usage}"
   exit 1
fi

if [ ! -f "${opt_asm}" ]; then
    echo "${opt_asm} file not found!";
    exit 1
fi

# test is user provided training set exists
busc=$(sudo ${wrapper} species | grep ${opt_busco})
if [[ ! ${busc} == ${opt_busc}* ]]; then
  (echo "busco training set not found, try running with -l"; exit 1)
fi

# set defaults
pfx=$(basename ${opt_asm%.fa*})
outfolder=${opt_name:-"fun"}
nthr=${opt_thr:-4}
mskthr=12
minlen=1000

########################
# funAnnotate Pipeline #
########################

# prepare assembly file
echo
echo "# preparing assembly file"
sudo ${wrapper} clean -i ${opt_asm} --minlen ${minlen} -o ${pfx}.genome.cleaned.fa
sudo ${wrapper} sort -i ${pfx}.genome.cleaned.fa -b scaffold -o ${pfx}.genome.cleaned.sorted.fa
sudo ${wrapper} mask -i ${pfx}.genome.cleaned.sorted.fa --cpus ${mskthr} -o Assembly.fa

# predict
echo
echo "# running funannotate predict"
time sudo ${wrapper} predict \
  -i Assembly.fa \
  -o ${outfolder} \
  --species ${opt_species} \
  --strain ${opt_strain} \
  --busco_seed_species ${opt_busco} \
  --cpus ${nthr}

# compute accessory annotations

# interproscan
mkdir -p iprscan_out

# get protein output from the prediction folder
protres=$(find ${outfolder}/predict_results -name '*.proteins.fa')

echo
echo "# running interproscan"
time interproscan.sh \
  -i ${protres} \
  -f gff3,xml,tsv \
  -dp \
  -cpu ${nthr} \
  -d iprscan_out

# eggnog-mapper
myenv=${eggn}
source /etc/profile.d/conda.sh
conda activate ${myenv} || \
  ( echo "# the conda environment ${myenv} was not found on this machine" ;
    echo "# please read the top part of the script!" \
    && exit 1 )

mkdir -p eggnog_out

echo
echo "# running eggnog mapper"
time emapper.py \
  -i ${protres} \
  --itype proteins \
  --output_dir eggnog_out \
  -o eggnog \
  --temp_dir . \
  --scratch_dir . \
  --cpu ${nthr} \
  --excel

conda deactivate

# final annotation
iprres=$(find iprscan_out -name '*.proteins.fa.xml')
eggres=eggnog_out/eggnog.emapper.annotations

echo
echo "# running funannotate annotate"
time sudo ${wrapper} annotate \
  -i ${outfolder} \
  --iprscan ${iprres} \
  --eggnog ${eggres} \
  --cpus ${nthr}

echo
echo "# copying final results to local folder"
cp -r ${outfolder}/annotate_results ./funannotate_results

exit 0

######################################################
################### Command info #####################
######################################################

https://funannotate.readthedocs.io/en/latest/commands.html

# $ funannotate
# 
#     Usage:       funannotate <command> <arguments>
#     version:     1.7.0
# 
#     Description: Funannotate is a genome prediction, annotation, and comparison pipeline.
# 
#     Commands:
#       clean       Find/remove small repetitive contigs
#       sort        Sort by size and rename contig headers
#       mask        Repeatmask genome assembly
# 
#       train       RNA-seq mediated training of Augustus/GeneMark
#       predict     Run gene prediction pipeline
#       fix         Fix annotation errors (generate new GenBank file)
#       update      RNA-seq/PASA mediated gene model refinement
#       remote      Partial functional annotation using remote servers
#       iprscan     InterProScan5 search (Docker or local)
#       annotate    Assign functional annotation to gene predictions
#       compare     Compare funannotated genomes
# 
#       util        Format conversion and misc utilities
#       setup       Setup/Install databases
#       test        Download/Run funannotate installation tests
#       check       Check Python, Perl, and External dependencies [--show-versions]
#       species     list pre-trained Augustus species
#       database    Manage databases
#       outgroups   Manage outgroups for funannotate compare
# 
#     Written by Jon Palmer (2016-2019) nextgenusfs@gmail.com
# 
# 
# $ funannotate clean
# 
#     Usage:       funannotate clean <arguments>
#     version:     1.7.0
# 
# 	Description: The script sorts contigs by size, starting with shortest contigs it uses minimap2
# 							 to find contigs duplicated elsewhere, and then removes duplicated contigs.
# 
# 	Arguments:
# 	  -i, --input    Multi-fasta genome file (Required)
# 	  -o, --out      Cleaned multi-fasta output file (Required)
# 	  -p, --pident   Percent identity of overlap. Default = 95
# 	  -c, --cov      Percent coverage of overlap. Default = 95
# 	  -m, --minlen   Minimum length of contig to keep. Default = 500
# 	  --exhaustive   Test every contig. Default is to stop at N50 value.
# 
# 
# $ funannotate sort
# 
# 	Usage:       funannotate sort <arguments>
# 	version:     1.7.0
# 
# 	Description: This script sorts the input contigs by size (longest->shortest) and then relabels
# 							 the contigs with a simple name (e.g. scaffold_1).  Augustus can have problems with
# 							 some complicated contig names.
# 
# 	Arguments:
# 	  -i, --input    Multi-fasta genome file. (Required)
# 	  -o, --out      Sorted by size and relabeled output file. (Required)
# 	  -b, --base     Base name to relabel contigs. Default: scaffold
# 	  --minlen       Shorter contigs are discarded. Default: 0
# 
# 
# $ funannotate mask
# 	  
# 	Usage:       funannotate mask <arguments>
# 	version:     1.7.0
# 
# 	Description: This script is a wrapper for repeat masking. Default is to run very simple
# 							 repeat masking with tantan. The script can also run RepeatMasker and/or
# 							 RepeatModeler. It will generate a softmasked genome. Tantan is probably not
# 							 sufficient for soft-masking an assembly, but with RepBase no longer being
# 							 available RepeatMasker/Modeler may not be functional for many users.
# 
# 	Arguments:
# 	  -i, --input                    Multi-FASTA genome file. (Required)
# 	  -o, --out                      Output softmasked FASTA file. (Required)
# 
# 	Optional:
# 	  -m, --method                   Method to use. Default: tantan [repeatmasker, repeatmodeler]
# 	  -s, --repeatmasker_species     Species to use for RepeatMasker
# 	  -l, --repeatmodeler_lib        Custom repeat database (FASTA format)
# 	  --cpus                         Number of cpus to use. Default: 2
# 	  --debug                        Keep intermediate files