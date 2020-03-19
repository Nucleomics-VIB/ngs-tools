#!/bin/bash
# script name: run_quast_nogenepred.sh 
#
# Stephane Plaisance (VIB-NC) 2020/03/16; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# requires quast installed and running
# required once: create a conda env to install the required apps
# adapt next line to point to the right conda.sh init script
# see conda activate script for details
source /etc/profile.d/conda.sh
conda activate atwork3 || \
  ( echo "# the conda environment 'atwork3' adding BUSCO was not found on this machine" ;
    echo "# please read the top part of the script!" \
    && exit 1 )

# put all assemblies in a folder named <assemblies>
# put the reference and gff3 annotations in a folder named <reference>
# put all ONT plain or gzipped reads in a folder named <ont_reads>

longreads="nanopore"
#longreads="pacbio"

outfolder="quast_results_large" 
input=( $(find assemblies -name "*.f*a" | sort -h) )

# create labels and remove trailing comma
labels=$(find assemblies -name "*.f*a" | sort -h | \
	sed -e 's/assemblies\///g' | tr "\n" "," | sed 's/.$//')

# format read inputs (if present)
readfolder="ont_reads"
declare -a reada
if [ -d "${readfolder}" ]; then
	reads=( $(find ${readfolder} -name "*.f*" | sort) )
	reada=( ${reads[@]/#/"--${longreads} "} )
fi

refg=$(find reference -name "*.f*a")
annot=$(find reference -name "*.gff?")

thr=24

# add gene predictions
#        --gene-finding \
#        --fungus \
#        --conserved-genes-finding \
#        --rna-finding \

# run with large for genomes >100Mb
cmd="quast.py -t ${thr} \
        --fragmented \
        --large \
        -r ${refg} \
        -g gene:${annot} \
        --min-identity 95 \
        -l ${labels} \
        --circos \
        -o ${outfolder} \
        ${reada[@]} \
        ${input[@]}"

echo "# ${cmd}"
eval ${cmd}

# return to base
conda deactivate
