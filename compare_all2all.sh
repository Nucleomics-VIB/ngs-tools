#!/bin/bash

# compare_all2all.sh produce a all to all assemblies comparison using WUB tools
#
# Requirements:
# run on a unix computer
# !! wub installed and virtual environment configured
# find it at: https://github.com/nanoporetech/wub
# two or more related fasta assemblies to be compared in the current folder
#
# Stephane Plaisance (VIB-NC+BITS) 2018/04/18; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# check parameters for your system
version="1.0, 2018/04/18"

# path to the Assemblytics scripts
default_path_to_wub="/opt/biotools/virtualenv/"

# activate wub virtenv (without this the script will not work!!)
. ${default_path_to_wub}/wub_env/bin/activate

# create timestamp for multiple runs
timestamp=$(date +%s)

# write all outputs to log file
exec > >(tee -i "compare_all2all_logfile_${timestamp}.txt")
exec 2>&1

# lastal defaults and 24 threads (see `lastal -h`)
lastal_opts="a:7,b:1,P:24"

# loop twice for all vs all results
for ref in *.fasta; do
        for qry in *.fasta; do
                echo "###############################################################################"
                echo "# comparing ${qry} to ${ref}"

                echo
                echo "# running dnadiff"
                cmd="compare_genomes_dnadiff.py -k \
                        -r ${qry%.fasta}_vs_${ref%.fasta}_dnadiff_results.txt \
                        ${ref} \
                        ${qry}"
                echo "# $cmd"
                eval $cmd

                echo
                echo "# running lastal"
                cmd="compare_genomes_lastal.py \
                        -r ${qry%.fasta}_vs_${ref%.fasta}_lastal_report.pdf \
                        -l ${lastal_opts} \
                        ${ref} \
                        ${qry}"
                echo "# $cmd"
                eval $cmd
        done
        echo "###############################################################################"
        echo
done
