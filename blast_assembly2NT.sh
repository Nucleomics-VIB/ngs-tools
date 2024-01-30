#!/bin/bash

# script: blast_assembly2NT.sh
# Compare denovo assembly contigs/scaffolds to NCBI NT
# find major match for each contig and report hit counts
# add acc title and taxonomy info to identify/document each contig
#
# SP@NC 2024-01-29; v1.0

# custom function to activate a conda env present on the server
activate_conda_env() {
  local myenv="$1"
  source /etc/profile.d/conda.sh
  conda activate "${myenv}" || \
  (echo "# the conda environment ${myenv} was not found on this machine" ;
    echo "# please read the top part of the script!" \
    && exit 1
  )
}

workdir=$PWD

outfolder=${workdir}/blast_results
mkdir -p ${outfolder}/logs

queryfasta=${1}
pfx=$(basename ${queryfasta%.fa*})

maxtarget=10

# RAM limits
nthr1=44
nthr2=24
nthr3=8

################################
# blastn assembly contigs to nt
################################

if [ ! -f "${outfolder}/${pfx}_blast_summary.txt" ]; then

# unset previous content between loops
unset genus
unset species
unset refid

myenv=blast
activate_conda_env "${myenv}"

blastn \
  -evalue 1e-100 \
  -db nt \
  -query ${queryfasta} \
  -num_threads ${nthr1} \
  -max_target_seqs ${maxtarget} \
  -outfmt 7 \
  -out ${outfolder}/${pfx}_blastn.txt \
  2> ${outfolder}/logs/${pfx}_blastn_log.txt

 sort the list and count
cut -f 1,2 ${outfolder}/${pfx}_blastn.txt \
  | grep -v "^#" \
  | uniq \
  > ${outfolder}/${pfx}_blast_counts.txt

conda deactivate

#######################
# identify blastn hits
#######################

myenv=entrez-direct
activate_conda_env "${myenv}"

# list unique NCBI acc.version frmo all blast results
cut -f 2 ${outfolder}/${pfx}_blast_counts.txt \
  | sort \
  | uniq > ${outfolder}/${pfx}_acc_list.txt

# add NCBI hit Title and taxonomy to hits
efetch -input ${outfolder}/${pfx}_acc_list.txt \
  -db nucleotide \
  -format docsum \
  | xtract -pattern DocumentSummary -element AccessionVersion,Title,Organism \
  > ${outfolder}/${pfx}_acc2info.txt

# fill arrays for retreival
declare -A acc2title
declare -A acc2tax

while IFS=$'\t' read -r -a myArray
do
  acc2title["${myArray[0]}"]="${myArray[1]}"
  acc2tax["${myArray[0]}"]="${myArray[2]}"
done < ${outfolder}/${pfx}_acc2info.txt

# parse original blast results and add columns
while read ctg id; do

title=${acc2title["${id}"]:-"na"}
chrom="$(echo ${title} | sed -r 's/.*chromosome[^0-9]*([0-9]+).*/\1/')"
taxonomy=${acc2tax["${id}"]:-"na"}

# Check if the chrom is an integer
if [[ $chrom =~ ^-?[0-9]+$ ]]; then
    chr=${chrom}
else
   chr="na"
fi

printf "%s\t%s\t%s\t%s\t%s\n" "${ctg}" "${id}" "${title}" "chr_${chr}" "${taxonomy}"

done < ${outfolder}/${pfx}_blast_counts.txt \
  > ${outfolder}/${pfx}_blast_summary.txt

conda deactivate

####################################################
# extract tophit for genus and species from summary
####################################################

# sort decreasing and keep first as candidate
tophit=( $(cut -f 5 ${outfolder}/${pfx}_blast_summary.txt \
  | sort \
  | uniq -c \
  | sort -k1,1nr \
  | sed 's/  */ /g' \
  | cut -d " " -f 3,4 \
  | head -1) )

# set variables for later use
genus=${tophit[0]}
species=${tophit[1]%%[[:space:]]}

# debug
if [ -z "${genus}" ] || [ -z "${species}" ]; then 
echo "# failed to identify a tophit: ${tophit[@]}"
exit 1
fi

echo ${tophit[@]} > ${outfolder}/${pfx}_tophit.txt

# end if ${outfolder}/${pfx}_blast_summary.txt not present
fi
