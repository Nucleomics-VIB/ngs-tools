#!/bin/bash

# script: blast_assembly2NT.sh
# Compare denovo assembly contigs/scaffolds to NCBI NT
# find major match for each contig and report hit counts
# add acc title, tlen, aliw and taxonomy info to identify/document each contig
#
# SP@NC 2024-01-29; v1.01

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

# sort the list and count
cut -f 1,2 ${outfolder}/${pfx}_blastn.txt \
  | grep -v "^#" \
  | uniq \
  > ${outfolder}/${pfx}_blast_counts.txt

conda deactivate

#################################
# compute target alignment width
#################################

blast2targetwidth.py ${outfolder}/${pfx}_blast_counts.txt \
  > ${outfolder}/${pfx}_blastn_targetwidth.txt

# put results into an array for later retrieval
declare -A acc2aliw

while IFS=$'\t' read -r -a myArray
do
  acc2aliw["${myArray[0]}"]="${myArray[1]}"
done < ${outfolder}/${pfx}_blastn_targetwidth.txt

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
  | xtract -pattern DocumentSummary -element AccessionVersion,Title,Slen,Organism \
  > ${outfolder}/${pfx}_acc2info.txt

# fill arrays for retreival
declare -A acc2title
declare -A acc2len
declare -A acc2tax

while IFS=$'\t' read -r -a myArray
do
  acc2title["${myArray[0]}"]="${myArray[1]}"
  acc2len["${myArray[0]}"]="${myArray[2]}"
  acc2tax["${myArray[0]}"]="${myArray[3]}"
done < ${outfolder}/${pfx}_acc2info.txt

##########################################
# add coverage and %coverage for each Acc
##########################################

while IFS=$'\t' read -r id 
do

  title=${acc2title["${id}"]:-"na"}
  chrom="$(echo ${title} | sed -r 's/.*chromosome[^0-9]*([0-9]+).*/\1/')"
  tlen=${acc2len["${id}"]:-"na"}
  aliw=${acc2aliw["${id}"]:-"na"}
  taxonomy=${acc2tax["${id}"]:-"na"}

  # Check if the chrom is an integer
  if [[ $chrom =~ ^-?[0-9]+$ ]]; then
    chr=${chrom}
  else
   chr="na"
  fi

  # compute % covered
  printf -v result "%.3f" "$(bc -l <<< "${aliw}*100/${tlen}")"
  if (( $(echo "$result < 1" | bc -l) )); then
    perc=$(printf "%.2f" "$result")
  else
    perc=$(printf "%.1f" "$result")
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${id}" "${title}" "chr_${chr}" "${tlen}" "${aliw}" "${perc}" "${taxonomy}"

done < ${outfolder}/${pfx}_acc_list.txt > ${outfolder}/${pfx}_acc_info.txt

#####################################
# parse blast counts and add columns
#####################################

# NOTE that "${aliw}" "${perc}" refer globally to the target accession 
# and not to that contig vs id alignment

while IFS=$'\t' read ctg id
do

  title=${acc2title["${id}"]:-"na"}
  chrom="$(echo ${title} | sed -r 's/.*chromosome[^0-9]*([0-9]+).*/\1/')"
  tlen=${acc2len["${id}"]:-"na"}
  aliw=${acc2aliw["${id}"]:-"na"}
  taxonomy=${acc2tax["${id}"]:-"na"}
  
  # Check if the chrom is an integer
  if [[ $chrom =~ ^-?[0-9]+$ ]]; then
      chr=${chrom}
  else
     chr="na"
  fi
  
  # compute % covered
  perc=$(echo "scale=1; ${aliw}*100/${tlen}" | bc)
  
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${ctg}" "${id}" "${title}" "chr_${chr}" "${tlen}" "${aliw}" "${perc}" "${taxonomy}"

done < ${outfolder}/${pfx}_blast_counts.txt \
  > ${outfolder}/${pfx}_blast_summary.txt

conda deactivate

####################################################
# extract tophit for genus and species from summary
####################################################

# sort decreasing and keep first as candidate
tophit=( $(cut -f 7 ${outfolder}/${pfx}_blast_summary.txt \
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
