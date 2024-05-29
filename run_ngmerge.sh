#!/bin/env bash

# author:Stephane Plaisance (VIB-NC), 2022-05-17
# parse a folder with _R1 and _R2 matched read files
# merge paired reads with NGmerge
# https://github.com/jsh58/NGmerge
# installed from /releases/tag/v0.3
# run with params or defaults
# version 1.1; 2024-05-21

# default values
ol=20
thr=1
readfolder="reads"
pmm=0.10

# help message for usage
usage() {
  echo "Usage: $0 [-o overlap] [-t threads] [-r readfolder] [-p percent_mismatch]" 1>&2
}

# getopts to parse command-line arguments
while getopts ":o:t:r:p:" opt; do
  case ${opt} in
    o )
      ol=${OPTARG}
      ;;
    t )
      thr=${OPTARG}
      ;;
    r )
      readfolder=${OPTARG}
      # validate read folder
      if [ ! -d "${readfolder}" ] || [ -z "$(ls -A ${readfolder}/*R1*)" ]; then
        echo "Error: Read folder does not exist or contains no files with R1 in the name."
        exit 1
      fi
      ;;
    p )
      pmm=${OPTARG}
      ;;
    \? )
      echo "Invalid option: -$OPTARG" 1>&2
      usage
      exit 1
      ;;
    : )
      echo "Invalid option: -$OPTARG requires an argument" 1>&2
      usage
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))

myenv=ngmerge
source /etc/profile.d/conda.sh
conda activate ${myenv} || \
  ( echo "# the conda environment ${myenv} was not found on this machine" \
    && exit 1 )

outfolder="merged_reads"
mkdir -p ${outfolder}

for s in $(find "${readfolder}" -name "*R1*"); do
  pfx=$(basename ${s%_R1*})

  echo "# merging reads for sample ${pfx}"

  NGmerge \
    -1 ${s} \
    -2 ${s/_R1/_R2} \
    -o ${outfolder}/${pfx}_merged_ol${ol} \
    -f ${outfolder}/${pfx}_failed_ol${ol} \
    -m ${ol} \
    -p ${pmm} \
    -n ${thr} \
    -z \
    -j ${outfolder}/${pfx}_merged_ol${ol}_ali.txt \
    -l ${outfolder}/${pfx}_merged_ol${ol}_log.txt \
    -v

done

conda deactivate

exit 0
