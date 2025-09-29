#!/bin/bash

# buildContrifugeDB.sh
# create centrifuge database from local downloads
# SP@NC; 2025-09-25 (+copilot)

env=centrifuge
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate "${myenv}" || {
  echo "# the conda environment ${myenv} was not found on this machine"
  echo "# please read the top part of the script!"
  exit 1
}

fasta_paths=""
thr=4

while getopts "f:t:" opt; do
  case ${opt} in
    f) fasta_paths="$OPTARG" ;;
    t) thr="$OPTARG" ;;
    *) echo "Usage: $0 [-f fasta_path1,fasta_path2] [-t threads]"
       exit 1 ;;
  esac
done
shift $((OPTIND -1))

if [[ -z "$fasta_paths" ]]; then
  echo "ERROR: -f option with comma-separated fasta file paths is required"
  echo "Usage: $0 -f fasta_path1,fasta_path2"
  exit 1
fi

# Validate threads as positive integer if needed
if ! [[ "$thr" =~ ^[0-9]+$ ]]; then
  echo "ERROR: Threads (-t) must be a positive integer"
  exit 1
fi

# remove trailing spaces and compute pfx
fasta_paths="${fasta_paths//[[:space:]]/}"
IFS=',' read -ra files <<< "$fasta_paths"
pfx=""


# Build prefix (sorted) and collect map files
bases=()
map_files=()
for file in "${files[@]}"; do
  base=$(basename "$file")
  base=${base%.fa}
  base=${base%.fasta}
  base=${base%.fna}
  bases+=("$base")
  dir=$(dirname "$file")
  map_file="${dir}/seqid2taxid_${base}.map"
  if [[ ! -f "$map_file" ]]; then
    echo "ERROR: Map file not found: $map_file"
    exit 1
  fi
  map_files+=("$map_file")
done

# Sort bases alphabetically and join with _
IFS=$'\n' sorted=($(sort <<<"${bases[*]}"))
unset IFS
pfx="${sorted[*]}"
pfx="${pfx// /_}"

outfolder="${pfx}_DB"
mkdir -p "$outfolder"
merged_map="${outfolder}/seqid2taxid_${pfx}.map"

# Concatenate all map files into the merged map file
cat /dev/null > "$merged_map"
for mf in "${map_files[@]}"; do
  cat "$mf" >> "$merged_map"
done

fasta="$fasta_paths"
dbidx="${outfolder}/${pfx}"
logfile="${dbidx}_buildlog.txt"

cmd="time centrifuge-build -p ${thr} \
  --conversion-table ${merged_map} \
  --taxonomy-tree ncbi_taxonomy/nodes.dmp \
  --name-table ncbi_taxonomy/names.dmp \
  ${fasta} \
  ${dbidx} 2>&1 | tee ${logfile}"

echo "# ${cmd}"
eval ${cmd}
