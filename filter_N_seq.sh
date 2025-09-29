#!/bin/bash

# filter_N_seq.sh
# Description: Split a fasta file into two files: 
# - one with sequences containing only N/./- (empty)
# - ne with all other sequences (noempty).
# Version: 1.0
# Date: 2025-09-29
# Author: SP@NC (Copilot-assisted)

# Usage: filter_N_seq.sh -i input.fasta
while getopts "i:" opt; do
  case $opt in
    i) infile="$OPTARG" ;;
    *) echo "Usage: $0 -i input.fasta"; exit 1 ;;
  esac
done

if [[ -z "$infile" ]]; then
  echo "Usage: $0 -i input.fasta"
  exit 1
fi

# Check input file exists and is readable
if [[ ! -r "$infile" ]]; then
  echo "ERROR: Input file '$infile' does not exist or is not readable."
  exit 2
fi

# Check output folder is writable
outdir=$(dirname "$infile")
if [[ ! -w "$outdir" ]]; then
  echo "ERROR: Output directory '$outdir' is not writable."
  exit 3
fi

prefix="${infile%.*}"
empty_out="${prefix}.empty.fasta"
noempty_out="${prefix}.noempty.fasta"

awk -v empty="$empty_out" -v noempty="$noempty_out" '
/^>/{
  if(seq) {
    if(seq ~ /^[-N\.]*$/) {
      print header > empty; print seq > empty
    } else {
      print header > noempty; print seq > noempty
    }
  }
  header=$0; seq=""
}
/^[^>]/ {seq=seq $0}
END {
  if(seq) {
    if(seq ~ /^[-N\.]*$/) {
      print header > empty; print seq > empty
    } else {
      print header > noempty; print seq > noempty
    }
  }
}' "$infile"