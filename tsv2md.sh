#!/bin/bash
#
# tsv2md.sh - Convert TSV data to Markdown table
# Author: Stéphane Plaisance
# Version: 1.2 - 2025-07-17

usage() {
  echo "Usage: $0 [-h] [-H] [file]"
  echo "  -h    Show this help message"
  echo "  -H    Input has header line (default: no header)"
  echo "  file  Input TSV file (optional, else read from stdin)"
  exit 1
}

has_header=0

while getopts ":hH" opt; do
  case $opt in
    h) usage ;;
    H) has_header=1 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
  esac
done
shift $((OPTIND -1))

input="$1"
if [[ -n "$input" ]]; then
  exec < "$input"
fi

if [[ $has_header -eq 1 ]]; then
  # Print header from first line
  read -r header
  echo "$header" | awk -F'\t' '{
    printf "|"
    for(i=1;i<=NF;i++) printf " %s |", $i
    printf "\n|"
    for(i=1;i<=NF;i++) printf " --- |"
    printf "\n"
  }'
else
  # Read first line to count columns, print empty header
  read -r first
  ncol=$(awk -F'\t' '{print NF}' <<< "$first")
  printf "|"
  for ((i=1;i<=ncol;i++)); do printf "   |"; done
  printf "\n|"
  for ((i=1;i<=ncol;i++)); do printf " --- |"; done
  printf "\n"
  # Print the first line as data
  echo "$first" | awk -F'\t' '{
    printf "|"
    for(i=1;i<=NF;i++) printf " %s |", $i
    printf "\n"
  }'
fi

# Print remaining lines as data
awk -F'\t' '{
  printf "|"
  for(i=1;i<=NF;i++) printf " %s |", $i
  printf "\n"
}'