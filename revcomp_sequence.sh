#!/bin/bash

# revcomp_sequence.sh - DNA sequence manipulation tool
# 
# Description:
#   This script takes a DNA sequence (potentially containing IUPAC ambiguity codes)
#   and reports the original sequence, reverse sequence, and reverse complement.
#   All IUPAC nucleotide ambiguity codes are properly handled.
#
# IUPAC nucleotide codes:
#   A = Adenine       M = A or C (aMino)
#   C = Cytosine      S = G or C (Strong)
#   G = Guanine       W = A or T (Weak)
#   T = Thymine       K = G or T (Keto)
#   U = Uracil (NA)   R = A or G (puRine)
#   Y = C or T (pYrimidine)
#   B = C, G, or T (not A)
#   D = A, G, or T (not C)
#   H = A, C, or T (not G)
#   V = A, C, or G (not T)
#   N = any base (A, C, G, or T)
#
# Usage: ./revcomp_sequence.sh "ACGTNRYSWKMBDHV"
#
# Author: SP@NC (email nucleomics@vib.be)
# Version: 1.0
# inspired by https://github.com/vmikk/NextITS/blob/main/bin/convert_IUPAC.sh (by @vmikk) 
# Date: June 17, 2025

# Define IUPAC complement mapping
declare -A complement_map
complement_map[A]=T
complement_map[T]=A
complement_map[G]=C
complement_map[C]=G
complement_map[U]=A
complement_map[R]=Y  # A/G -> T/C
complement_map[Y]=R  # C/T -> G/A
complement_map[S]=S  # G/C -> C/G (self-complementary)
complement_map[W]=W  # A/T -> T/A (self-complementary)
complement_map[K]=M  # G/T -> C/A
complement_map[M]=K  # A/C -> T/G
complement_map[B]=V  # C/G/T -> G/C/A
complement_map[D]=H  # A/G/T -> T/C/A
complement_map[H]=D  # A/C/T -> T/G/A
complement_map[V]=B  # A/C/G -> T/G/C
complement_map[N]=N  # any -> any
complement_map["-"]="-" # gap remains gap

# Function to reverse a string
reverse_sequence() {
    echo "$1" | rev
}

# Function to generate the complement of a sequence (supports IUPAC codes)
complement_sequence() {
    local seq="$1"
    local result=""
    
    # Process each character using the complement map
    for (( i=0; i<${#seq}; i++ )); do
        char="${seq:$i:1}"
        comp="${complement_map[${char^^}]}"
        if [ -n "$comp" ]; then
            result+="$comp"
        else
            # If unknown character, keep it unchanged
            result+="$char"
        fi
    done
    
    echo "$result"
}

# Function to generate reverse complement
reverse_complement() {
    local seq="$1"
    local complement=$(complement_sequence "$seq")
    local rev_complement=$(reverse_sequence "$complement")
    echo "$rev_complement"
}

# Function to clean the input sequence
clean_sequence() {
    local seq="$1"
    # First trim leading/trailing whitespace
    seq=$(echo "$seq" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
    # Convert to uppercase
    seq=$(echo "$seq" | tr '[:lower:]' '[:upper:]')
    # Replace U with T (RNA to DNA conversion)
    seq=$(echo "$seq" | tr 'U' 'T')
    # Replace internal spaces with hyphens
    seq=$(echo "$seq" | tr ' ' '-')
    # Remove other control characters
    seq=$(echo "$seq" | tr -d '\r\n\t\f\v')
    echo "$seq"
}

# Main script
if [ $# -eq 0 ]; then
    echo "Usage: $0 <DNA_sequence>"
    echo "Example: $0 CGACCWGCGGARGGATCATTA"
    exit 1
fi

# Clean and normalize the input sequence
raw_sequence="$1"
sequence=$(clean_sequence "$raw_sequence")
sequence="${sequence^^}"  # Convert to uppercase

if [ -z "$sequence" ]; then
    echo "Error: Empty sequence after cleaning. Please provide a valid DNA sequence."
    exit 1
fi

echo "Original sequence  5'→3': $sequence"
echo "Reverse sequence   3'→5': $(reverse_sequence "$sequence")"
echo "Reverse complement 5'→3': $(reverse_complement "$sequence")"

# Optional: Show whether the sequence contains IUPAC ambiguity codes
if echo "$sequence" | grep -q -E "[RYSWKMBDHVN]"; then
    echo -e "\nNote: This sequence contains IUPAC ambiguity codes."
fi

# Optional: Show if input was cleaned
if [[ "$raw_sequence" != "$sequence" && "${raw_sequence^^}" != "$sequence" ]]; then
    echo "Note: Input sequence was cleaned (spaces replaced with hyphens, control characters removed, U=>T)."
fi

exit 0