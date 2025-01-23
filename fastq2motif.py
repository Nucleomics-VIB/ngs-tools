#!/usr/bin/env python3

"""
fastq2motif.py

Created on: January 23, 2025
Version: 2.7
Author: SP@NC + AI Assistant

Description:
This script computes the DNA motif from the first <l (default 20)> bases of <s (default 1000)> reads in a FASTQ or FASTQ.gz file.
It processes the specified number of reads, outputs the consensus sequence and count matrix, generates a sequence logo,
and saves the motif matrix in MEME format.

Usage:
python fastq2motif.py -i <input_fastq> [-l <length (20)>] [-s <sample_size (1000)>]
"""

import sys
import gzip
import argparse
from Bio import SeqIO
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
import numpy as np
import os

def compute_motif(fastq_file, motif_length, sample_size):
    sequences = []
    
    open_func = gzip.open if fastq_file.endswith('.gz') else open
    mode = 'rt' if fastq_file.endswith('.gz') else 'r'
    
    with open_func(fastq_file, mode) as handle:
        for record in SeqIO.parse(handle, 'fastq'):
            if len(sequences) >= sample_size:
                break
            sequences.append(str(record.seq[:motif_length]))
    
    motif_counts = [Counter(seq[i] for seq in sequences) for i in range(motif_length)]
    consensus = ''.join(max(pos_counts, key=pos_counts.get) for pos_counts in motif_counts)
    
    return consensus, motif_counts

def create_logo(motif_counts, file_prefix):
    # Create a DataFrame with the correct structure for logomaker
    df = pd.DataFrame({
        base: [counts.get(base, 0) for counts in motif_counts]
        for base in 'ACGT'
    })
    
    # Rename the index to start from 1
    df.index = range(1, len(df) + 1)
    
    # Normalize counts to probabilities
    df = df.div(df.sum(axis=1), axis=0)
    
    # Replace any NaN or inf values with 0
    df = df.fillna(0).replace([np.inf, -np.inf], 0)
    
    # Create logo
    fig, ax = plt.subplots(figsize=(10, 3))
    logo = logomaker.Logo(df, ax=ax)
    
    # Style the logo
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    logo.style_xticks(rotation=90, fmt='%d', anchor=0)
    
    # Add labels and title
    ax.set_xlabel('Position')
    ax.set_ylabel('Probability')
    ax.set_title(f"DNA Motif Logo - {file_prefix}")
    
    # Save the logo
    plt.savefig("motif_logo.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Save the matrix in MEME format
    with open("motif_matrix.txt", "w") as f:
        f.write("MEME version 4\n\n")
        f.write("ALPHABET= ACGT\n")
        f.write("strands: + -\n\n")
        f.write("MOTIF motif_1\n")
        f.write(f"letter-probability matrix: alength= 4 w= {len(motif_counts)} nsites= {sum(motif_counts[0].values())}\n")
        for counts in motif_counts:
            total = sum(counts.values())
            probs = [counts.get(base, 0) / total for base in 'ACGT']
            f.write(" ".join(f"{prob:.6f}" for prob in probs) + "\n")

    print("\nMotif matrix saved as 'motif_matrix.txt' in MEME format")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute DNA motif from FASTQ file.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ or FASTQ.gz file")
    parser.add_argument("-l", "--length", type=int, default=20, help="Length of motif to analyze (default: 20)")
    parser.add_argument("-s", "--sample_size", type=int, default=1000, help="Number of reads to sample (default: 1000)")
    
    args = parser.parse_args()

    # Extract file prefix
    file_prefix = os.path.splitext(os.path.basename(args.input))[0]
    if file_prefix.endswith('.fastq') or file_prefix.endswith('.fq'):
        file_prefix = os.path.splitext(file_prefix)[0]

    print(f"Analyzing first {args.length} bases of {args.sample_size} reads")
    consensus, motif_counts = compute_motif(args.input, args.length, args.sample_size)

    print("\nConsensus sequence:")
    print(consensus)
    print("\nCount matrix:")
    for i, counts in enumerate(motif_counts, 1):
        print(f"Position {i}: {dict(counts)}")

    create_logo(motif_counts, file_prefix)
    print("\nMotif logo saved as 'motif_logo.png'")
