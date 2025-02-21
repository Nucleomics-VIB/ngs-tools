#!/usr/bin/env python3
# script: nexus2phy.py
# Description: This script converts a NEXUS file produced by vcf2phylip.py
#              to a binary PHYLIP file.
#              It reads the NEXUS file, extracts the taxa and sequences,
#              converts heterozygous sites (2) to presence (1), and writes
#              the data to a PHYLIP file used as input by iQtree2.
# Author: SP@NC (+GitHub Copilot)
# Date: 2025-02-12
# Version: 1.1

import os
import sys

def nexus_to_binary_phylip(nexus_file):
    """
    Converts a NEXUS file to a binary PHYLIP file.

    Args:
        nexus_file (str): Path to the input NEXUS file.

    Returns:
        None
    """
    # Deduce output filename by replacing .nexus with .phy
    phylip_file = os.path.splitext(nexus_file)[0] + ".phy"

    with open(nexus_file, 'r') as nexus, open(phylip_file, 'w') as phylip:
        lines = nexus.readlines()
        taxa = []
        matrix_start = False

        # Parse NEXUS file
        for line in lines:
            if line.strip().lower().startswith("matrix"):
                matrix_start = True
                continue
            if matrix_start and ";" in line:
                break
            if matrix_start:
                parts = line.strip().split()
                if len(parts) < 2:  # Skip invalid lines
                    continue
                taxa.append(line.strip())

        # Write PHYLIP header: number of taxa and characters
        num_taxa = len(taxa)
        num_sites = len(taxa[0].split()[1])
        phylip.write(f"{num_taxa} {num_sites}\n")

        # Write taxa and sequences in binary format (convert 2 -> 1)
        for taxon in taxa:
            name, sequence = taxon.split()
            binary_sequence = sequence.replace("2", "1")  # Convert 2 to 1 (heterozygous -> presence)
            phylip.write(f"{name.ljust(10)} {binary_sequence}\n")

    print(f"Converted {nexus_file} to {phylip_file}")

# Main function to handle command-line arguments
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 reformat_nexus.py <input_nexus_file>")
        sys.exit(1)

    nexus_file = sys.argv[1]
    if not os.path.exists(nexus_file):
        print(f"Error: File {nexus_file} does not exist.")
        sys.exit(1)

    nexus_to_binary_phylip(nexus_file)


"""
# the following iqtree2 command is as follows
# create tree using binary format
# input=<your input from above>; bs=2000; al=100; nthr=84; mem=256G \
#  iqtree2 -s ${input} \
#    -st BIN \
#    -m TEST+ASC \
#    -bb ${bs} \
#    -alrt ${al} \
#    -nt ${nthr} \
#    -mem ${mem} \
#    -bnni \
#    > ${vcf%.vcf.gz}_iqtree2.log 2>&1
"""
