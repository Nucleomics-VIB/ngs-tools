#!/bin/bash

# create conda env with required tools
conda create -n rDNA_database -c conda-forge -c bioconda ncbi-datasets-cli barrnap itsx hmmer bedtools samtools cd-hit blast kraken2 infernal
conda activate rDNA_database

# Download SSU and LSU rRNA HMMs from Rfam (run once)
export BIODATA=/opt/BIODATA
mkdir -p $BIODATA/rfam_models

#1 Download the full covariance model database from Rfam's FTP site
This will give you the concatenated file Rfam.cm containing all family models.
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz

#2. Extract Individual Models Using infernal's cmfetch

# Extract RF01960 (SSU) RF00002 (5.8S) and RF02543 (LSU) to separate files
cmfetch Rfam.cm RF01960 > RF01960.hmm
cmfetch Rfam.cm RF02543 > RF02543.hmm
cmfetch Rfam.cm RF00002 > RF00002.hmm

# get all fungi genomes (5232, takes time!)
# https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/large-download/
export NCBI_API_KEY=<mykey>
datasets download genome taxon 4751 --reference --include genome,seq-report --excluide mag --dehydrate

# decompress the archive
unzip ncbi_dataset.zip

# rehydrate
datasets rehydrate --directory ncbi_dataset

# extract classification data from NCBI downlod
# datasets summary genome accession GCF_038497785.1 | jq '.reports[0].organism | {organism_name, tax_id}'
# datasets summary genome accession GCF_038497785.1 | jq -r '.reports[0].organism.organism_name'
# acc=GCF_038497785.1
# datasets summary genome accession $acc | jq -r "\"$acc\t\" + .reports[0].organism.organism_name"

####################################################
# for each genome fasta (.fna) in ncbi_dataset/data

# Predict rDNA genes with Barrnap
outfolder=barrnap_results
mkdir -p ${outfolder}

for fasta in $(find . -name "*.fna"); do
    base=$(basename $(dirname "$fasta"))

    # 1. Extract classification and append to classification.txt
	datasets summary genome accession "$base" | jq -r "\"$base\t\" + .reports[0].organism.organism_name" >> classification.txt

	# 2. Predict rDNA genes with Barrnap
	barrnap --kingdom euk --threads 4 "$fasta" > "${outfolder}/${base}_barrnap_output.gff"

	# 3. Extract only contiguous SSU-LSU regions using Python script output
	python scripts/find_rDNA_region.py "${outfolder}/${base}_barrnap_output.gff" > "${outfolder}/${base}_contiguous_rDNA.bed"
    
	# Only extract fasta if BED file is not empty
	if [ -s "${outfolder}/${base}_contiguous_rDNA.bed" ]; then
		bedtools getfasta -fi "$fasta" -bed "${outfolder}/${base}_contiguous_rDNA.bed" -fo "${outfolder}/${base}_extracted_rDNA.fa"
	fi
done	

# Compile all extracted sequences from multiple genomes into one fasta
cat ncbi_dataset/data/*_extracted_rDNA.fa > fungal_rDNA_sequences.fa

# Optionally dereplicate
cd-hit -i fungal_rDNA_sequences.fa -o fungal_rDNA_nr.fa -c 0.99 -n 5 -M 16000 -T 8

# Build BLAST database
makeblastdb -in fungal_rDNA_nr.fa -dbtype nucl -out fungal_rDNA_db

# Or build Kraken2 database if taxonomic annotations are available
kraken2-build --add-to-library fungal_rDNA_nr.fa --db fungal_rDNA_kraken_db
kraken2-build --build --db fungal_rDNA_kraken_db



exit 0

# 3. Extract ITS regions with ITSx
# will not work as downloaded sequences are too long for this tool which needs amplicon to search against
#outfolder=itsx_results
#mkdir -p ${outfolder}
#ITSx -i "$fasta" -o "${base}_itsx_output" --preserve T --only_full T --cpu 4
# 4. Alternatively, run HMMER against rRNA models (downloaded from Rfam)
#outfolder=hmmsearch_results
#mkdir -p ${outfolder}
# Run HMMER for SSU and LSU separately
#hmmsearch --tblout "${outfolder}/${base}_ssu_hits.tbl" rfam_models/RF01960.hmm "$fasta"
#hmmsearch --tblout "${outfolder}/${base}_5.8S_hits.tbl" rfam_models/RF00002.hmm "$fasta"
#hmmsearch --tblout "${outfolder}/${base}_lsu_hits.tbl" rfam_models/RF02543.hmm "$fasta"
# Combine and parse HMMER results to find contiguous SSU-LSU regions
#python scripts/find_rDNA_region_hmmer.py "${outfolder}/${base}_ssu_hits.tbl" "${outfolder}/${base}_lsu_hits.tbl" > "${outfolder}/${base}_contiguous_rDNA_hmmer.bed"
# Only extract fasta if BED file is not empty
#if [ -s "${outfolder}/${base}_contiguous_rDNA_hmmer.bed" ]; then
#	bedtools getfasta -fi "$fasta" -bed "${outfolder}/${base}_contiguous_rDNA_hmmer.bed" -fo "${outfolder}/${base}_extracted_rRNA_hmmer.fa"
#fi
# 5. Extract sequences using coordinates from hmmsearch or barrnap output
# Example for Barrnap output (already done above):
# bedtools getfasta -fi "$fasta" -bed "${base}_barrnap_output.gff" -fo "${base}_extracted_rDNA.fa"
# Example for HMMER output (requires conversion of .tbl to BED format):
# (Assuming you have a script or command to convert "${base}_rRNA_hits.tbl" to BED format)
# bedtools getfasta -fi "$fasta" -bed "${base}_rRNA_hits.bed" -fo "${base}_extracted_rRNA_hmmer.fa"
