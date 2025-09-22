# Fungal rDNA Database Creation Pipeline

This project provides a reproducible pipeline for building a comprehensive fungal rDNA sequence database from NCBI genomes. The workflow is container-friendly and designed for use in a conda environment, following strict host/container path separation as described in `AI_INSTRUCTIONS.md`.

## Prerequisites
- [Conda](https://docs.conda.io/en/latest/)
- All commands are intended to run inside a container with `/data` as the workspace mount

## Setup
Create and activate the required conda environment:

```bash
conda create -n rDNA_database -c bioconda ncbi-datasets-cli barrnap itsx hmmer bedtools samtools cdhit blast kraken2
conda activate rDNA_database
```

## Pipeline Steps

1. **Download Fungal Genomes**
   - Downloads all fungal genomes from NCBI (taxon 4751).
   ```bash
datasets download genome taxon 4751 --reference --include genome,seq-report
```

2. **Decompress the Archive**
   ```bash
gunzip ncbi_dataset.zip
```

3. **Process Each Genome**
   For each genome FASTA file in `ncbi_dataset/data/*.fna`:
   - Predict rDNA genes with Barrnap
   - Extract rDNA sequences using BEDTools
   - Extract ITS regions with ITSx
   - Optionally, run HMMER against rRNA models
   - Extract sequences using coordinates from Barrnap or HMMER output

   Example loop:
   ```bash
for fasta in ncbi_dataset/data/*.fna; do
  base=$(basename "$fasta" .fna)
  barrnap --kingdom euk --threads 4 "$fasta" > "${base}_barrnap_output.gff"
  bedtools getfasta -fi "$fasta" -bed "${base}_barrnap_output.gff" -fo "${base}_extracted_rDNA.fa"
  ITSx -i "$fasta" -o "${base}_itsx_output" --preserve T --only_full T --cpu 4
  hmmsearch --tblout "${base}_rRNA_hits.tbl" fungal_rRNA_models.hmm "$fasta"
  # See commands.sh for details on extracting sequences from HMMER output
  # ...existing code...
done
```

4. **Compile All Extracted Sequences**
   ```bash
cat ncbi_dataset/data/*_extracted_rDNA.fa > fungal_rDNA_sequences.fa
```

5. **Dereplicate Sequences (Optional)**
   ```bash
cd-hit -i fungal_rDNA_sequences.fa -o fungal_rDNA_nr.fa -c 0.99 -n 5 -M 16000 -T 8
```

6. **Build BLAST Database**
   ```bash
makeblastdb -in fungal_rDNA_nr.fa -dbtype nucl -out fungal_rDNA_db
```

7. **Build Kraken2 Database (Optional)**
   ```bash
kraken2-build --add-to-library fungal_rDNA_nr.fa --db fungal_rDNA_kraken_db
kraken2-build --build --db fungal_rDNA_kraken_db
```

## Notes
- All file paths must be relative to the container image, not the host.
- Refer to `AI_INSTRUCTIONS.md` for architecture and flag conventions.
- Update `SESSION_DECISIONS.md` if you make changes to the workflow or flag definitions.

## References
- [Barrnap](https://github.com/tseemann/barrnap)
- [ITSx](https://github.com/ITSx/ITSx)
- [HMMER](http://hmmer.org/)
- [BEDTools](https://bedtools.readthedocs.io/en/latest/)
- [CD-HIT](https://github.com/weizhongli/cdhit)
- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- [Kraken2](https://ccb.jhu.edu/software/kraken2/)

---
For details and troubleshooting, see `commands.sh` and context files in this repository.
