[(Nucleomics-VIB)](https://github.com/Nucleomics-VIB)
![ngs-tools](ngstools.png) - NGS-Tools
==========

*All tools presented below have only been tested by me and may contain bugs, please let me know if you find some. Each tool relies on dependencies normally listed at the top of the code (cpan for perl, R packages, etc.)*

Please refer to the accompanying **[wiki](https://github.com/Nucleomics-VIB/ngs-tools/wiki)** for examples and workflows.

## Fasta and Fastq Utilities

### **fastastats.pl**

The perl script **[fastastats.pl](fastastats.pl)** computes simple length statistics on a multi-fasta file.
```bash
## Usage: fastastats.pl <-i fasta-file (can be compressed)>
# Additional optional filtering parameters are:
# <-m minsize in kb (0)>
# <-x maxsize in kb (1E+6 | 1Gb)>
# <-h to display this help>
```

### **fastaCleanHeader.pl**

The perl script **[fastaCleanHeader.pl](fastaCleanHeader.pl)** cleans complex fasta record names that may interfere with some applications.
```bash
## Usage: fastaCleanHeader.pl <-i fasta_file (required)>
# <-o output file name (default to "cleaned_"<infile>; optional)>
# <-c keep only the leftmost word (display_id field; optional)>
# <-d delimiter (default to '|'; optional)>
# <-z to compress results with bgzip>
# <-h to display this help>
```

### **fasta2merge.pl**

The perl script **[fasta2merge.pl](fasta2merge.pl)** creates a single fasta sequence from a multifasta file (required for ABACAS -r ref.fasta).
```bash
## Usage: fasta2merge.pl <-i fasta-file>
# script version 1.0, 2017_09_05
# Additional optional parameters are:
# <-o prefix for output (default to merged-<inputname>.fasta)>
# <-h to display this help>
```

### **fastabreak.pl**

The perl script **[fastabreak.pl](fastabreak.pl)** breaks a multifasta file into individual fasta files.
```bash
## Usage: fastabreak.pl <-i fasta-file>
# script version 1.0, 2017_09_05
# Additional optional parameters are:
# <-o output directory (default to current directory)>
# <-h to display this help>
```

### **fastarotate.pl**

The perl script **[fastarotate.pl](fastarotate.pl)** rotates sequences in a fasta file.
```bash
## Usage: fastarotate.pl <-i fasta-file>
# script version 1.0, 2017_09_05
# Additional optional parameters are:
# <-o output file (default to rotated_<inputname>.fasta)>
# <-h to display this help>
```

### **fastaSplit.pl**

The perl script **[fastaSplit.pl](fastaSplit.pl)** splits a multifasta file into smaller files.
```bash
## Usage: fastaSplit.pl <-i fasta-file>
# script version 1.0, 2017_09_05
# Additional optional parameters are:
# <-o output prefix (default to split_<inputname>)>
# <-n number of sequences per file (default to 1)>
# <-h to display this help>
```

### **fasta2consensus.sh**

The bash script **[fasta2consensus.sh](fasta2consensus.sh)** generates a consensus sequence from a multifasta file.
```bash
# Usage: fasta2consensus.sh -i <fasta file>
# script version 1.0, 2025_03_14
# [optional: -o <output file>]
# [optional: -h <this help text>]
```

### **fastp_benchmark.sh**

The bash script **[fastp_benchmark.sh](fastp_benchmark.sh)** benchmarks fastp for fastq file processing.
```bash
# Usage: fastp_benchmark.sh -i <fastq file>
# script version 1.0, 2025_03_14
# [optional: -o <output directory>]
# [optional: -t <number of threads>]
# [optional: -h <this help text>]
```

## Assembly and Alignment

### **run_Hifiasm_meta.sh**

The bash script **[run_Hifiasm_meta.sh](run_Hifiasm_meta.sh)** processes PacBio HiFi data through conversion (if needed) and assembly.
```bash
# Usage: run_Hifiasm_meta.sh [-q] [-b INPUT_BAM_DIR] [-f FASTQ_DIR] [-o OUTPUT_ASM_DIR] [-n THREADS] [-j MAX_JOBS] [-t ASM_THREADS]
# script version 1.0, 2025_03_14
# Options:
#   -q          start directly from FASTQ files (default: process BAM files)
#   -b DIR      input BAM directory (default: bam_data)
#   -f DIR      FASTQ input/output directory (default: fastq_data)
#   -o DIR      assembly output directory (default: asm_results)
#   -n INT      threads for bam2fastq (default: 4)
#   -j INT      max concurrent jobs (default: 4)
#   -t INT      threads per hifiasm_meta (default: 20)
```

### **run_freebayes.sh**

The bash script **[run_freebayes.sh](run_freebayes.sh)** calls variants with freebayes from mappings and a reference genome.
```bash
# Usage: run_freebayes.sh -i <bam file> -r <fasta reference>
# script version 1.0, 2016_09_28
# [optional: -m <minmapq|20>]
# [optional: -q <minbaseq|20>]
# [optional: -F <minaltfrac|0.01>]
# [optional: -C <minaltcnt|10>]
```

### **makeENSref.sh**

The bash script **[makeENSref.sh](makeENSref.sh)** creates a series of reference files from ENSembl FTP downloads.
```bash
# Usage: makeENSref.sh
# -o <organism (default to <homo_sapiens>)> 
# -b <build number (default to <GRCh38>)> 
# -r <release number (default to 88)>
# script version 1.0, 2017_04_05
# [-h for this help]
```

### **gepard_plot.sh**

The bash script **[gepard_plot.sh](gepard_plot.sh)** creates a xy-plot from two related assemblies. The Java GUI tools do the same but this script is applicable to multiple inputs in batch.
```bash
# Usage: gepard_plot.sh -x <reference assembly> -y <draft assembly> -p <path to gepard.jar and matrices>
# script version 1.1, 2017_09_05
# [optional: -o <result path/prefix>]
# [optional: -w <word size:10>]
# [optional: -W <window size:0>]
# [optional: -l <lower value% (default:0)>]
# [optional: -u <upper value% (default:100)>]
# [optional: -g <greyscale start value% (default:0)>]
# [optional: -J <java extra parameters (eg -Xmx1G, put between double quotes if it contains spaces)>
# [optional: -h <this help text>]
```

### **mauve_reorder.sh**

The bash script **[mauve_reorder.sh](mauve_reorder.sh)** reorders a draft-assembly based on a reference assembly at CLI. The Java GUI tools do the same but this script is applicable to multiple inputs in batch.
```bash
# Usage: mauve_reorder.sh -i <draft assembly> -r <reference assembly> -p <mauve path>
# script version 1.0, 2017_04_21
# [optional: -o <result folder>]
# [optional: -m <available RAM|1G>]
# [optional: -h <this help text>]
```

### **mappability.sh**

The bash script **[mappability.sh](mappability.sh)** creates a mappability track for a given single-read length and a reference genome. Such tracks used to be present in IGV from the web server and are now created locally.
```bash
# Usage: mappability.sh
# -i <reference assembly fasta>)> 
# -l <single read length for prediction (default to 100bps)>
# -p <prefix for output data (default to input file prefix)>
# -t <number of threads (default to 1)> 
# script version 1.0, 2017_05_19
# [-h for this help]
```

## Variant Calling and Analysis

### **GATK_Variant_Analysis.sh**

The bash script **[GATK_Variant_Analysis.sh](GATK_Variant_Analysis.sh)** performs variant analysis using GATK.
```bash
# Usage: GATK_Variant_Analysis.sh -i <input file> -r <reference genome>
# script version 1.0, 2025_03_14
# [optional: -o <output directory>]
# [optional: -t <number of threads>]
# [optional: -h <this help text>]
```

### **run_Deepvariant.sh**

The bash script **[run_Deepvariant.sh](run_Deepvariant.sh)** runs DeepVariant for variant calling.
```bash
# Usage: run_Deepvariant.sh -i <input file> -r <reference genome>
# script version 1.0, 2025_03_14
# [optional: -o <output directory>]
# [optional: -t <number of threads>]
# [optional: -h <this help text>]
```

### **varscan2mpileup2cns.sh**

The bash script **[varscan2mpileup2cns.sh](varscan2mpileup2cns.sh)** calls variants using Varscan2 and mpileup.
```bash
# Usage: varscan2mpileup2cns.sh -i <input file> -r <reference genome>
# script version 1.0, 2025_03_14
# [optional: -o <output directory>]
# [optional: -t <number of threads>]
# [optional: -h <this help text>]
```

### **varscan2mpileup2somatic.sh**

The bash script **[varscan2mpileup2somatic.sh](varscan2mpileup2somatic.sh)** calls somatic variants using Varscan2 and mpileup.
```bash
# Usage: varscan2mpileup2somatic.sh -i <input file> -r <reference genome>
# script version 1.0, 2025_03_14
# [optional: -o <output directory>]
# [optional: -t <number of threads>]
# [optional: -h <this help text>]
```

### **varscan2mpileup2trio.sh**

The bash script **[varscan2mpileup2trio.sh](varscan2mpileup2trio.sh)** calls trio variants using Varscan2 and mpileup.
```bash
# Usage: varscan2mpileup2trio.sh -i <input file> -r <reference genome>
# script version 1.0, 2025_03_14
# [optional: -o <output directory>]
# [optional: -t <number of threads>]
# [optional: -h <this help text>]
```

## Other Tools

### **asm2cvg.sh**

The bash script **[asm2cvg.sh](asm2cvg.sh)** computes coverage from an assembly file.
```bash
# Usage: asm2cvg.sh -i <assembly file>
# script version 1.0, 2025_03_14
# [optional: -o <output file>]
# [optional: -h <this help text>]
```

### **asm2dotplot.sh**

The bash script **[asm2dotplot.sh](asm2dotplot.sh)** creates a dot plot from an assembly file.
```bash
# Usage: asm2dotplot.sh -i <assembly file>
# script version 1.0, 2025_03_14
# [optional: -o <output file>]
# [optional: -h <this help text>]
```

### **bcftoolsParallel.sh**

The bash script **[bcftoolsParallel.sh](bcftoolsParallel.sh)** runs bcftools in parallel.
```bash
# Usage: bcftoolsParallel.sh -i <input file>
# script version 1.0, 2025_03_14
# [optional: -o <output directory>]
# [optional: -t <number of threads>]
# [optional: -h <this help text>]
```

### **bed2density.sh**

The bash script **[bed2density.sh](bed2density.sh)** converts a BED file to a density file.
```bash
# Usage: bed2density.sh -i <BED file>
# script version 1.0, 2025_03_14
# [optional: -o <output file>]
# [optional: -h <this help text>]
```

### **blast_assembly2NT.sh**

The bash script **[blast_assembly2NT.sh](blast_assembly2NT.sh)** performs BLAST on an assembly against the NT database.
```bash
# Usage: blast_assembly2NT.sh -i <assembly file>
# script version 1.0, 2025_03_14
# [optional: -o <output directory>]
# [optional: -t <number of threads>]
# [optional: -h <this help text>]
```

### **circ-perm.pl**

The perl script **[circ-perm.pl](circ-perm.pl)** performs circular permutation on sequences.
```bash
## Usage: circ-perm.pl <-i input file>
# script version 1.0, 2025_03_14
# Additional optional parameters are:
# <-o output file (default to circ_perm_<inputname>.txt)>
# <-h to display this help>
```

### **cmpasm2vcf.sh**

The bash script **[cmpasm2vcf.sh](cmpasm2vcf.sh)** compares assemblies and generates a VCF file.
```bash
# Usage: cmpasm2vcf.sh -i <input file>
# script version 1.0, 2025_03_14
# [optional: -o <output file>]
# [optional: -h <this help text>]
```

### **compare_all2all.sh**

The bash script **[compare_all2all.sh](compare_all2all.sh)** compares all sequences to all other sequences.
```bash
# Usage: compare_all2all.sh -i <input file>
# script version 1.0, 2025_03_14
# [optional: -o <output file>]
# [optional: -h <this help text>]
```

### **consensusQualFilter.sh**

The bash script **[consensusQualFilter.sh](consensusQualFilter.sh)** filters consensus sequences based on quality.
```bash
# Usage: consensusQualFilter.sh -i <input file>
# script version 1.0, 2025_03_14
# [optional: -o <output directory>]
# [optional: -t <number of threads>]
# [optional: -h <this help text>]
```

### **create_xenome-idx_all.sh**

The bash script **[create_xenome-idx_all.sh](create_xenome-idx_all.sh)** creates Xenome indices for all input files.
```bash
# Usage: create_xenome-idx_all.sh -i <input file>
# script version 1.0, 2025_03_14
# [optional: -o <output directory>]
# [optional: -t <number of threads>]
# [optional: -h <this help text>]
```

### **DownsampleBam_x.sh**

The bash script **[DownsampleBam_x.sh](DownsampleBam_x.sh)** downsamples BAM files.
```bash
# Usage: DownsampleBam_x.sh -i <input file>
# script version 1.0, 2025_03_14
# [optional: -o <output directory>]
# [optional: -t <number of threads>]
# [optional: -h <this help text>]
```

### **fastq2metrics.awk**

The awk script **[fastq2metrics.awk](fastq2metrics.awk)** calculates metrics from a fastq file.
```bash
# Usage: fastq2metrics.awk -i <fastq file>
# script version 1.0, 2025_03_14
# [optional: -o <output file>]
# [optional: -h <this help text>]
```

### **fastq2metrics.cpp**

The C++ script **[fastq2metrics.cpp](fastq2metrics.cpp)** calculates metrics from a fastq file.
```bash
# Usage: fastq2metrics.cpp -i <fastq file>
# script version 1.0, 2025_03_14
# [optional: -o <output file>]
# [optional: -h <this help text>]
```

### **fastq2metrics60.cpp**

The C++ script **[fastq2metrics60.cpp](fastq2metrics60.cpp)** calculates metrics from a fastq file.
```bash
# Usage: fastq2metrics60.cpp -i <fastq file>
# script version 1.0, 2025_03_14
# [optional: -o <output file>]
# [optional: -h <this help text>]
```

### **fastq2motif.py**

The python script **[fastq2motif.py](fastq2motif.py)** identifies motifs in fastq sequences.
```bash
# Usage: fastq2motif.py -i <fastq file>
# script version 1.0, 2025_03_14
# [optional: -o <output file>]
# [optional: -h <this help text>]
```

### **
