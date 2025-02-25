#!/usr/bin/env Rscript

## VCF_to_PCA.R
#
# This script creates a PCA plot from all variants in a VCF file.
# It converts the VCF to GDS format, performs PCA, and generates a plot.
#
# Usage:
#   Rscript VCF_to_PCA.R -i <input_vcf_file> -o <output_directory> [-f <output_format (png+, pdf)>]
#
# Example:
#   Rscript VCF_to_PCA.R -i input.vcf.gz -o /path/to/output/
#   Rscript VCF_to_PCA.R -i input.vcf.gz -o /path/to/output/ -f pdf
#
# Author: SP@NC (+AI)
# Date: 2025-02-25
# Version: 1.0

# Function to check and load required libraries
check_and_load_libraries <- function() {
  required_packages <- c("SNPRelate", "ggplot2", "ggrepel", "optparse")
  
  for (package in required_packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(paste("Error: Required package", package, "is not installed. Please install it using install.packages('", package, "')"), call. = FALSE)
    }
  }
  
  # Load the libraries
  library(SNPRelate)
  library(ggplot2)
  library(ggrepel)
  library(optparse)
}

# Parse command line arguments
parse_arguments <- function() {
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="Input VCF file", metavar="FILE"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output directory", metavar="DIRECTORY"),
    make_option(c("-f", "--format"), type="character", default="png",
                help="Output format (png or pdf) [default= %default]", metavar="FORMAT")
  )
  
  opt_parser <- OptionParser(option_list=option_list,
                             description="Generate PCA plot from VCF file")
  opt <- parse_args(opt_parser)
  
  # Validate input
  if (is.null(opt$input) || is.null(opt$output)) {
    print_help(opt_parser)
    stop("Input file and output directory must be specified", call.=FALSE)
  }
  
  if (!file.exists(opt$input)) {
    stop(paste("Error: Input file does not exist:", opt$input))
  }
  
  if (!dir.exists(opt$output)) {
    stop(paste("Error: Output directory does not exist:", opt$output))
  }
  
  if (!(opt$format %in% c("png", "pdf"))) {
    stop("Error: Invalid output format. Use 'png' or 'pdf'.")
  }
  
  return(opt)
}

# Main function
main <- function(input_file, output_dir, output_format) {
  # Set working directory
  setwd(output_dir)
  
  # Convert VCF to GDS format
  gds_file <- file.path(output_dir, "converted.gds")
  snpgdsVCF2GDS(input_file, gds_file, method="biallelic.only")
  
  # Perform PCA
  genofile <- snpgdsOpen(gds_file)
  pca <- snpgdsPCA(genofile, autosome.only=FALSE)
  
  # Extract results
  pc.percent <- pca$varprop * 100
  sample.id <- pca$sample.id
  pc.matrix <- data.frame(sample.id, pca$eigenvect[, 1:2])
  colnames(pc.matrix) <- c("sample.id", "PC1", "PC2")
  
  # Create PCA plot
  plot <- ggplot(pc.matrix, aes(x = PC1, y = PC2, label = sample.id)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, max.overlaps = Inf) +
    xlab(paste0("PC1 (", round(pc.percent[1], 2), "%)")) +
    ylab(paste0("PC2 (", round(pc.percent[2], 2), "%)")) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Save the plot
  output_file <- file.path(output_dir, paste0("PCA_plot_with_labels.", output_format))
  ggsave(output_file, plot = plot, width = 10, height = 8, dpi = 300)
  
  cat("PCA plot saved as:", output_file, "\n")
}

# Run the script
if (!interactive()) {
  # Check and load required libraries
  check_and_load_libraries()
  
  # Parse and validate arguments
  opt <- parse_arguments()
  
  # Run main function
  main(opt$input, opt$output, opt$format)
}
