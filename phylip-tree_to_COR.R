#!/usr/bin/env Rscript

## phylip-tree_to_COR.R
#
# This script reads a phylogenetic tree file, computes pairwise distances,
# and generates a correlation heatmap with clustered samples.
#
# Usage:
#   Rscript phylip-tree_to_COR.R -i <input_tree_file> -o <output_directory> [-f <output_format>]
#
# Example:
#   Rscript phylip-tree_to_COR.R -i input.treefile -o /path/to/output/
#   Rscript phylip-tree_to_COR.R -i input.treefile -o /path/to/output/ -f pdf
#
# Author: SP@NC (+AI)
# Date: 2025-02-25
# Version: 1.1

# Function to check and load required libraries
check_and_load_libraries <- function() {
  required_packages <- c("ape", "pheatmap", "optparse")
  
  for (package in required_packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(paste("Error: Required package", package, "is not installed. Please install it using install.packages('", package, "')"), call. = FALSE)
    }
  }
  
  # Load the libraries
  library(ape)
  library(pheatmap)
  library(optparse)
}

# Parse command line arguments
parse_arguments <- function() {
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="Input tree file", metavar="FILE"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output directory", metavar="DIRECTORY"),
    make_option(c("-f", "--format"), type="character", default="png",
                help="Output format (png or pdf) [default= %default]", metavar="FORMAT")
  )
  
  opt_parser <- OptionParser(option_list=option_list,
                             description="Generate correlation heatmap from phylogenetic tree")
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
  
  # Read the tree file
  tree <- read.tree(input_file)
  
  # Compute pairwise distances
  distances <- cophenetic(tree)
  
  # Create and save the heatmap
  output_file <- file.path(output_dir, paste0("correlation_heatmap.", output_format))
  
  pheatmap(distances, 
           clustering_distance_rows = as.dist(distances),
           clustering_distance_cols = as.dist(distances),
           main = "Pairwise Distances",
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 8,
           fontsize_col = 8,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           filename = output_file,
           width = 12,
           height = 10)
  
  cat("Correlation heatmap saved as:", output_file, "\n")
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
