#!/usr/bin/env Rscript

# Script to download and cache HGNC to UniProt mapping
# Run this script to pre-populate the mapping cache

# Load required libraries
cat("Loading required libraries...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

if (!requireNamespace("biomaRt", quietly = TRUE)) {
    cat("Installing biomaRt...\n")
    BiocManager::install("biomaRt", ask = FALSE, update = FALSE)
}

library(biomaRt)

# Source the network utils to get our functions
source("R/network_utils.R")

cat("=== HGNC to UniProt Mapping Download ===\n")
cat("This will download the complete mapping from Ensembl BioMart\n")
cat("and cache it locally for faster access.\n\n")

# Download the mapping
cat("Starting download...\n")
mapping_data <- download_hgnc_uniprot_mapping(force_refresh = TRUE)

if (nrow(mapping_data) > 0) {
    cat("\n=== Download Summary ===\n")
    cat("Total genes with UniProt IDs:", nrow(mapping_data), "\n")
    cat("SwissProt entries:", sum(!is.na(mapping_data$uniprot_swissprot) & mapping_data$uniprot_swissprot != ""), "\n")
    cat("TrEMBL entries:", sum(!is.na(mapping_data$uniprot_trembl) & mapping_data$uniprot_trembl != ""), "\n")
    
    # Show some examples
    cat("\n=== Sample Mappings ===\n")
    sample_genes <- c("TP53", "BRCA1", "EGFR", "MYC", "GAPDH")
    for (gene in sample_genes) {
        uniprot_id <- mapping_data$uniprot_id[mapping_data$hgnc_symbol == gene]
        if (length(uniprot_id) > 0 && !is.na(uniprot_id)) {
            cat(gene, "->", uniprot_id, "\n")
        } else {
            cat(gene, "-> Not found\n")
        }
    }
    
    cat("\n=== Files Created ===\n")
    cat("Cache file: data/hgnc_uniprot_mapping.rds\n")
    cat("TSV file: data/hgnc_uniprot_mapping.tsv\n")
    
    cat("\nDownload completed successfully!\n")
} else {
    cat("ERROR: Download failed or no data retrieved.\n")
    quit(status = 1)
}
