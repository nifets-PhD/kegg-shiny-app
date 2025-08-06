#!/usr/bin/env Rscript

# Script to download and cache comprehensive gene ID mapping
# Run this script to pre-populate the mapping cache with Entrez, HGNC, Ensembl, and UniProt IDs

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

cat("=== Comprehensive Gene ID Mapping Download ===\n")
cat("This will download the complete mapping including:\n")
cat("- HGNC symbols\n")
cat("- Entrez Gene IDs\n")
cat("- Ensembl Gene IDs\n")
cat("- UniProt IDs (SwissProt & TrEMBL)\n")
cat("from Ensembl BioMart and cache it locally.\n\n")

# Download the mapping
cat("Starting download...\n")
mapping_data <- download_hgnc_uniprot_mapping(force_refresh = TRUE)

if (nrow(mapping_data) > 0) {
    cat("\n=== Download Summary ===\n")
    cat("Total genes in comprehensive mapping:", nrow(mapping_data), "\n")
    cat("HGNC symbols:", sum(!is.na(mapping_data$hgnc_symbol)), "\n")
    cat("Entrez IDs:", sum(!is.na(mapping_data$entrez_id)), "\n")
    cat("Ensembl Gene IDs:", sum(!is.na(mapping_data$ensembl_gene_id)), "\n")
    cat("UniProt IDs:", sum(!is.na(mapping_data$uniprot_id)), "\n")
    cat("SwissProt entries:", sum(!is.na(mapping_data$uniprot_swissprot) & mapping_data$uniprot_swissprot != "", na.rm = TRUE), "\n")
    cat("TrEMBL entries:", sum(!is.na(mapping_data$uniprot_trembl) & mapping_data$uniprot_trembl != "", na.rm = TRUE), "\n")
    
    # Show some examples across all ID types
    cat("\n=== Sample Mappings ===\n")
    sample_genes <- c("TP53", "BRCA1", "EGFR", "MYC", "GAPDH")
    for (gene in sample_genes) {
        gene_row <- mapping_data[mapping_data$hgnc_symbol == gene & !is.na(mapping_data$hgnc_symbol), ]
        if (nrow(gene_row) > 0) {
            cat(sprintf("%-6s -> Entrez: %-8s Ensembl: %-18s UniProt: %-10s\n", 
                       gene, 
                       ifelse(is.na(gene_row$entrez_id[1]), "N/A", gene_row$entrez_id[1]),
                       ifelse(is.na(gene_row$ensembl_gene_id[1]), "N/A", gene_row$ensembl_gene_id[1]),
                       ifelse(is.na(gene_row$uniprot_id[1]), "N/A", gene_row$uniprot_id[1])))
        } else {
            cat(gene, "-> Not found\n")
        }
    }
    
    cat("\n=== Files Created ===\n")
    cat("Comprehensive cache: data/comprehensive_gene_mapping.rds\n")
    cat("Comprehensive TSV: data/comprehensive_gene_mapping.tsv\n")
    cat("Legacy cache: data/hgnc_uniprot_mapping.rds (backward compatibility)\n")
    cat("Legacy TSV: data/hgnc_uniprot_mapping.tsv\n")
    
    cat("\nDownload completed successfully!\n")
    cat("This mapping will now support all gene ID types in the Shiny app:\n")
    cat("- Entrez Gene IDs (recommended for KEGG)\n")
    cat("- HGNC gene symbols\n")
    cat("- Ensembl Gene IDs\n")
    cat("- UniProt IDs\n")
} else {
    cat("ERROR: Download failed or no data retrieved.\n")
    quit(status = 1)
}
