#!/usr/bin/env Rscript

# Refresh the HGNC to UniProt mapping cache
# This script can be run periodically (e.g., monthly) to update the mapping

cat("=== UniProt Mapping Refresh ===\n")
cat("This will refresh the HGNC to UniProt mapping cache\n\n")

# Source the network utils to get our functions
source("R/network_utils.R")

# Force refresh the mapping
cat("Refreshing mapping cache...\n")
mapping_data <- download_hgnc_uniprot_mapping(force_refresh = TRUE)

if (nrow(mapping_data) > 0) {
    cat("\n=== Refresh Summary ===\n")
    cat("Total genes with UniProt IDs:", nrow(mapping_data), "\n")
    cat("Cache updated successfully!\n")
    cat("Files updated:\n")
    cat("  - data/hgnc_uniprot_mapping.rds\n")
    cat("  - data/hgnc_uniprot_mapping.tsv\n")
} else {
    cat("ERROR: Refresh failed.\n")
    quit(status = 1)
}
