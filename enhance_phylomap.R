#!/usr/bin/env Rscript

#' Enhance phylomap.tsv with UniProt ID mappings
#' This script adds UniProt IDs to the phylomap for better gene ID compatibility

cat("Enhancing phylomap with UniProt ID mappings...\n")

# Load the current phylomap
cat("Loading current phylomap...\n")
phylomap <- read.table("data/phylomap.tsv", sep = "\t", header = TRUE, 
                       stringsAsFactors = FALSE, quote = '"')
cat("Loaded phylomap with", nrow(phylomap), "entries\n")

# Load the comprehensive gene ID mapping
cat("Loading gene ID mapping...\n")
if (file.exists("data/gene_id_mapping.rds")) {
    gene_mapping <- readRDS("data/gene_id_mapping.rds")
} else {
    stop("Gene ID mapping not found. Please run download_gene_mapping.R first.")
}
cat("Loaded gene mapping with", nrow(gene_mapping), "entries\n")

# Clean the gene mapping data
gene_mapping$hgnc_symbol[gene_mapping$hgnc_symbol == ""] <- NA
gene_mapping$uniprotswissprot[gene_mapping$uniprotswissprot == ""] <- NA
gene_mapping$uniprotswissprot[gene_mapping$uniprotswissprot == "NA"] <- NA

# Create a lookup table for HGNC -> UniProt mapping
# Prioritize entries with UniProt IDs and remove duplicates
hgnc_to_uniprot <- gene_mapping[!is.na(gene_mapping$hgnc_symbol) & 
                               !is.na(gene_mapping$uniprotswissprot), 
                               c("hgnc_symbol", "uniprotswissprot")]

# Remove duplicates, keeping first occurrence (should be highest priority from our deduplication)
hgnc_to_uniprot <- hgnc_to_uniprot[!duplicated(hgnc_to_uniprot$hgnc_symbol), ]

cat("Created HGNC->UniProt lookup with", nrow(hgnc_to_uniprot), "mappings\n")

# Match phylomap HGNC symbols to UniProt IDs
cat("Matching phylomap HGNC symbols to UniProt IDs...\n")

# Create the enhanced phylomap by adding UniProt column
phylomap_enhanced <- phylomap

# Add UniProt column by matching HGNC symbols
phylomap_enhanced$UniProt <- NA

# Match using HGNC symbols
matches <- match(phylomap_enhanced$HGNC, hgnc_to_uniprot$hgnc_symbol)
valid_matches <- !is.na(matches)
phylomap_enhanced$UniProt[valid_matches] <- hgnc_to_uniprot$uniprotswissprot[matches[valid_matches]]

# Count successful matches
uniprot_matches <- sum(!is.na(phylomap_enhanced$UniProt))
cat("Successfully matched", uniprot_matches, "out of", nrow(phylomap_enhanced), 
    "phylomap entries to UniProt IDs\n")
cat("Match rate:", round(100 * uniprot_matches / nrow(phylomap_enhanced), 1), "%\n")

# Show some examples
cat("\nExample matches:\n")
examples <- phylomap_enhanced[!is.na(phylomap_enhanced$UniProt) & 
                              phylomap_enhanced$HGNC %in% c("TP53", "BRCA1", "EGFR", "MYC"), 
                              c("HGNC", "UniProt", "Stratum", "Name")]
if (nrow(examples) > 0) {
    print(examples)
} else {
    # Show first few matches
    examples <- head(phylomap_enhanced[!is.na(phylomap_enhanced$UniProt), 
                                       c("HGNC", "UniProt", "Stratum", "Name")], 5)
    print(examples)
}

# Save the enhanced phylomap
cat("\nSaving enhanced phylomap...\n")

# Backup original
if (file.exists("data/phylomap.tsv") && !file.exists("data/phylomap_original.tsv")) {
    file.copy("data/phylomap.tsv", "data/phylomap_original.tsv")
    cat("Backed up original phylomap to phylomap_original.tsv\n")
}

# Save enhanced version
output_file <- "data/phylomap.tsv"
write.table(phylomap_enhanced, output_file, sep = "\t", row.names = FALSE, 
            quote = TRUE, na = "")
cat("Saved enhanced phylomap with UniProt IDs to", output_file, "\n")

# Summary statistics
cat("\n=== ENHANCED PHYLOMAP SUMMARY ===\n")
cat("Total entries:", nrow(phylomap_enhanced), "\n")
cat("Entries with HGNC symbols:", sum(!is.na(phylomap_enhanced$HGNC) & phylomap_enhanced$HGNC != ""), "\n")
cat("Entries with UniProt IDs:", sum(!is.na(phylomap_enhanced$UniProt)), "\n")
cat("Entries with Ensembl Protein IDs:", sum(!is.na(phylomap_enhanced$ENSP)), "\n")
cat("Entries with Ensembl Gene IDs:", sum(!is.na(phylomap_enhanced$ENSG)), "\n")

# Show phylostratum distribution
cat("\nPhylostratum distribution:\n")
stratum_counts <- table(phylomap_enhanced$Stratum)
print(stratum_counts)

cat("\nPhylomap enhancement completed successfully!\n")
