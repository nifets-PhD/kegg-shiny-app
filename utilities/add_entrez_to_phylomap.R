#!/usr/bin/env Rscript

# Script to add Entrez IDs to the phylomap.tsv file

cat("=== Adding Entrez IDs to Phylomap ===\n")

# Load the current phylomap
cat("Loading current phylomap...\n")
phylomap <- read.table("data/phylomap.tsv", sep="\t", header=TRUE, 
                       stringsAsFactors=FALSE, quote='"', comment.char="")
cat("Current phylomap has", nrow(phylomap), "entries\n")
cat("Current columns:", paste(colnames(phylomap), collapse=", "), "\n")

# The HGNC column should be available as "HGNC"
hgnc_col <- "HGNC"

# Load gene mapping data
cat("Loading gene ID mapping...\n")
gene_mapping <- readRDS("data/gene_id_mapping.rds")
cat("Gene mapping has", nrow(gene_mapping), "entries\n")

# Add Entrez IDs to phylomap
cat("Adding Entrez IDs via HGNC symbol matching...\n")

# Create a clean Entrez mapping from our gene data
entrez_mapping <- gene_mapping[!is.na(gene_mapping$hgnc_symbol) & 
                              !is.na(gene_mapping$entrezgene_id), 
                              c("hgnc_symbol", "entrezgene_id")]

# Remove duplicates, keeping the first match for each HGNC symbol
entrez_mapping <- entrez_mapping[!duplicated(entrez_mapping$hgnc_symbol), ]

cat("Created clean Entrez mapping with", nrow(entrez_mapping), "entries\n")

# Merge with phylomap using the correct column name
enhanced_phylomap <- merge(phylomap, entrez_mapping, 
                          by.x = hgnc_col, by.y = "hgnc_symbol", 
                          all.x = TRUE)

# Rename the column for clarity
colnames(enhanced_phylomap)[colnames(enhanced_phylomap) == "entrezgene_id"] <- "Entrez_ID"

# For entries that didn't get Entrez IDs via HGNC, try UniProt matching
# The UniProt column is "UniProt"
uniprot_col <- "UniProt"

missing_entrez <- is.na(enhanced_phylomap$Entrez_ID) & !is.na(enhanced_phylomap[[uniprot_col]])
cat("Attempting to fill", sum(missing_entrez), "missing Entrez IDs via UniProt matching...\n")

if (sum(missing_entrez) > 0) {
    # Create UniProt to Entrez mapping
    uniprot_to_entrez <- gene_mapping[!is.na(gene_mapping$uniprotswissprot) & 
                                     !is.na(gene_mapping$entrezgene_id), 
                                     c("uniprotswissprot", "entrezgene_id")]
    
    # Remove duplicates
    uniprot_to_entrez <- uniprot_to_entrez[!duplicated(uniprot_to_entrez$uniprotswissprot), ]
    
    # Match missing entries
    for (i in which(missing_entrez)) {
        uniprot_id <- enhanced_phylomap[[uniprot_col]][i]
        if (!is.na(uniprot_id) && uniprot_id != "") {
            entrez_match <- uniprot_to_entrez$entrezgene_id[uniprot_to_entrez$uniprotswissprot == uniprot_id]
            
            if (length(entrez_match) > 0) {
                enhanced_phylomap$Entrez_ID[i] <- entrez_match[1]
            }
        }
    }
}

# Final statistics
total_with_entrez <- sum(!is.na(enhanced_phylomap$Entrez_ID))
coverage_pct <- round(100 * total_with_entrez / nrow(enhanced_phylomap), 1)

cat("=== ENHANCEMENT RESULTS ===\n")
cat("Total phylomap entries:", nrow(enhanced_phylomap), "\n")
cat("Entries with Entrez IDs:", total_with_entrez, "\n")
cat("Entrez ID coverage:", coverage_pct, "%\n")
cat("Final columns:", paste(colnames(enhanced_phylomap), collapse=", "), "\n")

# Test key genes
cat("\n=== KEY GENE EXAMPLES ===\n")
test_genes <- c("TP53", "BRCA1", "EGFR", "MYC")
for (gene in test_genes) {
    match_row <- enhanced_phylomap[enhanced_phylomap[[hgnc_col]] == gene, ]
    if (nrow(match_row) > 0) {
        cat(sprintf("%s: Entrez=%s, UniProt=%s, Stratum=%s\n", 
                    gene, 
                    ifelse(is.na(match_row$Entrez_ID[1]), "NA", match_row$Entrez_ID[1]),
                    ifelse(is.na(match_row[[uniprot_col]][1]), "NA", match_row[[uniprot_col]][1]),
                    match_row$Stratum[1]))
    }
}

# Save the enhanced phylomap
cat("\nSaving enhanced phylomap to data/phylomap.tsv...\n")
write.table(enhanced_phylomap, "data/phylomap.tsv", sep="\t", row.names=FALSE, quote=FALSE)

cat("Enhancement completed successfully!\n")
