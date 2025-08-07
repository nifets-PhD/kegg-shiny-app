#!/usr/bin/env Rscript
# Download comprehensive gene ID mapping from Ensembl BioMart
# Run this script to update the gene mapping data file

library(biomaRt)
library(dplyr)

cat("Starting comprehensive gene ID mapping download from Ensembl...\n")

# Create data directory if it doesn't exist
if (!dir.exists("data")) {
    dir.create("data", recursive = TRUE)
}

# Function to download mapping with retry logic
download_gene_mapping <- function() {
    tryCatch({
        # Function to merge UniProt IDs for duplicate rows
        merge_uniprot <- function(uniprot_ids) {
            # Remove NA and empty strings
            clean_ids <- uniprot_ids[!is.na(uniprot_ids) & uniprot_ids != ""]
            
            if (length(clean_ids) == 0) {
                return(NA_character_)
            } else if (length(clean_ids) == 1) {
                return(clean_ids[1])
            } else {
                # If multiple non-empty UniProt IDs, join them with semicolon
                return(paste(unique(clean_ids), collapse = ";"))
            }
        }
        
        # Connect to Ensembl BioMart with retry logic
        cat("Connecting to Ensembl BioMart...\n")
        
        ensembl <- NULL
        max_retries <- 3
        for (retry_count in 1:max_retries) {
            ensembl <- tryCatch({
                biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
            }, error = function(e) {
                if (retry_count < max_retries) {
                    cat("Connection attempt", retry_count, "failed, retrying in 3 seconds...\n")
                    Sys.sleep(3)
                    NULL
                } else {
                    warning("Failed to connect to BioMart after ", max_retries, " attempts: ", e$message)
                    NULL
                }
            })
            
            if (!is.null(ensembl)) break
        }
        
        if (is.null(ensembl)) {
            stop("Could not establish connection to Ensembl BioMart")
        }
        
        # Get comprehensive mapping including Entrez, HGNC, Ensembl, and UniProt IDs
        # Split into smaller queries to avoid BioMart limits
        cat("Querying BioMart for comprehensive gene mappings (in batches)...\n")
        
        # First query: Basic gene identifiers
        cat("Batch 1: Getting basic gene identifiers...\n")
        mapping_basic <- biomaRt::getBM(
            attributes = c("hgnc_symbol", "entrezgene_id", "ensembl_gene_id"),
            mart = ensembl
        )
        
        # Second query: UniProt identifiers for genes that have HGNC symbols
        cat("Batch 2: Getting UniProt mappings...\n")
        hgnc_genes <- mapping_basic$hgnc_symbol[!is.na(mapping_basic$hgnc_symbol) & mapping_basic$hgnc_symbol != ""]
        
        if (length(hgnc_genes) > 0) {
            # Query in smaller batches to avoid timeout and attribute limits
            batch_size <- 1000  # Conservative batch size
            uniprot_results <- list()
            
            for (i in seq(1, length(hgnc_genes), batch_size)) {
                end_idx <- min(i + batch_size - 1, length(hgnc_genes))
                batch_genes <- hgnc_genes[i:end_idx]
                
                cat("  UniProt batch", ceiling(i/batch_size), "of", ceiling(length(hgnc_genes)/batch_size), 
                    ":", length(batch_genes), "genes\n")
                
                batch_uniprot <- tryCatch({
                    biomaRt::getBM(
                        attributes = c("hgnc_symbol", "uniprotswissprot"),  # Only SwissProt for reliability
                        filters = "hgnc_symbol",
                        values = batch_genes,
                        mart = ensembl
                    )
                }, error = function(e) {
                    cat("    Warning: Batch failed, skipping:", e$message, "\n")
                    data.frame(hgnc_symbol = character(0), uniprotswissprot = character(0), 
                              stringsAsFactors = FALSE)
                })
                
                if (nrow(batch_uniprot) > 0) {
                    uniprot_results[[length(uniprot_results) + 1]] <- batch_uniprot
                }
                
                # Small delay between batches to be nice to the server
                if (i + batch_size <= length(hgnc_genes)) {
                    Sys.sleep(1)
                }
            }
            
            # Combine UniProt results
            if (length(uniprot_results) > 0) {
                uniprot_mapping <- do.call(rbind, uniprot_results)
            } else {
                uniprot_mapping <- data.frame(hgnc_symbol = character(0), uniprotswissprot = character(0), 
                                            stringsAsFactors = FALSE)
            }
        } else {
            uniprot_mapping <- data.frame(hgnc_symbol = character(0), uniprotswissprot = character(0), 
                                        stringsAsFactors = FALSE)
        }
        
        # Merge all results
        cat("Merging batch results...\n")
        mapping_result <- merge(mapping_basic, uniprot_mapping, by = "hgnc_symbol", all.x = TRUE)
        
        cat("Downloaded", nrow(mapping_result), "gene records from Ensembl\n")
        
        if (nrow(mapping_result) > 0) {
            cat("Processing", nrow(mapping_result), "records...\n")
            
            # Merge duplicate rows that differ only in UniProt ID
            cat("Merging duplicate rows...\n")
            original_count <- nrow(mapping_result)
            mapping_result <- mapping_result %>%
                group_by(hgnc_symbol, entrezgene_id, ensembl_gene_id) %>%
                summarise(
                    uniprotswissprot = merge_uniprot(uniprotswissprot),
                    .groups = "drop"
                )
            
            cat("Merged", original_count - nrow(mapping_result), "duplicate rows\n")
            cat("Kept", nrow(mapping_result), "unique records\n")
            
            # Filter for records with at least one valid ID
            has_hgnc <- !is.na(mapping_result$hgnc_symbol) & mapping_result$hgnc_symbol != ""
            has_entrez <- !is.na(mapping_result$entrezgene_id)
            has_ensembl <- !is.na(mapping_result$ensembl_gene_id) & mapping_result$ensembl_gene_id != ""
            has_uniprot <- !is.na(mapping_result$uniprotswissprot) & mapping_result$uniprotswissprot != ""
            
            valid_records <- has_hgnc | has_entrez | has_ensembl | has_uniprot
            mapping_result <- mapping_result[valid_records, ]
            
            cat("Kept", nrow(mapping_result), "records with valid identifiers\n")
            
            # Save to RDS file
            output_file <- "data/gene_id_mapping.rds"
            saveRDS(mapping_result, output_file)
            cat("Saved gene ID mapping to", output_file, "\n")
            
            # Also save as TSV for manual inspection
            output_tsv <- "data/gene_id_mapping.tsv"
            write.table(mapping_result, output_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
            cat("Also saved as", output_tsv, "for inspection\n")
            
            # Print summary statistics
            cat("\n=== MAPPING SUMMARY ===\n")
            cat("Total records:", nrow(mapping_result), "\n")
            cat("Records with HGNC symbols:", sum(has_hgnc[valid_records]), "\n")
            cat("Records with Entrez IDs:", sum(has_entrez[valid_records]), "\n")
            cat("Records with Ensembl IDs:", sum(has_ensembl[valid_records]), "\n")
            cat("Records with UniProt IDs:", sum(has_uniprot[valid_records]), "\n")
            cat("======================\n")
            
            return(mapping_result)
        } else {
            stop("No valid gene mapping data retrieved from Ensembl")
        }
        
    }, error = function(e) {
        stop("Error downloading comprehensive gene mapping: ", e$message)
    })
}

# Run the download
mapping_data <- download_gene_mapping()
cat("Gene ID mapping download completed successfully!\n")
