# Sequence utilities for gene amino acid sequences
# Functions to fetch and manage protein sequences for HGNC genes

#' Fetch amino acid sequences for HGNC genes
#' @param gene_symbols character vector of HGNC gene symbols
#' @return data.frame with gene symbols and their amino acid sequences
fetch_gene_sequences <- function(gene_symbols) {
    tryCatch({
        if (!requireNamespace("biomaRt", quietly = TRUE)) {
            stop("biomaRt package is required. Please install it with: BiocManager::install('biomaRt')")
        }
        
        cat("Fetching amino acid sequences for", length(gene_symbols), "genes...\n")
        
        # Connect to Ensembl BioMart
        ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        
        # Query BioMart to get sequences
        # We'll get the peptide sequences (amino acid sequences)
        sequences_result <- biomaRt::getBM(
            attributes = c("hgnc_symbol", "ensembl_gene_id", "ensembl_peptide_id", "peptide"),
            filters = "hgnc_symbol",
            values = gene_symbols,
            mart = ensembl
        )
        
        cat("BioMart returned", nrow(sequences_result), "sequence records\n")
        
        # For genes with multiple isoforms, select the longest sequence
        if (nrow(sequences_result) > 0) {
            # Add sequence length
            sequences_result$sequence_length <- nchar(sequences_result$peptide)
            
            # Remove empty sequences
            sequences_result <- sequences_result[!is.na(sequences_result$peptide) & 
                                               sequences_result$peptide != "", ]
            
            if (nrow(sequences_result) > 0) {
                # For each gene, keep the longest sequence
                sequences_result <- sequences_result[order(sequences_result$hgnc_symbol, 
                                                         -sequences_result$sequence_length), ]
                sequences_final <- sequences_result[!duplicated(sequences_result$hgnc_symbol), ]
                
                cat("Selected longest sequences for", nrow(sequences_final), "unique genes\n")
                
                return(sequences_final)
            }
        }
        
        # Return empty dataframe if no sequences found
        return(data.frame(
            hgnc_symbol = character(0),
            ensembl_gene_id = character(0),
            ensembl_peptide_id = character(0),
            peptide = character(0),
            sequence_length = numeric(0),
            stringsAsFactors = FALSE
        ))
        
    }, error = function(e) {
        cat("Error fetching sequences:", e$message, "\n")
        return(data.frame(
            hgnc_symbol = character(0),
            ensembl_gene_id = character(0),
            ensembl_peptide_id = character(0),
            peptide = character(0),
            sequence_length = numeric(0),
            stringsAsFactors = FALSE
        ))
    })
}

#' Load and cache all sequences for phylomap genes
#' @return data.frame with all sequences for genes in phylomap
load_phylomap_sequences <- function() {
    # Define cache file path
    sequences_cache_file <- file.path("data", "phylomap_sequences.rds")
    
    # Check if cache exists and is recent (less than 7 days old)
    if (file.exists(sequences_cache_file)) {
        file_age <- difftime(Sys.time(), file.mtime(sequences_cache_file), units = "days")
        if (file_age < 7) {
            cat("Loading cached sequence data (", round(file_age, 1), "days old)...\n")
            sequences_data <- readRDS(sequences_cache_file)
            cat("Loaded", nrow(sequences_data), "cached sequences\n")
            return(sequences_data)
        } else {
            cat("Cache is", round(file_age, 1), "days old, refreshing...\n")
        }
    }
    
    # Load phylomap to get all gene symbols
    phylomap <- load_phylomap()
    if (is.null(phylomap)) {
        stop("Could not load phylomap data")
    }
    
    all_genes <- unique(phylomap$GeneID)
    cat("Found", length(all_genes), "unique genes in phylomap\n")
    
    # Fetch sequences in batches to avoid timeouts
    batch_size <- 500
    all_sequences <- data.frame()
    
    for (i in seq(1, length(all_genes), batch_size)) {
        end_idx <- min(i + batch_size - 1, length(all_genes))
        batch_genes <- all_genes[i:end_idx]
        
        cat("Processing batch", ceiling(i/batch_size), "of", ceiling(length(all_genes)/batch_size), 
            "(genes", i, "to", end_idx, ")...\n")
        
        batch_sequences <- fetch_gene_sequences(batch_genes)
        all_sequences <- rbind(all_sequences, batch_sequences)
        
        # Small delay between batches to be respectful to the server
        if (end_idx < length(all_genes)) {
            Sys.sleep(2)
        }
    }
    
    # Save to cache
    if (nrow(all_sequences) > 0) {
        cat("Saving", nrow(all_sequences), "sequences to cache...\n")
        saveRDS(all_sequences, sequences_cache_file)
        
        # Also save as TSV for human inspection
        tsv_file <- file.path("data", "phylomap_sequences.tsv")
        write.table(all_sequences, tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)
        cat("Also saved as TSV file:", tsv_file, "\n")
    }
    
    return(all_sequences)
}

#' Get sequence information for a specific gene
#' @param gene_symbol HGNC gene symbol
#' @param sequences_data cached sequence data (optional)
#' @return list with sequence information
get_gene_sequence_info <- function(gene_symbol, sequences_data = NULL) {
    if (is.null(sequences_data)) {
        sequences_data <- load_phylomap_sequences()
    }
    
    # Find the gene in the sequence data
    gene_seq <- sequences_data[sequences_data$hgnc_symbol == gene_symbol, ]
    
    if (nrow(gene_seq) == 0) {
        return(list(
            found = FALSE,
            message = paste("No sequence data available for", gene_symbol)
        ))
    }
    
    # Get the sequence
    sequence <- gene_seq$peptide[1]
    sequence_length <- gene_seq$sequence_length[1]
    ensembl_gene_id <- gene_seq$ensembl_gene_id[1]
    ensembl_peptide_id <- gene_seq$ensembl_peptide_id[1]
    
    # Calculate some basic sequence properties
    aa_composition <- table(strsplit(sequence, "")[[1]])
    
    # Most common amino acids
    common_aa <- names(sort(aa_composition, decreasing = TRUE))[1:3]
    
    # Basic properties
    has_signal_peptide <- grepl("^M", sequence)  # Starts with methionine
    
    return(list(
        found = TRUE,
        gene_symbol = gene_symbol,
        sequence = sequence,
        sequence_length = sequence_length,
        ensembl_gene_id = ensembl_gene_id,
        ensembl_peptide_id = ensembl_peptide_id,
        aa_composition = aa_composition,
        most_common_aa = paste(common_aa, collapse = ", "),
        starts_with_met = has_signal_peptide,
        preview = paste0(substr(sequence, 1, 50), if(nchar(sequence) > 50) "..." else "")
    ))
}

#' Format sequence information for display
#' @param sequence_info list from get_gene_sequence_info
#' @return formatted string for display
format_sequence_info <- function(sequence_info) {
    if (!sequence_info$found) {
        return(sequence_info$message)
    }
    
    paste0(
        "=== PROTEIN SEQUENCE INFORMATION ===", "\n",
        "Ensembl Gene ID: ", sequence_info$ensembl_gene_id, "\n",
        "Ensembl Peptide ID: ", sequence_info$ensembl_peptide_id, "\n",
        "Sequence Length: ", sequence_info$sequence_length, " amino acids", "\n",
        "Most Common AA: ", sequence_info$most_common_aa, "\n",
        "Starts with Met: ", ifelse(sequence_info$starts_with_met, "Yes", "No"), "\n",
        "Sequence Preview: ", sequence_info$preview, "\n"
    )
}

#' Initialize sequence data (run on app startup)
#' @return invisible success status
initialize_sequence_data <- function() {
    tryCatch({
        cat("Initializing sequence data...\n")
        sequences <- load_phylomap_sequences()
        cat("Successfully initialized sequence data for", nrow(sequences), "genes\n")
        return(invisible(TRUE))
    }, error = function(e) {
        warning("Failed to initialize sequence data: ", e$message)
        return(invisible(FALSE))
    })
}
