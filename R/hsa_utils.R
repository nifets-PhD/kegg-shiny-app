#' HSA (Human Gene) Utilities for KEGG Pathways
#' Functions for extracting, caching, and fetching data for HSA genes in pathways

library(biomaRt)
library(dplyr)
library(graph)  # For working with graphNEL objects

#' Extract HSA codes from KEGG pathway data
#' @param kegg_data either parsed KEGG XML data from parseKGML or graphNEL object
#' @return character vector of unique HSA codes
extract_hsa_codes <- function(kegg_data) {
    if (is.null(kegg_data)) {
        return(character(0))
    }
    
    hsa_codes <- character(0)
    
    # Handle different data types
    if ("graph" %in% class(kegg_data) || "graphNEL" %in% class(kegg_data)) {
        # If it's a graph object, extract from node names
        node_names <- nodes(kegg_data)
        cat("Extracting HSA codes from", length(node_names), "graph nodes\n")
        
        # Look for HSA patterns in node names
        for (node_name in node_names) {
            if (grepl("hsa:", node_name)) {
                # Extract HSA codes from node name
                codes <- unlist(strsplit(node_name, "\\s+"))
                hsa_matches <- codes[grepl("^hsa:", codes)]
                numeric_ids <- gsub("^hsa:", "", hsa_matches)
                hsa_codes <- c(hsa_codes, numeric_ids)
            }
        }
        
        # Debug: show sample node names
        hsa_nodes <- node_names[grepl("hsa:", node_names)]
        if (length(hsa_nodes) > 0) {
            cat("Sample HSA node names:\n")
            cat(paste(head(hsa_nodes, 3), collapse = "\n"), "\n")
        }
    } else if (is.list(kegg_data) && "nodes" %in% names(kegg_data)) {
        # If it's parsed XML data with nodes
        node_names <- kegg_data$nodes$name
        cat("Extracting HSA codes from", length(node_names), "XML nodes\n")
        
        # Filter for entries that contain "hsa:" 
        hsa_entries <- node_names[grepl("hsa:", node_names)]
        
        # Extract individual HSA codes
        for (entry in hsa_entries) {
            codes <- unlist(strsplit(entry, "\\s+"))
            hsa_matches <- codes[grepl("^hsa:", codes)]
            numeric_ids <- gsub("^hsa:", "", hsa_matches)
            hsa_codes <- c(hsa_codes, numeric_ids)
        }
    } else {
        cat("Unknown data type for HSA extraction:", class(kegg_data), "\n")
        return(character(0))
    }
    
    # Return unique codes
    unique_hsa_codes <- unique(hsa_codes)
    cat("Extracted", length(unique_hsa_codes), "unique HSA codes from pathway\n")
    
    # Debug: show sample codes
    if (length(unique_hsa_codes) > 0) {
        cat("Sample HSA codes:", paste(head(unique_hsa_codes, 5), collapse = ", "), "\n")
    }
    
    return(unique_hsa_codes)
}

#' Create HSA cache directory structure
#' @param base_dir base directory for cache
#' @return path to HSA cache directory
create_hsa_cache <- function(base_dir = ".") {
    hsa_cache_dir <- file.path(base_dir, "hsa_cache")
    
    if (!dir.exists(hsa_cache_dir)) {
        dir.create(hsa_cache_dir, recursive = TRUE)
        cat("Created HSA cache directory:", hsa_cache_dir, "\n")
    }
    
    return(hsa_cache_dir)
}

#' Load cached HSA data if available
#' @param hsa_codes character vector of HSA codes
#' @param cache_dir HSA cache directory
#' @return list with cached_data and missing_codes
load_cached_hsa_data <- function(hsa_codes, cache_dir) {
    cached_data <- data.frame()
    missing_codes <- character(0)
    
    for (code in hsa_codes) {
        cache_file <- file.path(cache_dir, paste0("hsa_", code, ".rds"))
        
        if (file.exists(cache_file)) {
            tryCatch({
                code_data <- readRDS(cache_file)
                cached_data <- rbind(cached_data, code_data)
            }, error = function(e) {
                cat("Error loading cached data for", code, ":", e$message, "\n")
                missing_codes <- c(missing_codes, code)
            })
        } else {
            missing_codes <- c(missing_codes, code)
        }
    }
    
    cat("Loaded", nrow(cached_data), "cached entries,", length(missing_codes), "codes need fetching\n")
    
    return(list(
        cached_data = cached_data,
        missing_codes = missing_codes
    ))
}

#' Fetch comprehensive gene data for HSA codes from BioMart
#' @param hsa_codes character vector of HSA numeric codes
#' @param cache_dir directory to cache results
#' @return data.frame with gene information
fetch_hsa_gene_data <- function(hsa_codes, cache_dir) {
    if (length(hsa_codes) == 0) {
        return(data.frame())
    }
    
    cat("Fetching comprehensive data for", length(hsa_codes), "HSA codes...\n")
    
    tryCatch({
        # Connect to Ensembl
        ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        
        cat("Testing BioMart connection with sample HSA codes...\n")
        
        # Test with just the first few codes
        test_codes <- head(hsa_codes, 5)
        
        # Basic connectivity test
        gene_data <- getBM(
            attributes = c("entrezgene_id", "hgnc_symbol"),
            filters = "entrezgene_id",
            values = test_codes,
            mart = ensembl
        )
        
        if (nrow(gene_data) > 0) {
            cat("BioMart connection successful! Now fetching all", length(hsa_codes), "HSA codes...\n")
        } else {
            cat("BioMart test failed - no results for test codes\n")
            return(data.frame())
        }
        
        cat("BioMart connection successful! Now fetching all", length(hsa_codes), "HSA codes...\n")
        
        # Fetch basic gene information for all codes
        gene_data <- getBM(
            attributes = c(
                "entrezgene_id", "hgnc_symbol", "external_gene_name", 
                "description", "gene_biotype", "chromosome_name", 
                "start_position", "end_position", "strand"
            ),
            filters = "entrezgene_id",
            values = hsa_codes,
            mart = ensembl
        )
        
        cat("BioMart returned", nrow(gene_data), "gene records\n")
        
        # Now fetch protein sequences separately using a different approach
        if (nrow(gene_data) > 0) {
            cat("Fetching protein sequences...\n")
            
            # Try to get sequences using HGNC symbols for genes that have them
            genes_with_symbols <- gene_data[!is.na(gene_data$hgnc_symbol) & gene_data$hgnc_symbol != "", ]
            
            if (nrow(genes_with_symbols) > 0) {
                tryCatch({
                    sequence_data <- getBM(
                        attributes = c("hgnc_symbol", "peptide"),
                        filters = "hgnc_symbol",
                        values = genes_with_symbols$hgnc_symbol,
                        mart = ensembl
                    )
                    
                    cat("Retrieved", nrow(sequence_data), "protein sequences\n")
                    
                    # Merge sequences back into gene_data
                    gene_data <- merge(gene_data, sequence_data, by = "hgnc_symbol", all.x = TRUE)
                    
                }, error = function(e) {
                    cat("Could not fetch sequences via HGNC symbols:", e$message, "\n")
                    gene_data$peptide <<- NA_character_
                })
            } else {
                gene_data$peptide <- NA_character_
            }
        }
        
        cat("BioMart returned", nrow(gene_data), "records\n")
        
        # Process and clean the data
        if (nrow(gene_data) > 0) {
            # Calculate peptide lengths and add missing columns
            gene_data <- gene_data %>%
                mutate(
                    peptide_length = nchar(peptide),
                    protein_id = NA_character_,     # Placeholder for protein ID
                    hsa_code = as.character(entrezgene_id)
                ) %>%
                group_by(entrezgene_id) %>%
                # If there are duplicates, select the one with the longest sequence
                slice_max(order_by = coalesce(peptide_length, 0), n = 1, with_ties = FALSE) %>%
                ungroup()
            
            # Cache individual entries
            for (i in 1:nrow(gene_data)) {
                row_data <- gene_data[i, ]
                hsa_id <- row_data$hsa_code
                
                if (!is.na(hsa_id)) {
                    cache_file <- file.path(cache_dir, paste0("hsa_", hsa_id, ".rds"))
                    tryCatch({
                        saveRDS(row_data, cache_file)
                    }, error = function(e) {
                        cat("Warning: Could not cache data for HSA", hsa_id, ":", e$message, "\n")
                    })
                }
            }
            
            cat("Successfully processed and cached", nrow(gene_data), "gene records\n")
        }
        
        return(gene_data)
        
    }, error = function(e) {
        cat("Error fetching HSA gene data:", e$message, "\n")
        return(data.frame())
    })
}

#' Get comprehensive HSA gene data for a pathway
#' @param kegg_data parsed KEGG pathway XML data
#' @param cache_dir cache directory (will be created if needed)
#' @return data.frame with comprehensive gene information
get_pathway_hsa_data <- function(kegg_data, cache_dir = NULL) {
    if (is.null(cache_dir)) {
        cache_dir <- create_hsa_cache()
    }
    
    # Extract HSA codes from pathway
    hsa_codes <- extract_hsa_codes(kegg_data)
    
    if (length(hsa_codes) == 0) {
        cat("No HSA codes found in pathway\n")
        return(data.frame())
    }
    
    # Load cached data and identify missing codes
    cache_result <- load_cached_hsa_data(hsa_codes, cache_dir)
    
    # Fetch missing data
    new_data <- data.frame()
    if (length(cache_result$missing_codes) > 0) {
        new_data <- fetch_hsa_gene_data(cache_result$missing_codes, cache_dir)
    }
    
    # Combine cached and new data
    all_data <- rbind(cache_result$cached_data, new_data)
    
    if (nrow(all_data) > 0) {
        cat("Retrieved complete data for", nrow(all_data), "genes in pathway\n")
        
        # Add some summary statistics
        cat("Gene types found:", paste(unique(all_data$gene_biotype), collapse = ", "), "\n")
        cat("Chromosomes represented:", length(unique(all_data$chromosome_name)), "\n")
        cat("Genes with sequences:", sum(!is.na(all_data$peptide) & all_data$peptide != ""), "\n")
    }
    
    return(all_data)
}

#' Create summary table of HSA genes for display
#' @param hsa_data data.frame from get_pathway_hsa_data
#' @return data.frame formatted for display
create_hsa_summary_table <- function(hsa_data) {
    if (nrow(hsa_data) == 0) {
        return(data.frame(
            HSA_ID = character(0),
            Gene_Symbol = character(0),
            Gene_Name = character(0),
            Chromosome = character(0),
            Gene_Type = character(0),
            Has_Sequence = logical(0)
        ))
    }
    
    summary_table <- hsa_data %>%
        select(
            HSA_ID = hsa_code,
            Gene_Symbol = hgnc_symbol,
            Gene_Name = external_gene_name,
            Description = description,
            Chromosome = chromosome_name,
            Gene_Type = gene_biotype,
            Sequence_Length = peptide_length
        ) %>%
        mutate(
            Gene_Symbol = ifelse(is.na(Gene_Symbol) | Gene_Symbol == "", HSA_ID, Gene_Symbol),
            Gene_Name = ifelse(is.na(Gene_Name) | Gene_Name == "", "Unknown", Gene_Name),
            Has_Sequence = !is.na(Sequence_Length) & Sequence_Length > 0,
            Description = substr(Description, 1, 100)  # Truncate long descriptions
        ) %>%
        arrange(Gene_Symbol)
    
    return(summary_table)
}

#' Get sequence information for a specific HSA code
#' @param hsa_code HSA numeric code
#' @param hsa_data complete HSA data from get_pathway_hsa_data
#' @return list with sequence information
get_hsa_sequence_info <- function(hsa_code, hsa_data) {
    if (nrow(hsa_data) == 0) {
        return(NULL)
    }
    
    gene_info <- hsa_data[hsa_data$hsa_code == hsa_code, ]
    
    if (nrow(gene_info) == 0) {
        return(NULL)
    }
    
    gene_info <- gene_info[1, ]  # Take first match
    
    sequence_info <- list(
        hsa_id = gene_info$hsa_code,
        gene_symbol = ifelse(is.na(gene_info$hgnc_symbol), gene_info$hsa_code, gene_info$hgnc_symbol),
        gene_name = gene_info$external_gene_name,
        description = gene_info$description,
        chromosome = gene_info$chromosome_name,
        gene_type = gene_info$gene_biotype,
        sequence = gene_info$peptide,
        sequence_length = gene_info$peptide_length,
        protein_id = gene_info$protein_id
    )
    
    return(sequence_info)
}
