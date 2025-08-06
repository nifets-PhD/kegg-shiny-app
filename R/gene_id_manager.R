# Centralized Gene ID Management System
# =====================================
# This module provides a clean, standardized approach to gene ID handling
# All gene IDs are internally converted to and stored as Entrez IDs
# This eliminates the need for complex conditional conversions throughout the app

#' Convert any gene ID type to Entrez IDs (internal standard)
#' 
#' @param gene_ids Vector of gene IDs
#' @param input_type Type of input IDs: "entrez", "symbol", "ensembl", "uniprot"
#' @param comprehensive_mapping Pre-loaded mapping data (optional)
#' @return List with entrez_ids, mapping_info, and conversion_stats
convert_to_internal_standard <- function(gene_ids, input_type, comprehensive_mapping = NULL) {
    if (length(gene_ids) == 0) {
        return(list(
            entrez_ids = character(0),
            mapping_info = data.frame(),
            conversion_stats = list(
                input_count = 0,
                converted_count = 0,
                conversion_rate = 0,
                failed_ids = character(0)
            )
        ))
    }
    
    # Load mapping if not provided
    if (is.null(comprehensive_mapping)) {
        comprehensive_mapping <- tryCatch({
            load_comprehensive_gene_mapping()
        }, error = function(e) {
            warning("Could not load gene mapping data: ", e$message)
            return(NULL)
        })
    }
    
    if (is.null(comprehensive_mapping)) {
        return(list(
            entrez_ids = character(0),
            mapping_info = data.frame(),
            conversion_stats = list(
                input_count = length(gene_ids),
                converted_count = 0,
                conversion_rate = 0,
                failed_ids = gene_ids
            )
        ))
    }
    
    # Clean input IDs
    gene_ids <- unique(trimws(as.character(gene_ids)))
    gene_ids <- gene_ids[gene_ids != "" & !is.na(gene_ids)]
    
    # Convert based on input type
    result <- switch(input_type,
        "entrez" = convert_entrez_to_entrez(gene_ids, comprehensive_mapping),
        "symbol" = convert_symbol_to_entrez(gene_ids, comprehensive_mapping),
        "ensembl" = convert_ensembl_to_entrez(gene_ids, comprehensive_mapping),
        "uniprot" = convert_uniprot_to_entrez(gene_ids, comprehensive_mapping),
        {
            # Default case - try to auto-detect or treat as symbols
            warning("Unknown gene ID type: ", input_type, ". Treating as symbols.")
            convert_symbol_to_entrez(gene_ids, comprehensive_mapping)
        }
    )
    
    return(result)
}

#' Convert Entrez IDs to Entrez IDs (validation only)
convert_entrez_to_entrez <- function(entrez_ids, mapping_data) {
    # Validate that these are actual Entrez IDs in our mapping
    valid_entrez <- mapping_data$entrezgene_id[!is.na(mapping_data$entrezgene_id)]
    
    # Convert to character for consistent comparison
    entrez_ids <- as.character(entrez_ids)
    valid_entrez <- as.character(valid_entrez)
    
    # Find which ones are valid
    found_mask <- entrez_ids %in% valid_entrez
    found_ids <- entrez_ids[found_mask]
    failed_ids <- entrez_ids[!found_mask]
    
    # Create mapping info for found IDs
    if (length(found_ids) > 0) {
        mapping_subset <- mapping_data[!is.na(mapping_data$entrezgene_id) & 
                                     as.character(mapping_data$entrezgene_id) %in% found_ids, ]
        
        mapping_info <- data.frame(
            input_id = found_ids,
            entrez_id = found_ids,
            symbol = sapply(found_ids, function(id) {
                matches <- mapping_subset[as.character(mapping_subset$entrezgene_id) == id, ]
                if (nrow(matches) > 0 && !is.na(matches$hgnc_symbol[1])) {
                    return(matches$hgnc_symbol[1])
                } else {
                    return(id)  # Use Entrez ID as fallback
                }
            }),
            input_type = "entrez",
            stringsAsFactors = FALSE
        )
    } else {
        mapping_info <- data.frame(
            input_id = character(0),
            entrez_id = character(0),
            symbol = character(0),
            input_type = character(0),
            stringsAsFactors = FALSE
        )
    }
    
    return(list(
        entrez_ids = found_ids,
        mapping_info = mapping_info,
        conversion_stats = list(
            input_count = length(entrez_ids),
            converted_count = length(found_ids),
            conversion_rate = if (length(entrez_ids) > 0) length(found_ids) / length(entrez_ids) else 0,
            failed_ids = failed_ids
        )
    ))
}

#' Convert gene symbols to Entrez IDs
convert_symbol_to_entrez <- function(symbols, mapping_data) {
    # Case-insensitive matching for symbols
    symbols_upper <- toupper(symbols)
    mapping_upper <- toupper(mapping_data$hgnc_symbol)
    
    # Find matches
    matches <- mapping_data[!is.na(mapping_data$hgnc_symbol) & 
                           !is.na(mapping_data$entrezgene_id) &
                           mapping_upper %in% symbols_upper, ]
    
    # Create mapping info
    mapping_info <- data.frame(
        input_id = character(0),
        entrez_id = character(0),
        symbol = character(0),
        input_type = character(0),
        stringsAsFactors = FALSE
    )
    
    entrez_ids <- character(0)
    failed_ids <- symbols
    
    if (nrow(matches) > 0) {
        # Create mapping for each input symbol
        for (symbol in symbols) {
            symbol_upper <- toupper(symbol)
            symbol_matches <- matches[toupper(matches$hgnc_symbol) == symbol_upper, ]
            
            if (nrow(symbol_matches) > 0) {
                # Take the first match (remove duplicates)
                best_match <- symbol_matches[1, ]
                entrez_id <- as.character(best_match$entrezgene_id)
                
                mapping_info <- rbind(mapping_info, data.frame(
                    input_id = symbol,
                    entrez_id = entrez_id,
                    symbol = best_match$hgnc_symbol,
                    input_type = "symbol",
                    stringsAsFactors = FALSE
                ))
                
                entrez_ids <- c(entrez_ids, entrez_id)
                failed_ids <- setdiff(failed_ids, symbol)
            }
        }
        
        # Remove duplicates from entrez_ids
        entrez_ids <- unique(entrez_ids)
    }
    
    return(list(
        entrez_ids = entrez_ids,
        mapping_info = mapping_info,
        conversion_stats = list(
            input_count = length(symbols),
            converted_count = length(entrez_ids),
            conversion_rate = if (length(symbols) > 0) length(entrez_ids) / length(symbols) else 0,
            failed_ids = failed_ids
        )
    ))
}

#' Convert Ensembl IDs to Entrez IDs
convert_ensembl_to_entrez <- function(ensembl_ids, mapping_data) {
    matches <- mapping_data[!is.na(mapping_data$ensembl_gene_id) & 
                           !is.na(mapping_data$entrezgene_id) &
                           mapping_data$ensembl_gene_id %in% ensembl_ids, ]
    
    mapping_info <- data.frame(
        input_id = character(0),
        entrez_id = character(0),
        symbol = character(0),
        input_type = character(0),
        stringsAsFactors = FALSE
    )
    
    entrez_ids <- character(0)
    failed_ids <- ensembl_ids
    
    if (nrow(matches) > 0) {
        for (ensembl_id in ensembl_ids) {
            ensembl_matches <- matches[matches$ensembl_gene_id == ensembl_id, ]
            
            if (nrow(ensembl_matches) > 0) {
                best_match <- ensembl_matches[1, ]
                entrez_id <- as.character(best_match$entrezgene_id)
                
                mapping_info <- rbind(mapping_info, data.frame(
                    input_id = ensembl_id,
                    entrez_id = entrez_id,
                    symbol = ifelse(is.na(best_match$hgnc_symbol), entrez_id, best_match$hgnc_symbol),
                    input_type = "ensembl",
                    stringsAsFactors = FALSE
                ))
                
                entrez_ids <- c(entrez_ids, entrez_id)
                failed_ids <- setdiff(failed_ids, ensembl_id)
            }
        }
        
        entrez_ids <- unique(entrez_ids)
    }
    
    return(list(
        entrez_ids = entrez_ids,
        mapping_info = mapping_info,
        conversion_stats = list(
            input_count = length(ensembl_ids),
            converted_count = length(entrez_ids),
            conversion_rate = if (length(ensembl_ids) > 0) length(entrez_ids) / length(ensembl_ids) else 0,
            failed_ids = failed_ids
        )
    ))
}

#' Convert UniProt IDs to Entrez IDs
convert_uniprot_to_entrez <- function(uniprot_ids, mapping_data) {
    matches <- mapping_data[!is.na(mapping_data$uniprotswissprot) & 
                           !is.na(mapping_data$entrezgene_id) &
                           mapping_data$uniprotswissprot %in% uniprot_ids, ]
    
    mapping_info <- data.frame(
        input_id = character(0),
        entrez_id = character(0),
        symbol = character(0),
        input_type = character(0),
        stringsAsFactors = FALSE
    )
    
    entrez_ids <- character(0)
    failed_ids <- uniprot_ids
    
    if (nrow(matches) > 0) {
        for (uniprot_id in uniprot_ids) {
            uniprot_matches <- matches[matches$uniprotswissprot == uniprot_id, ]
            
            if (nrow(uniprot_matches) > 0) {
                best_match <- uniprot_matches[1, ]
                entrez_id <- as.character(best_match$entrezgene_id)
                
                mapping_info <- rbind(mapping_info, data.frame(
                    input_id = uniprot_id,
                    entrez_id = entrez_id,
                    symbol = ifelse(is.na(best_match$hgnc_symbol), entrez_id, best_match$hgnc_symbol),
                    input_type = "uniprot",
                    stringsAsFactors = FALSE
                ))
                
                entrez_ids <- c(entrez_ids, entrez_id)
                failed_ids <- setdiff(failed_ids, uniprot_id)
            }
        }
        
        entrez_ids <- unique(entrez_ids)
    }
    
    return(list(
        entrez_ids = entrez_ids,
        mapping_info = mapping_info,
        conversion_stats = list(
            input_count = length(uniprot_ids),
            converted_count = length(entrez_ids),
            conversion_rate = if (length(uniprot_ids) > 0) length(entrez_ids) / length(uniprot_ids) else 0,
            failed_ids = failed_ids
        )
    ))
}

#' Convert Entrez IDs back to symbols for display purposes
#' 
#' @param entrez_ids Vector of Entrez IDs (internal standard)
#' @param comprehensive_mapping Pre-loaded mapping data (optional)
#' @return Named vector: entrez_id -> symbol
entrez_to_symbols <- function(entrez_ids, comprehensive_mapping = NULL) {
    if (length(entrez_ids) == 0) {
        return(character(0))
    }
    
    if (is.null(comprehensive_mapping)) {
        comprehensive_mapping <- tryCatch({
            load_comprehensive_gene_mapping()
        }, error = function(e) {
            return(NULL)
        })
    }
    
    if (is.null(comprehensive_mapping)) {
        # Fallback: use Entrez IDs as symbols
        symbols <- setNames(as.character(entrez_ids), as.character(entrez_ids))
        return(symbols)
    }
    
    # Convert entrez_ids to character for matching
    entrez_ids <- as.character(entrez_ids)
    
    # Find symbols for these Entrez IDs
    symbols <- sapply(entrez_ids, function(entrez_id) {
        matches <- comprehensive_mapping[!is.na(comprehensive_mapping$entrezgene_id) & 
                                       as.character(comprehensive_mapping$entrezgene_id) == entrez_id, ]
        if (nrow(matches) > 0 && !is.na(matches$hgnc_symbol[1]) && matches$hgnc_symbol[1] != "") {
            return(matches$hgnc_symbol[1])
        } else {
            return(entrez_id)  # Fallback to Entrez ID
        }
    })
    
    return(symbols)
}

#' Get genes for a specific context (pathway, expression matching, etc.)
#' All internal operations use Entrez IDs, this function handles the conversion
#' 
#' @param source_type Type of gene source: "user", "pathway", "network_interactions"
#' @param context_data Additional data needed (pathway nodes, edges, etc.)
#' @param user_entrez_ids Current user's gene set (in Entrez format)
#' @param comprehensive_mapping Gene mapping data
#' @return List with entrez_ids and display_info
get_genes_for_context <- function(source_type, context_data = NULL, user_entrez_ids = character(0), comprehensive_mapping = NULL) {
    
    switch(source_type,
        "user" = {
            list(
                entrez_ids = user_entrez_ids,
                display_info = list(
                    message = paste("Using", length(user_entrez_ids), "genes from your uploaded gene set"),
                    source = "user_upload"
                )
            )
        },
        
        "pathway" = {
            if (is.null(context_data$nodes)) {
                return(list(
                    entrez_ids = character(0),
                    display_info = list(
                        message = "No pathway loaded",
                        source = "pathway"
                    )
                ))
            }
            
            # Extract genes from pathway nodes (KEGG nodes have hgnc_symbol and kegg_id)
            pathway_genes <- extract_pathway_entrez_ids(context_data$nodes, comprehensive_mapping)
            
            list(
                entrez_ids = pathway_genes,
                display_info = list(
                    message = paste("Using", length(pathway_genes), "genes from the selected pathway"),
                    source = "pathway"
                )
            )
        },
        
        "network_interactions" = {
            if (is.null(context_data$selected_gene) || is.null(context_data$nodes) || is.null(context_data$edges)) {
                return(list(
                    entrez_ids = character(0),
                    display_info = list(
                        message = "No network interactions available",
                        source = "network"
                    )
                ))
            }
            
            # Get interacting genes (this function should return Entrez IDs)
            interaction_result <- extract_network_interactions_entrez(
                context_data$selected_gene, 
                context_data$nodes, 
                context_data$edges,
                comprehensive_mapping
            )
            
            interaction_result
        },
        
        {
            # Default: return empty
            list(
                entrez_ids = character(0),
                display_info = list(
                    message = "Unknown gene source type",
                    source = "unknown"
                )
            )
        }
    )
}

#' Extract Entrez IDs from pathway nodes
extract_pathway_entrez_ids <- function(pathway_nodes, comprehensive_mapping) {
    if (is.null(pathway_nodes) || nrow(pathway_nodes) == 0) {
        return(character(0))
    }
    
    entrez_ids <- character(0)
    
    # Try to get Entrez IDs from different columns
    if ("kegg_id" %in% names(pathway_nodes)) {
        # KEGG IDs are often Entrez IDs with "hsa:" prefix
        kegg_ids <- pathway_nodes$kegg_id[!is.na(pathway_nodes$kegg_id)]
        # Remove "hsa:" prefix if present
        clean_kegg_ids <- gsub("^hsa:", "", kegg_ids)
        # Validate these are numeric (Entrez IDs)
        numeric_ids <- clean_kegg_ids[grepl("^[0-9]+$", clean_kegg_ids)]
        entrez_ids <- c(entrez_ids, numeric_ids)
    }
    
    # Also try converting HGNC symbols to Entrez if available
    if ("hgnc_symbol" %in% names(pathway_nodes) && !is.null(comprehensive_mapping)) {
        symbols <- pathway_nodes$hgnc_symbol[!is.na(pathway_nodes$hgnc_symbol) & pathway_nodes$hgnc_symbol != ""]
        if (length(symbols) > 0) {
            symbol_result <- convert_symbol_to_entrez(symbols, comprehensive_mapping)
            entrez_ids <- c(entrez_ids, symbol_result$entrez_ids)
        }
    }
    
    return(unique(entrez_ids))
}

#' Extract Entrez IDs for network interactions
extract_network_interactions_entrez <- function(selected_gene, nodes, edges, comprehensive_mapping) {
    # This would be implemented based on your existing network interaction logic
    # but always return Entrez IDs
    
    # For now, return empty - this should be implemented based on your existing logic
    return(list(
        entrez_ids = character(0),
        display_info = list(
            message = "Network interaction extraction not yet implemented in centralized system",
            source = "network"
        )
    ))
}
