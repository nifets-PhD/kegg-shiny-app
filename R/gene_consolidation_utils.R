#' Gene Consolidation Utilities for KEGG Pathway Visualization
#' 
#' This module handles the consolidation of KEGG gene entries and groups into
#' unified gene nodes for network visualization.

#' Consolidate KEGG genes and groups into unified gene nodes
#' @param kegg_coords data.frame with KEGG coordinate data including groups
#' @param comprehensive_mapping data.frame with gene ID mappings
#' @return list with consolidated nodes and excluded_ids
consolidate_kegg_genes_and_groups <- function(kegg_coords, comprehensive_mapping = NULL) {
    
    if (is.null(kegg_coords) || nrow(kegg_coords) == 0) {
        return(list(
            nodes = data.frame(),
            excluded_component_ids = character(0)
        ))
    }
    
    cat("Starting gene consolidation process...\n")
    cat("Input data:", nrow(kegg_coords), "entries\n")
    
    # Separate different types of entries
    gene_entries <- kegg_coords[kegg_coords$entry_type == "gene", ]
    group_entries <- kegg_coords[kegg_coords$entry_type == "group", ]
    other_entries <- kegg_coords[!kegg_coords$entry_type %in% c("gene", "group"), ]
    
    cat("Found:", nrow(gene_entries), "genes,", nrow(group_entries), "groups,", nrow(other_entries), "other entries\n")
    
    # Track which gene entries are components of groups (will be excluded)
    excluded_component_ids <- character(0)
    consolidated_nodes <- list()
    
    # Process groups first to identify component genes
    if (nrow(group_entries) > 0) {
        cat("Processing", nrow(group_entries), "group entries...\n")
        
        for (i in seq_len(nrow(group_entries))) {
            group_entry <- group_entries[i, ]
            group_result <- process_group_entry(group_entry, gene_entries, comprehensive_mapping)
            
            if (!is.null(group_result$consolidated_node)) {
                consolidated_nodes[[length(consolidated_nodes) + 1]] <- group_result$consolidated_node
                excluded_component_ids <- c(excluded_component_ids, group_result$component_ids)
                
                cat("  Group", group_entry$entry_id, "consolidated", length(group_result$component_ids), "components\n")
            }
        }
    }
    
    # Process individual gene entries (excluding those that are group components)
    remaining_genes <- gene_entries[!gene_entries$entry_id %in% excluded_component_ids, ]
    
    if (nrow(remaining_genes) > 0) {
        cat("Processing", nrow(remaining_genes), "individual gene entries...\n")
        
        for (i in seq_len(nrow(remaining_genes))) {
            gene_entry <- remaining_genes[i, ]
            gene_result <- process_gene_entry(gene_entry, comprehensive_mapping)
            
            if (!is.null(gene_result)) {
                consolidated_nodes[[length(consolidated_nodes) + 1]] <- gene_result
            }
        }
    }
    
    # Keep other entries (compounds, maps) as-is
    if (nrow(other_entries) > 0) {
        cat("Keeping", nrow(other_entries), "non-gene entries as-is\n")
        
        for (i in seq_len(nrow(other_entries))) {
            other_entry <- other_entries[i, ]
            other_result <- process_other_entry(other_entry)
            
            if (!is.null(other_result)) {
                consolidated_nodes[[length(consolidated_nodes) + 1]] <- other_result
            }
        }
    }
    
    # Combine all consolidated nodes
    if (length(consolidated_nodes) > 0) {
        final_nodes <- do.call(rbind, lapply(consolidated_nodes, function(x) {
            data.frame(
                id = x$id,
                label = x$label,
                type = x$type,
                gene_name = x$gene_name,
                hgnc_symbol = x$hgnc_symbol,
                kegg_id = x$kegg_id,
                description = x$description,
                x = x$x,
                y = x$y,
                kegg_width = x$kegg_width,
                kegg_height = x$kegg_height,
                kegg_shape = x$kegg_shape,
                kegg_bgcolor = x$kegg_bgcolor,
                kegg_fgcolor = x$kegg_fgcolor,
                kegg_type = x$kegg_type,
                kegg_label = x$kegg_label,
                title = x$title,
                # New fields for consolidated genes
                all_hsa_ids = x$all_hsa_ids,
                all_gene_symbols = x$all_gene_symbols,
                gene_count = x$gene_count,
                is_consolidated = x$is_consolidated,
                stringsAsFactors = FALSE
            )
        }))
    } else {
        final_nodes <- data.frame()
    }
    
    cat("Consolidation complete:", nrow(final_nodes), "final nodes\n")
    cat("Excluded", length(excluded_component_ids), "component gene entries\n")
    
    return(list(
        nodes = final_nodes,
        excluded_component_ids = excluded_component_ids
    ))
}

#' Process a group entry to create a consolidated gene node
#' @param group_entry single row data.frame with group information
#' @param gene_entries data.frame with all gene entries
#' @param comprehensive_mapping data.frame with gene ID mappings
#' @return list with consolidated_node and component_ids
process_group_entry <- function(group_entry, gene_entries, comprehensive_mapping) {
    
    # Extract component IDs from the group
    # Note: In the XML, groups have <component id="xxx"/> elements
    # For now, we'll parse this from the description or handle it in the calling function
    
    # This is a placeholder - in practice, you'd need to extract component IDs
    # from the XML parsing process where group components are identified
    component_ids <- extract_group_components(group_entry)
    
    if (length(component_ids) == 0) {
        cat("  Warning: Group", group_entry$entry_id, "has no components\n")
        return(list(consolidated_node = NULL, component_ids = character(0)))
    }
    
    # Get the component gene entries
    component_genes <- gene_entries[gene_entries$entry_id %in% component_ids, ]
    
    if (nrow(component_genes) == 0) {
        cat("  Warning: Group", group_entry$entry_id, "components not found in gene entries\n")
        return(list(consolidated_node = NULL, component_ids = character(0)))
    }
    
    # Extract all HSA IDs and gene symbols from components
    all_hsa_ids <- character(0)
    all_gene_symbols <- character(0)
    
    for (i in seq_len(nrow(component_genes))) {
        gene_data <- extract_gene_data_from_entry(component_genes[i, ], comprehensive_mapping)
        all_hsa_ids <- c(all_hsa_ids, gene_data$hsa_ids)
        all_gene_symbols <- c(all_gene_symbols, gene_data$symbols)
    }
    
    # Remove duplicates and clean up
    all_hsa_ids <- unique(all_hsa_ids[!is.na(all_hsa_ids) & all_hsa_ids != ""])
    all_gene_symbols <- unique(all_gene_symbols[!is.na(all_gene_symbols) & all_gene_symbols != ""])
    
    # For groups, show all gene symbols separated by newlines (as was done before)
    display_label <- paste(all_gene_symbols, collapse = "\n")
    
    # Create consolidated node
    consolidated_node <- list(
        id = paste0("kegg_", group_entry$entry_id),
        label = display_label,  # Use concise label for clean visualization
        type = "gene",  # Groups become gene nodes
        gene_name = paste("Group:", paste(component_ids, collapse = ",")),
        hgnc_symbol = paste(all_gene_symbols, collapse = "; "),
        kegg_id = if(length(all_hsa_ids) > 0) all_hsa_ids[1] else group_entry$entry_id,
        description = paste("Gene complex/group containing:", paste(all_gene_symbols, collapse = ", ")),
        x = group_entry$x,
        y = group_entry$y,
        kegg_width = group_entry$width,
        kegg_height = group_entry$height,
        kegg_shape = group_entry$shape,
        kegg_bgcolor = "#BFFFBF",  # Gene color instead of group color
        kegg_fgcolor = group_entry$fgcolor,
        kegg_type = "gene",
        kegg_label = paste(all_gene_symbols, collapse = ", "),
        title = create_group_tooltip(all_gene_symbols, all_hsa_ids),
        # Extended fields for consolidated genes
        all_hsa_ids = paste(all_hsa_ids, collapse = ","),
        all_gene_symbols = paste(all_gene_symbols, collapse = ","),
        gene_count = length(all_gene_symbols),
        is_consolidated = TRUE
    )
    
    return(list(
        consolidated_node = consolidated_node,
        component_ids = component_ids
    ))
}

#' Process an individual gene entry (potentially with multiple HSA IDs)
#' @param gene_entry single row data.frame with gene information
#' @param comprehensive_mapping data.frame with gene ID mappings
#' @return list representing the consolidated gene node
process_gene_entry <- function(gene_entry, comprehensive_mapping) {
    
    # Extract gene data (HSA IDs and symbols)
    gene_data <- extract_gene_data_from_entry(gene_entry, comprehensive_mapping)
    
    if (length(gene_data$hsa_ids) == 0) {
        cat("  Warning: Gene entry", gene_entry$entry_id, "has no valid HSA IDs\n")
        return(NULL)
    }
    
    # Determine if this is a multi-gene entry
    is_multi_gene <- length(gene_data$hsa_ids) > 1
    
    # Always use the original XML label for display (keep visualization clean)
    display_label <- gene_entry$label
    
    # Create consolidated node
    consolidated_node <- list(
        id = paste0("kegg_", gene_entry$entry_id),
        label = display_label,  # Use original XML label for clean visualization
        type = "gene",
        gene_name = gene_entry$kegg_name,
        hgnc_symbol = paste(gene_data$symbols, collapse = "; "),
        kegg_id = gene_data$hsa_ids[1],  # Primary HSA ID
        description = if (is_multi_gene) {
            paste("Gene complex/group containing:", paste(gene_data$symbols, collapse = ", "))
        } else {
            paste("Gene:", gene_data$symbols[1])
        },
        x = gene_entry$x,
        y = gene_entry$y,
        kegg_width = gene_entry$width,
        kegg_height = gene_entry$height,
        kegg_shape = gene_entry$shape,
        kegg_bgcolor = gene_entry$bgcolor,
        kegg_fgcolor = gene_entry$fgcolor,
        kegg_type = gene_entry$entry_type,
        kegg_label = gene_entry$label,
        title = create_gene_tooltip(gene_data$symbols, gene_data$hsa_ids),
        # Extended fields for consolidated genes
        all_hsa_ids = paste(gene_data$hsa_ids, collapse = ","),
        all_gene_symbols = paste(gene_data$symbols, collapse = ","),
        gene_count = length(gene_data$symbols),
        is_consolidated = is_multi_gene
    )
    
    return(consolidated_node)
}

#' Process non-gene entries (compounds, maps) 
#' @param other_entry single row data.frame with entry information
#' @return list representing the node
process_other_entry <- function(other_entry) {
    
    # Handle compound entries specially to get proper names
    if (other_entry$entry_type == "compound") {
        # Split multiple compound IDs if present (e.g., "cpd:C00668 cpd:C01172")
        compound_id_string <- other_entry$kegg_name
        compound_ids <- trimws(strsplit(compound_id_string, "\\s+")[[1]])
        
        # Use the existing get_compound_info function for proper name lookup
        compound_info <- get_compound_info(compound_ids)
        
        # Use the compound name for label if available, otherwise use the original label
        display_label <- if (!is.na(compound_info$first_name)) {
            compound_info$first_name
        } else {
            other_entry$label
        }
        
        # Create enhanced description for compounds
        description <- if (length(compound_info$all_names) > 0) {
            paste("Compound:", paste(compound_info$all_names, collapse = "; "))
        } else {
            paste("Compound:", other_entry$label)
        }
        
        consolidated_node <- list(
            id = paste0("kegg_", other_entry$entry_id),
            label = display_label,
            type = other_entry$entry_type,
            gene_name = other_entry$kegg_name,
            hgnc_symbol = display_label,  # Use compound name for consistency
            kegg_id = other_entry$entry_id,
            description = description,
            x = other_entry$x,
            y = other_entry$y,
            kegg_width = other_entry$width,
            kegg_height = other_entry$height,
            kegg_shape = other_entry$shape,
            kegg_bgcolor = other_entry$bgcolor,
            kegg_fgcolor = other_entry$fgcolor,
            kegg_type = other_entry$entry_type,
            kegg_label = other_entry$label,
            title = compound_info$tooltip,
            # Extended fields (empty for non-genes)
            all_hsa_ids = "",
            all_gene_symbols = "",
            gene_count = 0,
            is_consolidated = FALSE
        )
        
    } else {
        # Handle other entry types (maps, etc.) as before
        consolidated_node <- list(
            id = paste0("kegg_", other_entry$entry_id),
            label = other_entry$label,
            type = other_entry$entry_type,
            gene_name = other_entry$kegg_name,
            hgnc_symbol = other_entry$label,
            kegg_id = other_entry$entry_id,
            description = other_entry$description,
            x = other_entry$x,
            y = other_entry$y,
            kegg_width = other_entry$width,
            kegg_height = other_entry$height,
            kegg_shape = other_entry$shape,
            kegg_bgcolor = other_entry$bgcolor,
            kegg_fgcolor = other_entry$fgcolor,
            kegg_type = other_entry$entry_type,
            kegg_label = other_entry$label,
            title = paste0(other_entry$entry_type, ": ", other_entry$label),
            # Extended fields (empty for non-genes)
            all_hsa_ids = "",
            all_gene_symbols = "",
            gene_count = 0,
            is_consolidated = FALSE
        )
    }
    
    return(consolidated_node)
}

#' Extract HSA IDs and gene symbols from a gene entry
#' @param gene_entry single row data.frame with gene information
#' @param comprehensive_mapping data.frame with gene ID mappings
#' @return list with hsa_ids and symbols vectors
extract_gene_data_from_entry <- function(gene_entry, comprehensive_mapping) {
    
    # Extract HSA IDs from the KEGG name field (e.g., "hsa:10000 hsa:207 hsa:208")
    hsa_ids <- regmatches(gene_entry$kegg_name, gregexpr("\\d+", gene_entry$kegg_name))[[1]]
    
    # Extract gene symbols from the graphics name field (e.g., "AKT3, MPPH, MPPH2, PKB-GAMMA...")
    symbols_text <- gene_entry$label
    
    if (!is.na(symbols_text) && symbols_text != "") {
        # Split by comma and take unique symbols
        symbols <- unique(trimws(strsplit(symbols_text, ",")[[1]]))
        # Remove ellipsis and clean up
        symbols <- gsub("\\.\\.\\.$", "", symbols)
        symbols <- symbols[symbols != "" & !is.na(symbols)]
    } else {
        symbols <- character(0)
    }
    
    # If we have comprehensive mapping, try to get better symbols for the HSA IDs
    if (!is.null(comprehensive_mapping) && length(hsa_ids) > 0) {
        improved_symbols <- get_symbols_from_entrez(hsa_ids, comprehensive_mapping)
        if (length(improved_symbols) > 0) {
            symbols <- improved_symbols
        }
    }
    
    # Ensure we have at least some symbols
    if (length(symbols) == 0 && length(hsa_ids) > 0) {
        symbols <- paste0("Gene", hsa_ids)
    }
    
    return(list(
        hsa_ids = hsa_ids,
        symbols = symbols
    ))
}

#' Get gene symbols from Entrez IDs using comprehensive mapping
#' @param entrez_ids character vector of Entrez IDs
#' @param comprehensive_mapping data.frame with gene ID mappings
#' @return character vector of gene symbols
get_symbols_from_entrez <- function(entrez_ids, comprehensive_mapping) {
    
    symbols <- character(0)
    
    for (entrez_id in entrez_ids) {
        matches <- comprehensive_mapping[
            !is.na(comprehensive_mapping$entrezgene_id) & 
            as.character(comprehensive_mapping$entrezgene_id) == entrez_id, 
        ]
        
        if (nrow(matches) > 0) {
            symbol <- matches$hgnc_symbol[1]
            if (!is.na(symbol) && symbol != "") {
                symbols <- c(symbols, symbol)
            }
        }
    }
    
    return(unique(symbols))
}

#' Create tooltip for gene nodes
#' @param symbols character vector of gene symbols
#' @param hsa_ids character vector of HSA IDs
#' @return character string for tooltip
create_gene_tooltip <- function(symbols, hsa_ids) {
    
    if (length(symbols) == 1) {
        return(paste0("Gene: ", symbols[1], "\nHSA ID: ", hsa_ids[1]))
    } else {
        return(paste0("Multi-gene entry (", length(symbols), " genes)\n",
                     paste(paste0(symbols, " (", hsa_ids, ")"), collapse = "\n")))
    }
}

#' Create tooltip for group nodes
#' @param symbols character vector of gene symbols
#' @param hsa_ids character vector of HSA IDs
#' @return character string for tooltip
create_group_tooltip <- function(symbols, hsa_ids) {
    
    return(paste0("Gene complex/group (", length(symbols), " genes)\n",
                 paste(paste0(symbols, " (", hsa_ids, ")"), collapse = "\n")))
}

#' Extract component IDs from a group entry
#' @param group_entry single row data.frame with group information
#' @return character vector of component entry IDs
extract_group_components <- function(group_entry) {
    
    # Extract component IDs from the group_components field
    if (!is.null(group_entry$group_components) && !is.na(group_entry$group_components) && 
        group_entry$group_components != "") {
        component_ids <- strsplit(group_entry$group_components, ",")[[1]]
        component_ids <- trimws(component_ids)
        component_ids <- component_ids[component_ids != "" & !is.na(component_ids)]
        return(component_ids)
    }
    
    return(character(0))
}

#' Enhanced gene highlighting for consolidated nodes
#' @param nodes data.frame with consolidated node information
#' @param highlight_genes character vector of gene identifiers to highlight
#' @return nodes data.frame with highlighting applied
apply_consolidated_gene_highlighting <- function(nodes, highlight_genes) {
    
    if (is.null(nodes) || nrow(nodes) == 0 || length(highlight_genes) == 0) {
        return(nodes)
    }
    
    cat("Applying gene highlighting to consolidated nodes...\n")
    
    # Convert highlight genes to character and uppercase for matching
    highlight_genes <- toupper(as.character(highlight_genes))
    
    # Initialize highlighting column
    nodes$is_highlighted <- FALSE
    
    # Only consider gene nodes
    gene_nodes <- nodes$type == "gene" & !is.na(nodes$type)
    
    if (sum(gene_nodes) == 0) {
        return(nodes)
    }
    
    # Check each gene node for matches
    for (i in which(gene_nodes)) {
        node_matched <- FALSE
        
        # Check HSA IDs (Entrez IDs)
        if (!is.na(nodes$all_hsa_ids[i]) && nodes$all_hsa_ids[i] != "") {
            node_hsa_ids <- strsplit(nodes$all_hsa_ids[i], ",")[[1]]
            numeric_highlights <- highlight_genes[grepl("^\\d+$", highlight_genes)]
            
            if (any(numeric_highlights %in% node_hsa_ids)) {
                node_matched <- TRUE
                cat("  Matched node", nodes$id[i], "by HSA ID\n")
            }
        }
        
        # Check gene symbols
        if (!node_matched && !is.na(nodes$all_gene_symbols[i]) && nodes$all_gene_symbols[i] != "") {
            node_symbols <- toupper(strsplit(nodes$all_gene_symbols[i], ",")[[1]])
            
            if (any(highlight_genes %in% node_symbols)) {
                node_matched <- TRUE
                cat("  Matched node", nodes$id[i], "by gene symbol\n")
            }
        }
        
        # Check main HGNC symbol field
        if (!node_matched && !is.na(nodes$hgnc_symbol[i]) && nodes$hgnc_symbol[i] != "") {
            node_symbols <- toupper(strsplit(nodes$hgnc_symbol[i], ";")[[1]])
            node_symbols <- trimws(node_symbols)
            
            if (any(highlight_genes %in% node_symbols)) {
                node_matched <- TRUE
                cat("  Matched node", nodes$id[i], "by HGNC symbol\n")
            }
        }
        
        nodes$is_highlighted[i] <- node_matched
    }
    
    highlighted_count <- sum(nodes$is_highlighted)
    cat("Highlighted", highlighted_count, "out of", sum(gene_nodes), "gene nodes\n")
    
    return(nodes)
}

#' Check if a gene list contains any genes from a consolidated node
#' @param gene_list character vector of gene identifiers
#' @param consolidated_node list representing a consolidated gene node
#' @return logical indicating if there's a match
gene_list_matches_node <- function(gene_list, consolidated_node) {
    
    if (length(gene_list) == 0 || consolidated_node$type != "gene") {
        return(FALSE)
    }
    
    # Convert to uppercase for matching
    gene_list_upper <- toupper(as.character(gene_list))
    
    # Check HSA IDs
    if (!is.na(consolidated_node$all_hsa_ids) && consolidated_node$all_hsa_ids != "") {
        node_hsa_ids <- strsplit(consolidated_node$all_hsa_ids, ",")[[1]]
        numeric_genes <- gene_list[grepl("^\\d+$", gene_list)]
        
        if (any(numeric_genes %in% node_hsa_ids)) {
            return(TRUE)
        }
    }
    
    # Check gene symbols
    if (!is.na(consolidated_node$all_gene_symbols) && consolidated_node$all_gene_symbols != "") {
        node_symbols <- toupper(strsplit(consolidated_node$all_gene_symbols, ",")[[1]])
        
        if (any(gene_list_upper %in% node_symbols)) {
            return(TRUE)
        }
    }
    
    return(FALSE)
}
