# Helper function: Get pathway genes with UniProt IDs
get_pathway_genes_with_uniprot <- function(current_node_id = NULL, nodes_data, comprehensive_mapping) {
    pathway_genes_with_uniprot <- c()
    if (!is.null(nodes_data) && nrow(nodes_data) > 0 && !is.null(comprehensive_mapping)) {
        for (i in seq_len(nrow(nodes_data))) {
            node <- nodes_data[i, ]
            if (!is.null(node$hgnc_symbol) && !is.na(node$hgnc_symbol) && node$hgnc_symbol != "" &&
                (is.null(current_node_id) || node$id != current_node_id)) {  # Exclude current node if specified
                # Look up UniProt for this node
                node_mapping <- comprehensive_mapping[
                    !is.na(comprehensive_mapping$hgnc_symbol) & 
                    comprehensive_mapping$hgnc_symbol == node$hgnc_symbol & 
                    !is.na(comprehensive_mapping$uniprotswissprot) & 
                    comprehensive_mapping$uniprotswissprot != "", 
                ]
                if (nrow(node_mapping) > 0) {
                    pathway_genes_with_uniprot <- c(pathway_genes_with_uniprot, node$hgnc_symbol)
                }
            }
        }
    }
    return(pathway_genes_with_uniprot)
}

# Helper function: Get gene info for consolidated nodes  
get_gene_info_for_consolidated <- function(gene_symbol, hsa_id = NULL, comprehensive_mapping = NULL, pathway_genes_with_uniprot = NULL, selected_node = NULL) {
    # Initialize output
    html_output <- ""
    uniprot_html <- ""
    
    # Clean inputs
    gene_symbol <- trimws(gene_symbol)
    if (!is.null(hsa_id)) hsa_id <- trimws(hsa_id)
    
    # Extract Entrez ID from HSA format
    entrez_id <- if (!is.null(hsa_id) && grepl("^hsa:", hsa_id)) {
        sub("^hsa:", "", hsa_id)
    } else {
        hsa_id
    }
    
    # Get comprehensive gene information from mapping
    gene_info <- NULL
    if (!is.null(comprehensive_mapping) && !is.null(gene_symbol) && gene_symbol != "") {
        tryCatch({
            # Look up by HGNC symbol first
            gene_matches <- comprehensive_mapping[
                !is.na(comprehensive_mapping$hgnc_symbol) & 
                comprehensive_mapping$hgnc_symbol == gene_symbol, 
            ]
            
            # If no symbol match, try by Entrez ID
            if (nrow(gene_matches) == 0 && !is.null(entrez_id) && entrez_id != "") {
                gene_matches <- comprehensive_mapping[
                    !is.na(comprehensive_mapping$entrezgene_id) & 
                    comprehensive_mapping$entrezgene_id == entrez_id, 
                ]
            }
            
            if (nrow(gene_matches) > 0) {
                gene_info <- gene_matches[1, ]  # Take first match if multiple
            }
        }, error = function(e) {
            cat("Error getting gene info for", gene_symbol, ":", e$message, "\n")
        })
    }
    
    # Build basic gene info
    if (!is.null(gene_info)) {
        html_output <- paste0(
            "<strong>Gene Symbol:</strong> ", gene_symbol, "<br>",
            if (!is.null(gene_info$gene_name) && !is.na(gene_info$gene_name)) 
                paste0("<strong>Gene Name:</strong> ", gene_info$gene_name, "<br>") else "",
            if (!is.null(entrez_id) && entrez_id != "") 
                paste0("<strong>Entrez ID:</strong> ", entrez_id, "<br>") else "",
            if (!is.null(gene_info$ensembl_gene_id) && !is.na(gene_info$ensembl_gene_id)) 
                paste0("<strong>Ensembl ID:</strong> ", gene_info$ensembl_gene_id, "<br>") else "",
            if (!is.null(gene_info$uniprotswissprot) && !is.na(gene_info$uniprotswissprot)) 
                paste0("<strong>UniProt ID:</strong> ", gene_info$uniprotswissprot, "<br>") else ""
        )
        
        # Build UniProt 3D structure links
        if (!is.null(gene_info$uniprotswissprot) && !is.na(gene_info$uniprotswissprot) && gene_info$uniprotswissprot != "") {
            uniprot_html <- paste0(
                "<a href='https://www.uniprot.org/uniprot/", gene_info$uniprotswissprot, "' target='_blank' style='color: #007bff; text-decoration: underline;'>",
                "View ", gene_info$uniprotswissprot, " in UniProt ↗</a><br>",
                "<a href='https://alphafold.ebi.ac.uk/entry/", gene_info$uniprotswissprot, "' target='_blank' style='color: #007bff; text-decoration: underline;'>",
                "View AlphaFold Structure ↗</a><br>"
            )
        }
    } else {
        html_output <- paste0(
            "<strong>Gene Symbol:</strong> ", gene_symbol, "<br>",
            if (!is.null(entrez_id) && entrez_id != "") 
                paste0("<strong>Entrez ID:</strong> ", entrez_id, "<br>") else "",
            "<em>Additional gene information not available</em><br>"
        )
    }
    
    return(list(
        html = html_output,
        uniprot_html = uniprot_html,
        gene_info = gene_info
    ))
}

# Helper function: Extract gene symbol from node (only for gene nodes)
extract_gene_symbol <- function(selected_node, gene_name, comprehensive_mapping = NULL) {
    if ("title" %in% names(selected_node) && !is.null(selected_node$title) && 
        !is.na(selected_node$title) && selected_node$title != "") {
        title_text <- trimws(selected_node$title)
        if (grepl("^Gene: ", title_text)) {
            gene_part <- strsplit(title_text, "\n")[[1]][1]
            return(gsub("^Gene: ", "", gene_part))
        } else {
            return(title_text)
        }
    } else if (!is.null(gene_name) && !is.na(gene_name) && grepl("^hsa:", gene_name)) {
        entrez_id <- sub("^hsa:", "", gene_name)
        if (!is.null(comprehensive_mapping) && entrez_id %in% comprehensive_mapping$entrezgene_id) {
            match_idx <- which(comprehensive_mapping$entrezgene_id == entrez_id)[1]
            return(comprehensive_mapping$hgnc_symbol[match_idx])
        } else {
            return(entrez_id)
        }
    } else {
        return(gene_name)
    }
}

# Helper function: Build relationship info (common for all node types)
build_relationship_info <- function(incoming_edges, outgoing_edges, nodes_data) {
    if (nrow(incoming_edges) == 0 && nrow(outgoing_edges) == 0) {
        return("")
    }
    
    info <- paste0(
        "<br><strong>=== NETWORK RELATIONSHIPS ===</strong><br>",
        "Incoming connections: ", nrow(incoming_edges), "<br>",
        "Outgoing connections: ", nrow(outgoing_edges), "<br>"
    )
    
    # Show sample incoming relationships
    if (nrow(incoming_edges) > 0) {
        sample_incoming <- head(incoming_edges, 3)
        for (i in seq_len(nrow(sample_incoming))) {
            edge <- sample_incoming[i, ]
            rel_type <- if (!is.null(edge$relation_type) && !is.na(edge$relation_type)) edge$relation_type else "Unknown"
            subtype <- if (!is.null(edge$subtype) && !is.na(edge$subtype) && edge$subtype != "") edge$subtype else ""
            
            source_node <- nodes_data[nodes_data$id == edge$from, ]
            source_name <- edge$from
            if (nrow(source_node) > 0) {
                if (!is.null(source_node$hgnc_symbol) && !is.na(source_node$hgnc_symbol) && source_node$hgnc_symbol != "") {
                    source_name <- source_node$hgnc_symbol
                } else if (!is.null(source_node$label) && !is.na(source_node$label) && source_node$label != "") {
                    source_name <- source_node$label
                }
            }
            
            info <- paste0(info, "← ", source_name, " (", rel_type, 
                          if (subtype != "") paste0(": ", subtype) else "", ")", "<br>")
        }
    }
    
    # Show sample outgoing relationships
    if (nrow(outgoing_edges) > 0) {
        sample_outgoing <- head(outgoing_edges, 3)
        for (i in seq_len(nrow(sample_outgoing))) {
            edge <- sample_outgoing[i, ]
            rel_type <- if (!is.null(edge$relation_type) && !is.na(edge$relation_type)) edge$relation_type else "Unknown"
            subtype <- if (!is.null(edge$subtype) && !is.na(edge$subtype) && edge$subtype != "") edge$subtype else ""
            
            target_node <- nodes_data[nodes_data$id == edge$to, ]
            target_name <- edge$to
            if (nrow(target_node) > 0) {
                if (!is.null(target_node$hgnc_symbol) && !is.na(target_node$hgnc_symbol) && target_node$hgnc_symbol != "") {
                    target_name <- target_node$hgnc_symbol
                } else if (!is.null(target_node$label) && !is.na(target_node$label) && target_node$label != "") {
                    target_name <- target_node$label
                }
            }
            
            info <- paste0(info, "→ ", target_name, " (", rel_type, 
                          if (subtype != "") paste0(": ", subtype) else "", ")", "<br>")
        }
    }
    
    return(info)
}

# Helper function: Build phylostratum info (only for genes)
build_phylostratum_info <- function(gene_symbol, gene_name, selected_node) {
    phylo_stratum <- NULL
    
    # Try to get from node first
    if (!is.null(selected_node$phylostratum) && !is.na(selected_node$phylostratum)) {
        phylo_stratum <- selected_node$phylostratum
    } else {
        # Fallback: lookup phylostratum
        phylo_stratum <- tryCatch({
            if (!exists("phylomap_data") || is.null(phylomap_data)) {
                source("R/kegg_utils.R", local = TRUE)
                phylomap_data <<- load_phylomap(verbose = FALSE)
            }
            
            if (!is.null(phylomap_data)) {
                # Try by gene symbol first
                if (!is.null(gene_symbol) && !is.na(gene_symbol) && gene_symbol != "") {
                    gene_strata <- map_genes_to_phylostrata(gene_symbol, "symbol", phylomap = phylomap_data, verbose = FALSE)
                    if (!is.na(gene_strata[1])) {
                        cat("Found phylostratum for", gene_symbol, ":", gene_strata[1], "\n")
                        return(gene_strata[1])
                    }
                }
                
                # Try by Entrez ID
                if (!is.null(gene_name) && grepl("^hsa:", gene_name)) {
                    entrez_id <- sub("^hsa:", "", gene_name)
                    gene_strata <- map_genes_to_phylostrata(entrez_id, "entrez", phylomap = phylomap_data, verbose = FALSE)
                    if (!is.na(gene_strata[1])) {
                        cat("Found phylostratum for Entrez", entrez_id, ":", gene_strata[1], "\n")
                        return(gene_strata[1])
                    }
                }
            }
            return(NULL)
        }, error = function(e) NULL)
    }
    
    if (is.null(phylo_stratum) || is.na(phylo_stratum)) {
        return("")
    }
    
    # Get phylostratum display with color
    legend_data <- tryCatch({
        if (!exists("phylo_legend_data") || is.null(phylo_legend_data)) {
            source("R/kegg_utils.R", local = TRUE)
            phylo_legend_data <<- load_phylostratum_legend()
        }
        phylo_legend_data
    }, error = function(e) NULL)
    
    if (!is.null(legend_data)) {
        phylo_row <- legend_data[legend_data$Rank == phylo_stratum, ]
        if (nrow(phylo_row) > 0) {
            phylo_name <- phylo_row$Name[1]
            phylo_color <- if (!is.null(phylo_row$Color) && !is.na(phylo_row$Color[1])) phylo_row$Color[1] else "#000000"
            return(paste0(
                "<br><strong>=== PHYLOSTRATUM ===</strong><br>",
                "<div style='background-color: ", phylo_color, "; color: white; padding: 4px 8px; border-radius: 4px; display: inline-block; margin: 2px 0;'>",
                "<strong>Phylostratum:</strong> ", phylo_stratum, " - ", phylo_name,
                "</div><br>"
            ))
        }
    }
    
    return(paste0(
        "<br><strong>=== PHYLOSTRATUM ===</strong><br>",
        "<strong>Phylostratum:</strong> ", phylo_stratum, "<br>"
    ))
}

# GENE NODE HANDLER (unified for single and multi-gene)
build_gene_node_info <- function(selected_node, incoming_edges, outgoing_edges, nodes_data, comprehensive_mapping) {
    # Determine if consolidated or single gene
    is_consolidated <- !is.null(selected_node$is_consolidated) && 
                     !is.na(selected_node$is_consolidated) && 
                     selected_node$is_consolidated
    
    if (is_consolidated) {
        # Multi-gene node
        gene_count <- selected_node$gene_count
        all_hsa_ids <- strsplit(selected_node$all_hsa_ids, ",")[[1]]
        all_gene_symbols <- strsplit(selected_node$all_gene_symbols, ",")[[1]]
        node_label <- "CONSOLIDATED GENE NODE"
    } else {
        # Single gene node - treat as consolidated with count=1
        gene_count <- 1
        gene_symbol <- extract_gene_symbol(selected_node, selected_node$gene_name, comprehensive_mapping)
        all_hsa_ids <- c(selected_node$gene_name)
        all_gene_symbols <- c(gene_symbol)
        node_label <- "GENE NODE"
    }
    
    # Get pathway genes for EMERALD (common logic)
    pathway_genes_with_uniprot <- get_pathway_genes_with_uniprot(selected_node$id, nodes_data, comprehensive_mapping)
    
    # Process each gene (unified logic for single and multi)
    genes_info_html <- ""
    for (i in seq_along(all_gene_symbols)) {
        symbol <- trimws(all_gene_symbols[i])
        hsa_id <- if (i <= length(all_hsa_ids)) trimws(all_hsa_ids[i]) else ""
        
        gene_info_result <- get_gene_info_for_consolidated(symbol, hsa_id, comprehensive_mapping, pathway_genes_with_uniprot, selected_node)
        
        # Add phylostratum info for first gene only (to avoid repetition)
        phylo_info <- if (i == 1) build_phylostratum_info(symbol, hsa_id, selected_node) else ""
        
        genes_info_html <- paste0(genes_info_html, 
            "<div style='margin-bottom: 15px; border-left: 3px solid #007cba; padding-left: 10px; background-color: #f8f9fa;'>",
            if (gene_count > 1) paste0("<strong>Gene ", i, ":</strong> ", symbol, "<br>") else "",
            gene_info_result$html,
            if (gene_info_result$uniprot_html != "") paste0("<br><strong>=== 3D PROTEIN STRUCTURES ===</strong><br>", gene_info_result$uniprot_html) else "",
            phylo_info,
            "</div>"
        )
    }
    
    # Build relationship info
    relationship_info <- build_relationship_info(incoming_edges, outgoing_edges, nodes_data)
    
    # Create final HTML
    html_content <- paste0(
        "<div style='font-family: monospace; font-size: 14px; line-height: 1.6;'>",
        "<strong>=== ", node_label, " ===</strong><br>",
        "<strong>Node Type:</strong> ", if (gene_count > 1) "Multi-gene entry" else "Single gene", "<br>",
        "<strong>Gene Count:</strong> ", gene_count, if (gene_count > 1) " genes" else " gene", "<br>",
        "<strong>Display Label:</strong> ", selected_node$label, "<br>",
        "<strong>Gene Symbol(s):</strong> ", paste(all_gene_symbols, collapse=", "), "<br>",
        "<strong>KEGG Entry(s):</strong> <code>", paste(all_hsa_ids, collapse=", "), "</code><br>",
        
        "<br><strong>=== GENE INFORMATION ===</strong><br>",
        genes_info_html,
        
        relationship_info,
        
        "<br><strong>=== PATHWAY CONTEXT ===</strong><br>",
        "In KEGG, ", if (gene_count > 1) "these genes are" else "this gene is", " referenced as: <code>", paste(all_hsa_ids, collapse=", "), "</code><br>",
        "Common name/symbol: <code>", paste(all_gene_symbols, collapse=", "), "</code><br>",
        "<br>",
        "<em>Tip: Search for '", paste(all_hsa_ids, collapse="' or '"), "' in KEGG database<br>",
        "or '", paste(all_gene_symbols, collapse="' or '"), "' in gene databases</em>",
        
        # JavaScript for EMERALD (simplified, single function for both cases)
        "<script>
        function selectSecondGeneForGene(firstGene, firstUniProt, secondGeneSymbol, containerId) {
            if (!secondGeneSymbol) {
                document.getElementById(containerId + '_container').innerHTML = '';
                return;
            }
            Shiny.setInputValue('emerald_request', {
                firstGene: firstGene,
                firstUniProt: firstUniProt,
                secondGene: secondGeneSymbol,
                containerId: containerId
            }, {priority: 'event'});
        }
        </script>",
        
        "</div>"
    )
    
    return(html_content)
}

# COMPOUND NODE HANDLER
build_compound_node_info <- function(selected_node, incoming_edges, outgoing_edges, nodes_data) {
    # Extract compound-specific information
    compound_names <- character(0)
    if (!is.null(selected_node$description) && !is.na(selected_node$description)) {
        desc_parts <- strsplit(selected_node$description, "Compound: ")[[1]]
        if (length(desc_parts) > 1) {
            compound_names <- trimws(strsplit(desc_parts[2], ";")[[1]])
        }
    }
    
    # Extract compound IDs
    compound_ids <- character(0)
    if (!is.null(selected_node$gene_name) && !is.na(selected_node$gene_name)) {
        compound_ids <- regmatches(selected_node$gene_name, gregexpr("C\\d{5}", selected_node$gene_name))[[1]]
    }
    
    # Build compound info section
    compound_info <- ""
    if (length(compound_names) > 0) {
        compound_info <- paste0(
            "<strong>Compound Names:</strong><br>",
            paste("• ", compound_names, collapse = "<br>"), "<br><br>"
        )
    }
    
    # Build KEGG links
    compound_links <- ""
    if (length(compound_ids) > 0) {
        compound_links <- paste0(
            "<br><strong>=== KEGG DATABASE LINKS ===</strong><br>",
            paste(sapply(compound_ids, function(cid) {
                paste0("<strong>", cid, ":</strong> ",
                      "<a href='https://www.kegg.jp/entry/", cid, "' target='_blank' style='color: #007bff; text-decoration: underline;'>",
                      "View in KEGG ↗</a>")
            }), collapse = "<br>"), "<br>"
        )
    } else if (!is.null(selected_node$kegg_id) && selected_node$kegg_id != "") {
        compound_links <- paste0(
            "<br><strong>=== KEGG DATABASE LINK ===</strong><br>",
            "<a href='https://www.kegg.jp/entry/", selected_node$kegg_id, "' target='_blank' style='color: #007bff; text-decoration: underline;'>",
            "View ", selected_node$kegg_id, " in KEGG ↗</a><br>"
        )
    }
    
    relationship_info <- build_relationship_info(incoming_edges, outgoing_edges, nodes_data)
    
    return(paste0(
        "<div style='font-family: monospace; font-size: 14px; line-height: 1.6;'>",
        "<strong>=== COMPOUND INFORMATION ===</strong><br>",
        compound_info,
        "<strong>KEGG Compound ID:</strong> <code>", selected_node$kegg_id, "</code><br>",
        "<strong>Original KEGG Entry:</strong> <code>", selected_node$gene_name, "</code><br>",
        "<strong>Display Label:</strong> ", selected_node$label, "<br>",
        relationship_info,
        compound_links,
        "<br><strong>=== PATHWAY CONTEXT ===</strong><br>",
        "This compound participates in pathway reactions.<br>",
        "<br>",
        "<em>Tip: Click the KEGG links above to view detailed compound information</em>",
        "</div>"
    ))
}

# MAP NODE HANDLER
build_map_node_info <- function(selected_node, incoming_edges, outgoing_edges, nodes_data) {
    # Extract pathway ID
    pathway_id <- ""
    if (!is.null(selected_node$gene_name) && selected_node$gene_name != "") {
        pathway_match <- regmatches(selected_node$gene_name, regexpr("hsa\\d+", selected_node$gene_name))
        if (length(pathway_match) > 0) {
            pathway_id <- pathway_match[1]
        }
    }
    
    # Build pathway links
    pathway_links <- ""
    if (pathway_id != "") {
        pathway_links <- paste0(
            "<br><strong>=== PATHWAY NAVIGATION ===</strong><br>",
            "<strong>KEGG Pathway Link:</strong> ",
            "<a href='https://www.kegg.jp/kegg-bin/show_pathway?", pathway_id, "' target='_blank' style='color: #007bff; text-decoration: underline;'>",
            "View in KEGG ↗</a><br>",
            "<strong>Load in Visualizer:</strong> ",
            "<button id='load_pathway_btn_", selected_node$id, "' onclick='loadPathwayFromNode(\"", pathway_id, "\")' ",
            "style='background-color: #007bff; color: white; border: none; padding: 5px 10px; border-radius: 3px; cursor: pointer; margin-top: 5px;'>",
            "Load ", selected_node$label, "</button><br>"
        )
    }
    
    relationship_info <- build_relationship_info(incoming_edges, outgoing_edges, nodes_data)
    
    return(paste0(
        "<div style='font-family: monospace; font-size: 14px; line-height: 1.6;'>",
        "<strong>=== PATHWAY MAP INFORMATION ===</strong><br>",
        if (pathway_id != "") paste0("<strong>Pathway ID:</strong> <code>", pathway_id, "</code><br>") else "",
        "<strong>Original KEGG Entry:</strong> <code>", selected_node$gene_name, "</code><br>",
        "<strong>Display Label:</strong> ", selected_node$label, "<br>",
        relationship_info,
        pathway_links,
        "<br><strong>=== PATHWAY CONTEXT ===</strong><br>",
        "This node represents a reference to another KEGG pathway.<br>",
        if (pathway_id != "") paste0("You can load the ", pathway_id, " pathway to explore it in detail.<br>") else "",
        "<br>",
        "<em>Tip: Click 'Load Pathway' button to navigate to the referenced pathway</em>",
        
        # JavaScript for pathway loading
        "<script>
        function loadPathwayFromNode(pathwayId) {
            Shiny.setInputValue('load_pathway_from_node', pathwayId, {priority: 'event'});
            document.getElementById('load_pathway_btn_", selected_node$id, "').innerHTML = 'Loading...';
            document.getElementById('load_pathway_btn_", selected_node$id, "').disabled = true;
        }
        </script>",
        "</div>"
    ))
}

# GENERIC NODE HANDLER
build_generic_node_info <- function(selected_node, incoming_edges, outgoing_edges, nodes_data) {
    relationship_info <- build_relationship_info(incoming_edges, outgoing_edges, nodes_data)
    
    return(paste0(
        "<div style='font-family: monospace; font-size: 14px; line-height: 1.6;'>",
        "<strong>=== NODE INFORMATION ===</strong><br>",
        "<strong>KEGG ID:</strong> <code>", selected_node$kegg_id, "</code><br>",
        "<strong>Original KEGG Entry:</strong> <code>", selected_node$gene_name, "</code><br>",
        "<strong>Display Label:</strong> ", selected_node$label, "<br>",
        relationship_info,
        "<br><strong>=== PATHWAY CONTEXT ===</strong><br>",
        "In KEGG, this node is referenced as: <code>", selected_node$kegg_id, "</code><br>",
        "<br>",
        "<em>Tip: Search for '", selected_node$kegg_id, "' in KEGG database</em>",
        "</div>"
    ))
}
