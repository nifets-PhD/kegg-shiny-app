# Helper functions for network visualization

create_network_visualization <- function(nodes, edges, layout = "forceAtlas2Based", 
                                       node_color = "#97C2FC", node_size = 25,
                                       edge_color = "#848484", edge_width = 2,
                                       repulsion_strength = 2000, spring_length = 200, 
                                       spring_strength = 0.05, enable_stabilization = TRUE,
                                       hierarchical_direction = "UD", hierarchical_spacing = 200) {
    
    # Validate inputs
    if (is.null(nodes) || nrow(nodes) == 0) {
        nodes <- data.frame(
            id = "placeholder",
            label = "No data available",
            stringsAsFactors = FALSE
        )
        edges <- data.frame(
            from = character(0),
            to = character(0),
            stringsAsFactors = FALSE
        )
    }
    
    if (is.null(edges)) {
        edges <- data.frame(
            from = character(0),
            to = character(0),
            stringsAsFactors = FALSE
        )
    }
    
    # Prepare nodes and edges
    vis_nodes <- prepare_nodes(nodes, color = node_color, size = node_size)
    vis_edges <- prepare_edges(edges, color = edge_color, width = edge_width, 
                              dense_network = nrow(nodes) > 50)
    
    # Create base network with essential settings only
    network <- visNetwork(vis_nodes, vis_edges) %>%
        visInteraction(
            navigationButtons = TRUE, 
            hover = TRUE,
            zoomView = TRUE,
            dragView = TRUE
        ) %>%
        visOptions(
            highlightNearest = list(enabled = TRUE, hover = TRUE, degree = 1),
            nodesIdSelection = list(
                enabled = TRUE,
                values = nodes$id[nodes$type == "gene" | is.na(nodes$type)]
            )
        ) %>%
        visEvents(select = "function(nodes) {
            Shiny.onInputChange('pathway_network_selected', nodes.nodes);
        }")
    
    # Apply layout-specific settings
    if (layout == "kegg_layout") {
        # Use KEGG's original coordinates if available
        if (!is.null(nodes$x) && !is.null(nodes$y) && any(!is.na(nodes$x))) {
            
            # Apply KEGG coordinates directly (no scaling)
            vis_nodes$x <- nodes$x
            vis_nodes$y <- nodes$y
            
            # Apply KEGG styling if available
            if (!is.null(nodes$kegg_bgcolor)) {
                for (i in seq_len(nrow(vis_nodes))) {
                    if (!is.na(nodes$kegg_bgcolor[i])) {
                        vis_nodes$color[i] <- nodes$kegg_bgcolor[i]
                    }
                    if (!is.na(nodes$kegg_shape[i])) {
                        # Map KEGG shapes to visNetwork shapes
                        kegg_shape <- nodes$kegg_shape[i]
                        vis_nodes$shape[i] <- switch(kegg_shape,
                                                   "rectangle" = "box",
                                                   "circle" = "circle",
                                                   "roundrectangle" = "box",
                                                   "ellipse" = "ellipse",
                                                   "box"  # default
                        )
                    }
                    if (!is.na(nodes$kegg_width[i])) {
                        # Scale width appropriately for visNetwork
                        vis_nodes$size[i] <- max(15, min(50, nodes$kegg_width[i] * 0.8))
                    }
                }
            }
            
            # Use fixed positions - this is key for maintaining KEGG layout
            network <- network %>%
                visNodes(
                    physics = FALSE, 
                    fixed = list(x = TRUE, y = TRUE),
                    chosen = FALSE
                ) %>%
                visPhysics(enabled = FALSE) %>%
                visLayout(improvedLayout = FALSE, randomSeed = NULL) %>%
                visOptions(
                    autoResize = TRUE,
                    height = "100%",
                    width = "100%",
                    clickToUse = FALSE
                )
            
            cat("Applied authentic KEGG layout with original coordinates\n")
            
            # Print coordinate info for debugging
            coord_range_x <- range(vis_nodes$x, na.rm = TRUE)
            coord_range_y <- range(vis_nodes$y, na.rm = TRUE)
            cat("KEGG coordinate ranges - X: [", coord_range_x[1], "-", coord_range_x[2], 
                "] Y: [", coord_range_y[1], "-", coord_range_y[2], "]\n")
                
        } else {
            # Fall back to hierarchical if no coordinates
            cat("No KEGG coordinates available, falling back to hierarchical layout\n")
            network <- network %>%
                visHierarchicalLayout(
                    direction = "UD",
                    levelSeparation = 150,
                    nodeSpacing = 100
                ) %>%
                visPhysics(enabled = FALSE)
        }
        
    } else if (layout == "hierarchical") {
        network <- network %>%
            visHierarchicalLayout(
                direction = hierarchical_direction,
                levelSeparation = hierarchical_spacing,
                nodeSpacing = 150
            ) %>%
            visPhysics(enabled = FALSE)
            
    } else if (layout == "forceAtlas2Based") {
        network <- network %>%
            visPhysics(
                enabled = TRUE,
                solver = "forceAtlas2Based",
                forceAtlas2Based = list(
                    gravitationalConstant = -repulsion_strength,
                    centralGravity = 0.01,
                    springLength = spring_length,
                    springConstant = spring_strength,
                    damping = 0.4,
                    avoidOverlap = 0.5
                ),
                stabilization = list(
                    enabled = enable_stabilization,
                    iterations = 1000,
                    updateInterval = 25
                )
            )
            
    } else if (layout == "physics") {
        network <- network %>%
            visPhysics(
                enabled = TRUE,
                solver = "barnesHut",
                barnesHut = list(
                    gravitationalConstant = -repulsion_strength,
                    centralGravity = 0.3,
                    springLength = spring_length,
                    springConstant = spring_strength,
                    damping = 0.09
                ),
                stabilization = list(
                    enabled = enable_stabilization,
                    iterations = 1000
                )
            )
            
    } else {
        # Circular or random layout
        network <- network %>%
            visLayout(randomSeed = 123) %>%
            visPhysics(enabled = FALSE)
    }
    
    return(network)
}

prepare_nodes <- function(nodes, color = "#97C2FC", size = 25) {
    if (is.null(nodes) || nrow(nodes) == 0) {
        # Return empty dataframe with correct structure for visNetwork
        return(data.frame(
            id = character(0),
            label = character(0),
            color = character(0),
            size = numeric(0),
            stringsAsFactors = FALSE
        ))
    }
    
    # Create a copy to avoid modifying the original
    vis_nodes <- data.frame(
        id = nodes$id,
        label = nodes$label,
        color = color,
        size = size,
        borderWidth = 2,
        borderWidthSelected = 4,
        font = list(color = "#000000", size = 12),
        stringsAsFactors = FALSE
    )
    
    # Use KEGG styling if available
    if (!is.null(nodes$kegg_bgcolor)) {
        for (i in seq_len(nrow(vis_nodes))) {
            if (!is.na(nodes$kegg_bgcolor[i])) {
                vis_nodes$color[i] <- nodes$kegg_bgcolor[i]
            }
            
            # Set font color based on KEGG foreground color
            if (!is.null(nodes$kegg_fgcolor) && !is.na(nodes$kegg_fgcolor[i])) {
                vis_nodes$font[[i]] <- list(color = nodes$kegg_fgcolor[i], size = 11)
            }
            
            # Set shape based on KEGG shape
            if (!is.null(nodes$kegg_shape) && !is.na(nodes$kegg_shape[i])) {
                kegg_shape <- nodes$kegg_shape[i]
                vis_nodes$shape[i] <- switch(kegg_shape,
                                           "rectangle" = "box",
                                           "circle" = "circle", 
                                           "roundrectangle" = "box",
                                           "ellipse" = "ellipse",
                                           "box"  # default
                )
            }
            
            # Set size based on KEGG dimensions
            if (!is.null(nodes$kegg_width) && !is.na(nodes$kegg_width[i])) {
                # Scale KEGG width to reasonable visNetwork size
                kegg_width <- nodes$kegg_width[i]
                vis_nodes$size[i] <- max(15, min(60, kegg_width * 0.6))
            }
        }
    }
    
    # Add any additional columns from original nodes (like hgnc_symbol, kegg_id, etc.)
    additional_cols <- setdiff(names(nodes), c("id", "label", "color", "size", "borderWidth", 
                                              "borderWidthSelected", "font", "shape"))
    for (col in additional_cols) {
        vis_nodes[[col]] <- nodes[[col]]
    }
    
    return(vis_nodes)
}

prepare_edges <- function(edges, color = "#848484", width = 2, dense_network = FALSE) {
    if (is.null(edges) || nrow(edges) == 0) {
        # Return empty dataframe with correct structure for visNetwork
        return(data.frame(
            from = character(0),
            to = character(0),
            color = character(0),
            width = numeric(0),
            arrows = character(0),
            stringsAsFactors = FALSE
        ))
    }
    
    # Create simple edge styling to avoid list assignment issues
    vis_edges <- data.frame(
        from = edges$from,
        to = edges$to,
        color = color,
        width = width,
        arrows = "to",
        stringsAsFactors = FALSE
    )
    
    # Add opacity for dense networks
    if (dense_network) {
        vis_edges$color <- paste0(color, "CC")  # Add alpha transparency
    }
    
    return(vis_edges)
}

apply_gene_annotations <- function(nodes, annotations, annotation_column, color_scheme = "viridis") {
    # Try to match annotations to nodes based on multiple ID types
    # First try HGNC symbols, then KEGG IDs, then original gene names
    
    # Get the first column of annotations (assumed to be gene identifiers)
    annotation_genes <- annotations[[1]]
    
    # Try matching with HGNC symbols first
    matched_annotations <- annotations[match(nodes$hgnc_symbol, annotation_genes), ]
    
    # For unmatched nodes, try KEGG IDs
    unmatched <- is.na(matched_annotations[[1]])
    if (any(unmatched)) {
        kegg_matches <- annotations[match(nodes$kegg_id[unmatched], annotation_genes), ]
        matched_annotations[unmatched, ] <- kegg_matches
    }
    
    # For still unmatched nodes, try original gene names
    still_unmatched <- is.na(matched_annotations[[1]])
    if (any(still_unmatched)) {
        gene_matches <- annotations[match(nodes$gene_name[still_unmatched], annotation_genes), ]
        matched_annotations[still_unmatched, ] <- gene_matches
    }
    
    # Get annotation values
    annotation_values <- matched_annotations[[annotation_column]]
    
    # Create color scale
    if (is.numeric(annotation_values)) {
        # Continuous color scale
        color_scale <- create_continuous_color_scale(annotation_values, color_scheme)
    } else {
        # Discrete color scale
        color_scale <- create_discrete_color_scale(annotation_values, color_scheme)
    }
    
    # Apply colors to nodes
    nodes$color <- color_scale$colors
    nodes$annotation_value <- annotation_values
    
    # Count successful matches
    successful_matches <- sum(!is.na(annotation_values))
    cat("Successfully matched", successful_matches, "out of", nrow(nodes), "nodes with annotations\n")
    
    return(list(
        nodes = nodes,
        color_scale = color_scale
    ))
}

create_continuous_color_scale <- function(values, scheme = "viridis") {
    # Remove NA values for scaling
    clean_values <- values[!is.na(values)]
    
    if (length(clean_values) == 0) {
        return(list(colors = rep("#97C2FC", length(values))))
    }
    
    # Normalize values to 0-1 range
    normalized <- (values - min(clean_values, na.rm = TRUE)) / 
                  (max(clean_values, na.rm = TRUE) - min(clean_values, na.rm = TRUE))
    
    # Generate colors
    colors <- switch(scheme,
        "viridis" = viridisLite::viridis(100)[round(normalized * 99) + 1],
        "plasma" = viridisLite::plasma(100)[round(normalized * 99) + 1],
        "RdBu" = RColorBrewer::brewer.pal(11, "RdBu")[round(normalized * 10) + 1],
        "RdYlBu" = RColorBrewer::brewer.pal(11, "RdYlBu")[round(normalized * 10) + 1],
        "Spectral" = RColorBrewer::brewer.pal(11, "Spectral")[round(normalized * 10) + 1],
        viridisLite::viridis(100)[round(normalized * 99) + 1]
    )
    
    # Handle NA values
    colors[is.na(values)] <- "#CCCCCC"
    
    return(list(
        colors = colors,
        scale_type = "continuous",
        min_value = min(clean_values, na.rm = TRUE),
        max_value = max(clean_values, na.rm = TRUE),
        scheme = scheme
    ))
}

create_discrete_color_scale <- function(values, scheme = "viridis") {
    unique_values <- unique(values[!is.na(values)])
    n_colors <- length(unique_values)
    
    if (n_colors == 0) {
        return(list(colors = rep("#97C2FC", length(values))))
    }
    
    # Generate colors
    palette_colors <- switch(scheme,
        "viridis" = viridisLite::viridis(n_colors),
        "plasma" = viridisLite::plasma(n_colors),
        "RdBu" = RColorBrewer::brewer.pal(min(n_colors, 11), "RdBu"),
        "RdYlBu" = RColorBrewer::brewer.pal(min(n_colors, 11), "RdYlBu"),
        "Spectral" = RColorBrewer::brewer.pal(min(n_colors, 11), "Spectral"),
        viridisLite::viridis(n_colors)
    )
    
    # Map values to colors
    color_mapping <- setNames(palette_colors, unique_values)
    colors <- color_mapping[values]
    
    # Handle NA values
    colors[is.na(values)] <- "#CCCCCC"
    
    return(list(
        colors = colors,
        scale_type = "discrete",
        mapping = color_mapping,
        scheme = scheme
    ))
}

create_color_legend <- function(color_scale, annotation_name) {
    if (color_scale$scale_type == "continuous") {
        # Create continuous legend
        values <- seq(color_scale$min_value, color_scale$max_value, length.out = 100)
        colors <- switch(color_scale$scheme,
            "viridis" = viridisLite::viridis(100),
            "plasma" = viridisLite::plasma(100),
            "RdBu" = colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(100),
            "RdYlBu" = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(100),
            "Spectral" = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(100),
            viridisLite::viridis(100)
        )
        
        plot_ly(
            x = values,
            y = rep(1, length(values)),
            type = "scatter",
            mode = "markers",
            marker = list(
                color = colors,
                size = 20,
                line = list(width = 0)
            ),
            showlegend = FALSE
        ) %>%
            layout(
                title = paste("Color Scale:", annotation_name),
                xaxis = list(title = annotation_name),
                yaxis = list(visible = FALSE),
                margin = list(l = 50, r = 50, t = 50, b = 50)
            )
    } else {
        # Create discrete legend
        values <- names(color_scale$mapping)
        colors <- as.character(color_scale$mapping)
        
        plot_ly(
            x = values,
            y = rep(1, length(values)),
            type = "scatter",
            mode = "markers",
            marker = list(
                color = colors,
                size = 20,
                line = list(width = 0)
            ),
            text = values,
            showlegend = FALSE
        ) %>%
            layout(
                title = paste("Color Scale:", annotation_name),
                xaxis = list(title = annotation_name),
                yaxis = list(visible = FALSE),
                margin = list(l = 50, r = 50, t = 50, b = 50)
            )
    }
}

#' Download and cache complete HGNC to UniProt mapping
#' @param force_refresh logical, whether to force refresh the cache
#' @return data.frame with HGNC symbols and UniProt IDs
download_hgnc_uniprot_mapping <- function(force_refresh = FALSE) {
    # Define cache file path
    mapping_cache_file <- file.path("data", "hgnc_uniprot_mapping.rds")
    
    # Create data directory if it doesn't exist
    if (!dir.exists("data")) {
        dir.create("data", recursive = TRUE)
    }
    
    # Check if cache exists and is recent (less than 30 days old)
    if (file.exists(mapping_cache_file) && !force_refresh) {
        file_age <- difftime(Sys.time(), file.mtime(mapping_cache_file), units = "days")
        if (file_age < 30) {
            cat("Loading cached HGNC-UniProt mapping (", round(file_age, 1), "days old)...\n")
            mapping_data <- readRDS(mapping_cache_file)
            cat("Loaded", nrow(mapping_data), "cached mappings\n")
            return(mapping_data)
        } else {
            cat("Cache is", round(file_age, 1), "days old, refreshing...\n")
        }
    }
    
    tryCatch({
        if (!requireNamespace("biomaRt", quietly = TRUE)) {
            stop("biomaRt package is required. Please install it with: BiocManager::install('biomaRt')")
        }
        
        cat("Downloading complete HGNC to UniProt mapping from Ensembl...\n")
        
        # Connect to Ensembl BioMart
        cat("Connecting to Ensembl BioMart...\n")
        ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        
        # Get all HGNC symbols and their UniProt mappings
        # This downloads the complete mapping for all human genes
        cat("Querying BioMart for gene mappings...\n")
        mapping_result <- biomaRt::getBM(
            attributes = c("hgnc_symbol", "uniprotswissprot", "uniprotsptrembl", 
                          "ensembl_gene_id", "description"),
            mart = ensembl
        )
        
        cat("Downloaded", nrow(mapping_result), "gene records from Ensembl\n")
        
        if (nrow(mapping_result) > 0) {
            cat("Processing", nrow(mapping_result), "records...\n")
            
            # First, filter for valid HGNC symbols
            has_hgnc <- !is.na(mapping_result$hgnc_symbol) & mapping_result$hgnc_symbol != ""
            cat("Records with HGNC symbols:", sum(has_hgnc), "\n")
            
            # Then filter for at least one UniProt ID
            has_swissprot <- !is.na(mapping_result$uniprotswissprot) & mapping_result$uniprotswissprot != ""
            has_trembl <- !is.na(mapping_result$uniprotsptrembl) & mapping_result$uniprotsptrembl != ""
            has_uniprot <- has_swissprot | has_trembl
            
            cat("Records with SwissProt IDs:", sum(has_swissprot, na.rm = TRUE), "\n")
            cat("Records with TrEMBL IDs:", sum(has_trembl, na.rm = TRUE), "\n")
            cat("Records with any UniProt ID:", sum(has_uniprot, na.rm = TRUE), "\n")
            
            # Combine filters
            mapping_filtered <- mapping_result[has_hgnc & has_uniprot, ]
            cat("Records after filtering:", nrow(mapping_filtered), "\n")
            
            # Create a clean mapping with preferred UniProt ID
            mapping_clean <- data.frame(
                hgnc_symbol = mapping_filtered$hgnc_symbol,
                uniprot_id = ifelse(
                    !is.na(mapping_filtered$uniprotswissprot) & mapping_filtered$uniprotswissprot != "",
                    mapping_filtered$uniprotswissprot,  # Prefer Swiss-Prot
                    mapping_filtered$uniprotsptrembl    # Fall back to TrEMBL
                ),
                uniprot_swissprot = mapping_filtered$uniprotswissprot,
                uniprot_trembl = mapping_filtered$uniprotsptrembl,
                ensembl_gene_id = mapping_filtered$ensembl_gene_id,
                description = mapping_filtered$description,
                stringsAsFactors = FALSE
            )
            
            # Remove entries where uniprot_id is still NA or empty
            mapping_clean <- mapping_clean[
                !is.na(mapping_clean$uniprot_id) & mapping_clean$uniprot_id != "",
            ]
            
            # For genes with multiple UniProt IDs, keep the first one (usually the primary)
            mapping_clean <- mapping_clean[!duplicated(mapping_clean$hgnc_symbol), ]
            
            cat("Created clean mapping for", nrow(mapping_clean), "genes with UniProt IDs\n")
            
            # Save to cache
            saveRDS(mapping_clean, mapping_cache_file)
            cat("Saved mapping to cache:", mapping_cache_file, "\n")
            
            # Also save as TSV for human inspection
            tsv_file <- file.path("data", "hgnc_uniprot_mapping.tsv")
            write.table(mapping_clean, tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)
            cat("Also saved as TSV file:", tsv_file, "\n")
            
            return(mapping_clean)
        } else {
            stop("No mapping data retrieved from Ensembl")
        }
        
    }, error = function(e) {
        cat("Error downloading HGNC-UniProt mapping:", e$message, "\n")
        
        # Return empty dataframe if download fails
        return(data.frame(
            hgnc_symbol = character(0),
            uniprot_id = character(0),
            uniprot_swissprot = character(0),
            uniprot_trembl = character(0),
            ensembl_gene_id = character(0),
            description = character(0),
            stringsAsFactors = FALSE
        ))
    })
}

#' Get UniProt ID for HGNC gene symbol using cached mapping
#' @param hgnc_symbol HGNC gene symbol
#' @param mapping_data optional cached mapping data
#' @return UniProt ID or NULL if not found
get_uniprot_id <- function(hgnc_symbol, mapping_data = NULL) {
    # Load mapping data if not provided
    if (is.null(mapping_data)) {
        mapping_data <- download_hgnc_uniprot_mapping()
    }
    
    # Find the UniProt ID for this gene
    match_idx <- which(mapping_data$hgnc_symbol == hgnc_symbol)
    
    if (length(match_idx) > 0) {
        uniprot_id <- mapping_data$uniprot_id[match_idx[1]]  # Take first match
        if (!is.na(uniprot_id) && uniprot_id != "") {
            return(uniprot_id)
        }
    }
    
    return(NULL)
}

#' Create HTML for UniProt structure link
#' @param uniprot_id UniProt ID
#' @return HTML string for the link
create_uniprot_structure_link <- function(uniprot_id) {
    if (is.null(uniprot_id) || is.na(uniprot_id) || uniprot_id == "") {
        return("")
    }
    
    uniprot_url <- paste0("https://www.uniprot.org/uniprotkb/", uniprot_id, "/entry#structure")
    return(paste0('<a href="', uniprot_url, '" target="_blank" style="color: #0066cc; text-decoration: underline;">Go to UniProt Structure</a>'))
}

#' Initialize HGNC-UniProt mapping (run on app startup)
#' @param force_refresh logical, whether to force refresh the cache
#' @return invisible success status
initialize_uniprot_mapping <- function(force_refresh = FALSE) {
    tryCatch({
        cat("Initializing HGNC-UniProt mapping...\n")
        mapping_data <- download_hgnc_uniprot_mapping(force_refresh = force_refresh)
        cat("Successfully initialized UniProt mapping for", nrow(mapping_data), "genes\n")
        return(invisible(TRUE))
    }, error = function(e) {
        warning("Failed to initialize UniProt mapping: ", e$message)
        return(invisible(FALSE))
    })
}
