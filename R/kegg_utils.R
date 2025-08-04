# Helper functions for KEGG data processing
get_kegg_node_coordinates <- function(pathway_id) {
    tryCatch({
        # KEGG provides coordinate data in their KGML files
        # This function extracts x,y coordinates and styling information
        
        cache_dir <- "kegg_cache"
        kgml_file <- file.path(cache_dir, paste0(pathway_id, ".xml"))
        
        cat("Looking for KGML file:", kgml_file, "\n")
        
        if (!file.exists(kgml_file)) {
            cat("KGML file does not exist:", kgml_file, "\n")
            return(NULL)
        } else {
            cat("KGML file found, size:", file.size(kgml_file), "bytes\n")
        }
        
        # Parse XML to extract coordinate information
        if (requireNamespace("xml2", quietly = TRUE)) {
            cat("xml2 package available, parsing XML...\n")
            xml_doc <- xml2::read_xml(kgml_file)
            
            # Extract entry nodes with graphics information
            entries <- xml2::xml_find_all(xml_doc, "//entry")
            
            cat("Found", length(entries), "entries in XML\n")
            
            if (length(entries) > 0) {
                coords_list <- list()
                
                for (i in seq_along(entries)) {
                    entry <- entries[[i]]
                    entry_id <- xml2::xml_attr(entry, "id")
                    entry_name <- xml2::xml_attr(entry, "name")
                    entry_type <- xml2::xml_attr(entry, "type")
                    
                    cat("Processing entry", i, "- ID:", entry_id, "Name:", entry_name, "Type:", entry_type, "\n")
                    
                    # Find graphics element
                    graphics <- xml2::xml_find_first(entry, ".//graphics")
                    
                    if (!is.na(graphics)) {
                        # Extract all graphics attributes
                        x <- as.numeric(xml2::xml_attr(graphics, "x"))
                        y <- as.numeric(xml2::xml_attr(graphics, "y"))
                        width <- as.numeric(xml2::xml_attr(graphics, "width"))
                        height <- as.numeric(xml2::xml_attr(graphics, "height"))
                        shape_type <- xml2::xml_attr(graphics, "type")
                        fgcolor <- xml2::xml_attr(graphics, "fgcolor")
                        bgcolor <- xml2::xml_attr(graphics, "bgcolor")
                        label_text <- xml2::xml_attr(graphics, "name")
                        
                        # Extract first entry from graphics name (before first comma) for gene aliases
                        if (!is.na(label_text) && entry_type == "gene") {
                            # Split by comma and take the first entry, clean up any trailing "..."
                            first_alias <- strsplit(label_text, ",", fixed = TRUE)[[1]][1]
                            first_alias <- trimws(first_alias)  # Remove whitespace
                            first_alias <- gsub("\\.\\.\\.$", "", first_alias)  # Remove trailing "..."
                            label_text <- first_alias
                        }
                        
                        cat("  Graphics - X:", x, "Y:", y, "Shape:", shape_type, "BgColor:", bgcolor, "\n")
                        
                        if (!is.na(x) && !is.na(y)) {
                            coords_list[[paste0("entry_", entry_id)]] <- list(
                                entry_id = entry_id,
                                kegg_name = entry_name,
                                entry_type = entry_type,
                                x = x,
                                y = y,
                                width = ifelse(is.na(width), 46, width),
                                height = ifelse(is.na(height), 17, height),
                                shape = ifelse(is.na(shape_type), "rectangle", shape_type),
                                fgcolor = ifelse(is.na(fgcolor), "#000000", fgcolor),
                                bgcolor = ifelse(is.na(bgcolor), "#BFFFBF", bgcolor),
                                label = ifelse(is.na(label_text), entry_name, label_text)
                            )
                            cat("  Added to coords_list\n")
                        } else {
                            cat("  Skipping - missing coordinates\n")
                        }
                    } else {
                        cat("  No graphics element found\n")
                    }
                }
                
                if (length(coords_list) > 0) {
                    coords_df <- do.call(rbind, lapply(coords_list, function(x) {
                        data.frame(
                            entry_id = x$entry_id,
                            kegg_name = x$kegg_name,
                            entry_type = x$entry_type,
                            x = x$x,
                            y = x$y,
                            width = x$width,
                            height = x$height,
                            shape = x$shape,
                            fgcolor = x$fgcolor,
                            bgcolor = x$bgcolor,
                            label = x$label,
                            stringsAsFactors = FALSE
                        )
                    }))
                    
                    cat("SUCCESS: Extracted detailed KEGG layout for", nrow(coords_df), "entries\n")
                    cat("Coordinate ranges - X: [", min(coords_df$x), "-", max(coords_df$x), 
                        "] Y: [", min(coords_df$y), "-", max(coords_df$y), "]\n")
                    
                    return(coords_df)
                } else {
                    cat("No valid coordinates found in any entry\n")
                }
            } else {
                cat("No entries found in XML\n")
            }
        } else {
            cat("xml2 package not available - install with install.packages('xml2') for coordinate extraction\n")
        }
        
        return(NULL)
        
    }, error = function(e) {
        cat("ERROR in extracting KEGG coordinates:", e$message, "\n")
        return(NULL)
    })
}

get_cache_info <- function() {
    cache_dir <- "kegg_cache"
    if (!dir.exists(cache_dir)) {
        return(list(
            exists = FALSE,
            file_count = 0,
            total_size = 0,
            files = character(0)
        ))
    }
    
    files <- list.files(cache_dir, full.names = TRUE)
    file_info <- file.info(files)
    
    return(list(
        exists = TRUE,
        file_count = length(files),
        total_size = sum(file_info$size, na.rm = TRUE),
        files = basename(files),
        cache_dir = cache_dir
    ))
}

clear_cache <- function() {
    cache_dir <- "kegg_cache"
    if (dir.exists(cache_dir)) {
        files <- list.files(cache_dir, full.names = TRUE)
        unlink(files)
        cat("Cleared", length(files), "cache files\n")
        return(length(files))
    }
    return(0)
}

convert_kegg_to_hgnc <- function(kegg_ids) {
    tryCatch({
        cat("Converting KEGG IDs to HGNC symbols...\n")
        
        # Clean KEGG IDs - remove "hsa:" prefix if present
        clean_kegg_ids <- gsub("^hsa:", "", kegg_ids)
        clean_kegg_ids <- gsub("^.*:", "", clean_kegg_ids) # Remove any other prefixes
        
        # Filter out non-numeric IDs (keep only gene IDs)
        numeric_ids <- clean_kegg_ids[grepl("^\\d+$", clean_kegg_ids)]
        
        if (length(numeric_ids) == 0) {
            cat("No valid KEGG gene IDs found\n")
            return(data.frame(
                kegg_id = kegg_ids,
                hgnc_symbol = kegg_ids,
                entrez_id = NA,
                stringsAsFactors = FALSE
            ))
        }
        
        cat("Found", length(numeric_ids), "valid KEGG gene IDs\n")
        
        # Connect to Ensembl BioMart
        ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        
        # Query BioMart to get HGNC symbols
        # Note: KEGG IDs are actually Entrez Gene IDs
        conversion_result <- getBM(
            attributes = c("entrezgene_id", "hgnc_symbol", "external_gene_name"),
            filters = "entrezgene_id",
            values = numeric_ids,
            mart = ensembl
        )
        
        cat("BioMart returned", nrow(conversion_result), "mappings\n")
        
        # Create mapping dataframe
        mapping_df <- data.frame(
            original_id = kegg_ids,
            kegg_id = clean_kegg_ids,
            hgnc_symbol = NA,
            entrez_id = NA,
            stringsAsFactors = FALSE
        )
        
        # Fill in the mappings
        for (i in seq_len(nrow(mapping_df))) {
            kegg_id <- mapping_df$kegg_id[i]
            if (kegg_id %in% conversion_result$entrezgene_id) {
                match_row <- conversion_result[conversion_result$entrezgene_id == kegg_id, ]
                if (nrow(match_row) > 0) {
                    # Prefer hgnc_symbol, fall back to external_gene_name
                    symbol <- match_row$hgnc_symbol[1]
                    if (is.na(symbol) || symbol == "") {
                        symbol <- match_row$external_gene_name[1]
                    }
                    mapping_df$hgnc_symbol[i] <- symbol
                    mapping_df$entrez_id[i] <- kegg_id
                }
            }
        }
        
        # For unmapped IDs, use the original ID
        unmapped <- is.na(mapping_df$hgnc_symbol) | mapping_df$hgnc_symbol == ""
        mapping_df$hgnc_symbol[unmapped] <- mapping_df$original_id[unmapped]
        
        cat("Successfully mapped", sum(!unmapped), "out of", length(kegg_ids), "IDs\n")
        
        return(mapping_df)
        
    }, error = function(e) {
        cat("Error in KEGG to HGNC conversion:", e$message, "\n")
        # Return original IDs if conversion fails
        return(data.frame(
            original_id = kegg_ids,
            kegg_id = kegg_ids,
            hgnc_symbol = kegg_ids,
            entrez_id = NA,
            stringsAsFactors = FALSE
        ))
    })
}

get_kegg_pathways <- function() {
    tryCatch({
        # Get human pathways
        pathways <- keggList("pathway", "hsa")
        
        # Convert to data frame
        pathway_df <- data.frame(
            pathway_id = names(pathways),
            pathway_name = as.character(pathways),
            stringsAsFactors = FALSE
        )
        
        # Clean up pathway names - remove the "path:" prefix and " - Homo sapiens (human)" suffix
        pathway_df$pathway_name <- gsub("^path:hsa\\d+\\s+", "", pathway_df$pathway_name)
        pathway_df$pathway_name <- gsub("\\s+-\\s+Homo sapiens \\(human\\)$", "", pathway_df$pathway_name)
        
        # Add descriptions (use the pathway name)
        pathway_df$description <- pathway_df$pathway_name
        
        # Debug: print first few rows
        cat("First few pathways loaded:\n")
        print(head(pathway_df, 3))
        
        return(pathway_df)
        
    }, error = function(e) {
        stop(paste("Error fetching KEGG pathways:", e$message))
    })
}

search_pathways <- function(pathways_list, search_term) {
    if (is.null(pathways_list) || search_term == "") {
        return(pathways_list)
    }
    
    matches <- grepl(search_term, pathways_list$pathway_name, ignore.case = TRUE)
    return(pathways_list[matches, ])
}

filter_pathways_by_category <- function(pathways_list, category) {
    if (is.null(pathways_list)) {
        return(NULL)
    }
    
    # Filter pathways based on category prefix
    category_matches <- grepl(paste0("^hsa", category), pathways_list$pathway_id)
    return(pathways_list[category_matches, ])
}

load_kegg_pathway <- function(pathway_id) {
    tryCatch({
        # Debug print
        cat("Original pathway_id:", pathway_id, "\n")
        
        # Clean pathway ID
        pathway_id <- gsub("^path:", "", pathway_id)
        
        # Ensure pathway_id is not NA and is properly formatted  
        if (is.na(pathway_id) || pathway_id == "" || pathway_id == "hsaNA") {
            stop(paste("Invalid pathway ID:", pathway_id))
        }
        
        # Ensure pathway_id starts with "hsa" and is properly formatted
        if (!grepl("^hsa\\d+", pathway_id)) {
            stop(paste("Invalid pathway ID format:", pathway_id))
        }
        
        cat("Cleaned pathway_id:", pathway_id, "\n")
        
        # Create cache directory if it doesn't exist
        cache_dir <- "kegg_cache"
        if (!dir.exists(cache_dir)) {
            dir.create(cache_dir, recursive = TRUE)
            cat("Created cache directory:", cache_dir, "\n")
        }
        
        # Define cache file paths
        kgml_cache_file <- file.path(cache_dir, paste0(pathway_id, ".xml"))
        rds_cache_file <- file.path(cache_dir, paste0(pathway_id, ".rds"))
        
        # Check if we have a processed RDS cache file
        if (file.exists(rds_cache_file)) {
            cat("Loading pathway from RDS cache:", rds_cache_file, "\n")
            cached_data <- readRDS(rds_cache_file)
            
            # Verify cache integrity
            if (is.list(cached_data) && 
                !is.null(cached_data$nodes) && 
                !is.null(cached_data$edges)) {
                cat("Successfully loaded cached pathway data\n")
                return(cached_data)
            } else {
                cat("Cache file corrupted, removing and re-downloading\n")
                unlink(rds_cache_file)
                if (file.exists(kgml_cache_file)) unlink(kgml_cache_file)
            }
        }
        
        # Check if we have raw KGML cache file
        temp_file <- kgml_cache_file
        download_needed <- TRUE
        
        if (file.exists(kgml_cache_file) && file.size(kgml_cache_file) > 0) {
            cat("Using cached KGML file:", kgml_cache_file, "\n")
            download_needed <- FALSE
        }
        
        # Download if needed
        if (download_needed) {
            cat("Downloading pathway data to cache...\n")
            max_retries <- 3
            retry_count <- 0
            success <- FALSE
            
            while (retry_count < max_retries && !success) {
                tryCatch({
                    cat("Attempt", retry_count + 1, "to download pathway data...\n")
                    
                    # Set a longer timeout and try to download
                    options(timeout = 120)  # 2 minutes timeout
                    pathway_kgml <- retrieveKGML(pathway_id, organism = "hsa", destfile = temp_file)
                    
                    if (file.exists(temp_file) && file.size(temp_file) > 0) {
                        success <- TRUE
                        cat("Download successful! Cached to:", temp_file, "\n")
                    } else {
                        retry_count <- retry_count + 1
                        Sys.sleep(2)  # Wait 2 seconds before retry
                    }
                }, error = function(e) {
                    cat("Download attempt failed:", e$message, "\n")
                    retry_count <<- retry_count + 1
                    if (retry_count < max_retries) {
                        Sys.sleep(2)  # Wait 2 seconds before retry
                    }
                })
            }
            
            if (!success) {
                # If download fails, create a mock pathway for demonstration
                cat("Creating mock pathway data due to download failure\n")
                
                mock_nodes <- data.frame(
                    id = paste0("node_", 1:10),
                    label = paste("Gene", 1:10),
                    type = "gene",
                    gene_name = paste("GENE", 1:10),
                    hgnc_symbol = paste("GENE", 1:10),
                    kegg_id = paste0("mock_", 1:10),
                    description = paste("Mock gene", 1:10, "from pathway", pathway_id),
                    stringsAsFactors = FALSE
                )
                
                mock_edges <- data.frame(
                    from = c("node_1", "node_2", "node_3", "node_4", "node_5", 
                            "node_6", "node_7", "node_8", "node_9"),
                    to = c("node_2", "node_3", "node_4", "node_5", "node_6",
                          "node_7", "node_8", "node_9", "node_10"),
                    stringsAsFactors = FALSE
                )
                
                mock_data <- list(
                    graph = NULL,
                    nodes = mock_nodes,
                    edges = mock_edges
                )
                
                # Cache the mock data
                saveRDS(mock_data, rds_cache_file)
                cat("Cached mock data to:", rds_cache_file, "\n")
                
                return(mock_data)
            }
        }
        
        # Parse the KGML file
        cat("Parsing KGML file...\n")
        pathway_graph <- parseKGML2Graph(temp_file)
        
        if (is.null(pathway_graph)) {
            stop("Failed to parse pathway graph")
        }
        
        # Extract nodes and edges
        nodes_data <- extract_nodes_from_graph(pathway_graph)
        edges_data <- extract_edges_from_graph(pathway_graph)
        
        # Try to get KEGG coordinate data for better layout
        kegg_coords <- get_kegg_node_coordinates(pathway_id)
        if (!is.null(kegg_coords)) {
            cat("SUCCESS: KEGG coordinates extracted, creating nodes from KEGG XML data\n")
            
            # Create nodes directly from KEGG coordinate data instead of relying on KEGGgraph
            nodes_data <- data.frame(
                id = paste0("kegg_", kegg_coords$entry_id),
                label = kegg_coords$label,
                type = kegg_coords$entry_type,
                gene_name = kegg_coords$kegg_name,
                hgnc_symbol = kegg_coords$label,
                kegg_id = kegg_coords$entry_id,
                description = paste("KEGG", kegg_coords$entry_type, "entry"),
                x = kegg_coords$x,
                y = kegg_coords$y,
                kegg_width = kegg_coords$width,
                kegg_height = kegg_coords$height,
                kegg_shape = kegg_coords$shape,
                kegg_bgcolor = kegg_coords$bgcolor,
                kegg_fgcolor = kegg_coords$fgcolor,
                kegg_type = kegg_coords$entry_type,
                kegg_label = kegg_coords$label,
                stringsAsFactors = FALSE
            )
            
            # For genes, we already have the correct gene symbols from XML graphics names
            # Skip HGNC conversion since we have better data directly from KEGG
            gene_rows <- nodes_data$type == "gene"
            if (any(gene_rows)) {
                cat("Using gene symbols directly from KEGG XML graphics names for", sum(gene_rows), "genes\n")
                
                # Keep the gene symbols from XML graphics names as the main labels
                # The kegg_coords$label already contains the first gene symbol from graphics name
                # nodes_data$label and nodes_data$hgnc_symbol are already set correctly above
                
                # Set kegg_id to the first numeric ID from the entry name for reference
                for (i in which(gene_rows)) {
                    entry_name <- nodes_data$gene_name[i]
                    # Extract first numeric ID from entry name like "hsa:2475 hsa:57521 hsa:84335"
                    numeric_ids <- regmatches(entry_name, gregexpr("\\d+", entry_name))[[1]]
                    if (length(numeric_ids) > 0) {
                        nodes_data$kegg_id[i] <- numeric_ids[1]
                    }
                }
            }
            
            cat("Created", nrow(nodes_data), "nodes from KEGG data:", 
                sum(nodes_data$type == "gene"), "genes,",
                sum(nodes_data$type == "compound"), "compounds,",
                sum(nodes_data$type == "map"), "pathway maps\n")
            
            # Also get edges from KEGG XML
            kegg_edges <- tryCatch({
                xml_file <- file.path("kegg_cache", paste0(pathway_id, ".xml"))
                doc <- xml2::read_xml(xml_file)
                
                # Get all relation elements
                relations <- xml2::xml_find_all(doc, "//relation")
                cat("Found", length(relations), "relations in KEGG XML\n")
                
                if (length(relations) > 0) {
                    edges_list <- list()
                    
                    for (i in seq_along(relations)) {
                        relation <- relations[[i]]
                        
                        entry1 <- xml2::xml_attr(relation, "entry1")
                        entry2 <- xml2::xml_attr(relation, "entry2")
                        rel_type <- xml2::xml_attr(relation, "type")
                        
                        # Get subtype information
                        subtypes <- xml2::xml_find_all(relation, ".//subtype")
                        subtype_names <- sapply(subtypes, function(st) xml2::xml_attr(st, "name"))
                        subtype_values <- sapply(subtypes, function(st) xml2::xml_attr(st, "value"))
                        
                        edges_list[[i]] <- data.frame(
                            from = paste0("kegg_", entry1),
                            to = paste0("kegg_", entry2),
                            relation_type = ifelse(is.na(rel_type), "unknown", rel_type),
                            subtype = paste(subtype_names, collapse = "; "),
                            subtype_value = paste(subtype_values, collapse = "; "),
                            stringsAsFactors = FALSE
                        )
                    }
                    
                    edges_df <- do.call(rbind, edges_list)
                    cat("Successfully extracted", nrow(edges_df), "edge relations\n")
                    edges_df
                } else {
                    data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE)
                }
            }, error = function(e) {
                cat("Error extracting KEGG edges:", e$message, "\n")
                data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE)
            })
            
            edges_data <- kegg_edges
            
        } else {
            cat("No KEGG coordinates available, using KEGGgraph extraction\n")
            
            # Fallback to original method
            nodes_data <- extract_nodes_from_graph(pathway_graph)
            edges_data <- extract_edges_from_graph(pathway_graph)
        }
        
        # Ensure we have at least some nodes
        if (nrow(nodes_data) == 0) {
            nodes_data <- data.frame(
                id = c("node1", "node2"),
                label = c("Gene 1", "Gene 2"),
                type = "gene",
                gene_name = c("Gene 1", "Gene 2"),
                hgnc_symbol = c("Gene 1", "Gene 2"),
                kegg_id = c("1", "2"),
                description = "KEGG pathway component",
                stringsAsFactors = FALSE
            )
            edges_data <- data.frame(
                from = "node1",
                to = "node2",
                stringsAsFactors = FALSE
            )
        }
        
        # Create the final data structure
        pathway_data <- list(
            graph = pathway_graph,
            nodes = nodes_data,
            edges = edges_data,
            kegg_coordinates = kegg_coords
        )
        
        # Cache the processed data
        cat("Caching processed pathway data to:", rds_cache_file, "\n")
        saveRDS(pathway_data, rds_cache_file)
        
        cat("Pathway loading complete. Files cached for future use.\n")
        return(pathway_data)
        
    }, error = function(e) {
        stop(paste("Error loading pathway:", e$message))
    })
}

extract_nodes_from_graph <- function(graph) {
    if (is.null(graph) || length(nodes(graph)) == 0) {
        return(data.frame(
            id = character(0),
            label = character(0),
            type = character(0),
            gene_name = character(0),
            hgnc_symbol = character(0),
            kegg_id = character(0),
            description = character(0),
            stringsAsFactors = FALSE
        ))
    }
    
    node_names <- nodes(graph)
    
    # Create nodes dataframe with proper handling of node attributes
    nodes_df <- data.frame(
        id = node_names,
        label = node_names,
        type = "gene",
        gene_name = node_names,
        hgnc_symbol = node_names,
        kegg_id = node_names,
        description = "KEGG gene",
        stringsAsFactors = FALSE
    )
    
    # Try to get node attributes if available
    tryCatch({
        node_data <- nodeData(graph)
        if (!is.null(node_data) && length(node_data) > 0) {
            for (i in seq_along(node_names)) {
                node_id <- node_names[i]
                if (node_id %in% names(node_data)) {
                    attrs <- node_data[[node_id]]
                    if (!is.null(attrs$label)) {
                        nodes_df$label[i] <- attrs$label
                    }
                    if (!is.null(attrs$type)) {
                        nodes_df$type[i] <- attrs$type
                    }
                }
            }
        }
    }, error = function(e) {
        # Continue with default values if node data extraction fails
    })
    
    # Convert KEGG IDs to HGNC symbols
    if (nrow(nodes_df) > 0) {
        cat("Converting node IDs to HGNC symbols...\n")
        hgnc_mapping <- convert_kegg_to_hgnc(nodes_df$gene_name)
        
        # Update nodes with HGNC symbols
        nodes_df$hgnc_symbol <- hgnc_mapping$hgnc_symbol
        nodes_df$kegg_id <- hgnc_mapping$kegg_id
        
        # Create informative labels showing both KEGG and HGNC IDs
        # Prioritize KEGG ID for pathway context, show HGNC as secondary
        nodes_df$label <- ifelse(
            !is.na(hgnc_mapping$hgnc_symbol) & hgnc_mapping$hgnc_symbol != hgnc_mapping$kegg_id,
            paste0(hgnc_mapping$kegg_id, "\n(", hgnc_mapping$hgnc_symbol, ")"),
            hgnc_mapping$kegg_id
        )
        
        # Update description with comprehensive information
        nodes_df$description <- paste0(
            "KEGG ID: ", hgnc_mapping$kegg_id, "\n",
            "HGNC Symbol: ", hgnc_mapping$hgnc_symbol, "\n",
            "Original ID: ", hgnc_mapping$original_id
        )
    }
    
    return(nodes_df)
}

extract_edges_from_graph <- function(graph) {
    if (is.null(graph) || length(nodes(graph)) == 0) {
        return(data.frame(
            from = character(0),
            to = character(0),
            stringsAsFactors = FALSE
        ))
    }
    
    edges_df <- data.frame(
        from = character(0),
        to = character(0),
        stringsAsFactors = FALSE
    )
    
    tryCatch({
        edge_list <- edges(graph)
        
        if (!is.null(edge_list) && length(edge_list) > 0) {
            for (node in names(edge_list)) {
                targets <- edge_list[[node]]
                if (!is.null(targets) && length(targets) > 0) {
                    # Filter out any NA or invalid targets
                    valid_targets <- targets[!is.na(targets) & targets != ""]
                    
                    if (length(valid_targets) > 0) {
                        node_edges <- data.frame(
                            from = rep(node, length(valid_targets)),
                            to = valid_targets,
                            stringsAsFactors = FALSE
                        )
                        edges_df <- rbind(edges_df, node_edges)
                    }
                }
            }
        }
    }, error = function(e) {
        # Return empty dataframe if edge extraction fails
        warning(paste("Could not extract edges:", e$message))
    })
    
    return(edges_df)
}

# Simplified KEGG-specific network visualization
create_kegg_network_visualization <- function(nodes, edges, show_labels = TRUE, show_edges = TRUE, coloring_mode = "kegg_default") {
    
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
    
    # Apply phylostratum coloring if requested
    if (coloring_mode == "phylostratum") {
        tryCatch({
            nodes <- apply_phylostratum_coloring(nodes)
        }, error = function(e) {
            warning("Failed to apply phylostratum coloring: ", e$message)
            cat("Using default KEGG coloring instead\n")
        })
    }
    
    # Prepare nodes with KEGG styling
    vis_nodes <- prepare_kegg_nodes(nodes, show_labels)
    vis_edges <- prepare_kegg_edges(edges, show_edges)
    

    
        # Create network with KEGG layout
    network <- visNetwork(vis_nodes, vis_edges) %>%
        visNodes(
            font = list(size = 11, strokeWidth = 0),  # Remove text stroke that might interfere
            borderWidth = 1,
            shadow = FALSE  # Disable shadows that might affect text rendering
        ) %>%
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
                values = vis_nodes$id[nodes$type == "gene" | is.na(nodes$type)]
            )
        ) %>%
        visEvents(select = "function(nodes) {
            Shiny.onInputChange('pathway_network_selected', nodes.nodes);
        }")
    
    # Apply KEGG coordinates if available
    if (!is.null(nodes$x) && !is.null(nodes$y) && any(!is.na(nodes$x))) {
        # Use fixed positions - this maintains KEGG layout
        network <- network %>%
            visNodes(physics = FALSE, fixed = list(x = TRUE, y = TRUE)) %>%
            visPhysics(enabled = FALSE) %>%
            visLayout(improvedLayout = FALSE, randomSeed = NULL)
            
        # Print coordinate info for debugging if needed (commented out)
        # coord_range_x <- range(vis_nodes$x, na.rm = TRUE)
        # coord_range_y <- range(vis_nodes$y, na.rm = TRUE)
        # cat("KEGG coordinates - X: [", coord_range_x[1], "-", coord_range_x[2], 
        #     "] Y: [", coord_range_y[1], "-", coord_range_y[2], "]\n")
    } else {
        cat("No KEGG coordinates available\n")
        network <- network %>%
            visPhysics(enabled = FALSE) %>%
            visLayout(randomSeed = 123)
    }
    
    return(network)
}

# Prepare nodes with authentic KEGG styling
prepare_kegg_nodes <- function(nodes, show_labels = TRUE) {
    if (is.null(nodes) || nrow(nodes) == 0) {
        return(data.frame(
            id = character(0),
            label = character(0),
            x = numeric(0),
            y = numeric(0),
            color = character(0),
            size = numeric(0),
            shape = character(0),
            stringsAsFactors = FALSE
        ))
    }
    
    # Start with basic node structure
    vis_nodes <- data.frame(
        id = nodes$id,
        label = if(show_labels) nodes$label else "",
        color = "#BFFFBF",  # Default KEGG gene color
        size = 25,
        shape = "box",      # Default KEGG shape
        borderWidth = 1,
        borderWidthSelected = 2,
        stringsAsFactors = FALSE
    )
    
    # Initialize font as a list column - each row gets a named list
    vis_nodes$font <- replicate(nrow(vis_nodes), list(color = "#000000", size = 11, face = "arial"), simplify = FALSE)
    
    # Apply KEGG coordinates if available
    if (!is.null(nodes$x) && !is.null(nodes$y)) {
        vis_nodes$x <- nodes$x
        vis_nodes$y <- nodes$y
    }
    
    # Add title (tooltip) if available
    if (!is.null(nodes$title)) {
        vis_nodes$title <- nodes$title
    }
    
    # Apply KEGG styling if available
    for (i in seq_len(nrow(vis_nodes))) {
        # Set font color if specified (from phylostratum coloring)
        if (!is.null(nodes$font.color) && !is.na(nodes$font.color[i])) {
            vis_nodes$font[[i]]$color <- nodes$font.color[i]
        }
        
        # Apply phylostratum color first (highest priority)
        if (!is.null(nodes$color.background) && !is.na(nodes$color.background[i])) {
            vis_nodes$color[i] <- nodes$color.background[i]
        } 
        # Apply KEGG background color if no phylostratum color
        else if (!is.null(nodes$kegg_bgcolor) && !is.na(nodes$kegg_bgcolor[i])) {
            vis_nodes$color[i] <- nodes$kegg_bgcolor[i]
        }
        
        # Apply border width if border color is specified (but keep color simple)
        if (!is.null(nodes$color.border) && !is.na(nodes$color.border[i])) {
            vis_nodes$borderWidth[i] <- 2
        }
        
        # Apply KEGG shape
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
        
        # Apply KEGG size based on width
        if (!is.null(nodes$kegg_width) && !is.na(nodes$kegg_width[i])) {
            # Scale KEGG width to reasonable visNetwork size
            kegg_width <- nodes$kegg_width[i]
            vis_nodes$size[i] <- max(20, min(50, kegg_width * 0.8))
        }
        
        # Set font color and size from KEGG foreground color (only if no phylostratum color is set)
        if (!is.null(nodes$kegg_fgcolor) && !is.na(nodes$kegg_fgcolor[i])) {
            # Only apply KEGG foreground color if phylostratum coloring hasn't set a font color
            if (is.null(nodes$font.color) || is.na(nodes$font.color[i])) {
                vis_nodes$font[[i]]$color <- nodes$kegg_fgcolor[i]
                vis_nodes$font[[i]]$size <- 10  # Smaller size for KEGG foreground color
            }
        }
    }
    
    return(vis_nodes)
}

# Prepare edges with KEGG styling and relationship types
prepare_kegg_edges <- function(edges, show_edges = TRUE) {
    if (is.null(edges) || nrow(edges) == 0 || !show_edges) {
        return(data.frame(
            from = character(0),
            to = character(0),
            stringsAsFactors = FALSE
        ))
    }
    
    # Start with basic edge structure
    vis_edges <- data.frame(
        from = edges$from,
        to = edges$to,
        color = "#666666",
        width = 2,
        arrows = list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "arrow")),  # Smaller arrows
        smooth = list(enabled = TRUE, type = "continuous"),
        label = "",  # Add label column for phosphorylation
        font = list(size = 12, align = "middle"),  # Font for edge labels
        stringsAsFactors = FALSE
    )
    
    # Apply KEGG relationship styling if available
    if (!is.null(edges$relation_type) || !is.null(edges$subtype)) {
        for (i in seq_len(nrow(vis_edges))) {
            edge_color <- "#666666"
            edge_width <- 2
            edge_style <- "solid"
            arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "arrow"))  # Default smaller arrow
            edge_label <- ""
            
            # Style based on relation type
            if (!is.null(edges$relation_type) && !is.na(edges$relation_type[i])) {
                relation_type <- edges$relation_type[i]
                switch(relation_type,
                    "PPrel" = {  # Protein-Protein relation
                        edge_color <- "#2E8B57"  # Sea green
                        edge_width <- 2
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "arrow"))
                    },
                    "PCrel" = {  # Protein-Compound relation
                        edge_color <- "#FF6347"  # Tomato red
                        edge_width <- 2
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "arrow"))
                    },
                    "ECrel" = {  # Enzyme-Compound relation
                        edge_color <- "#4169E1"  # Royal blue
                        edge_width <- 2
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "arrow"))
                    },
                    "GErel" = {  # Gene Expression relation
                        edge_color <- "#8A2BE2"  # Blue violet
                        edge_width <- 2
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "arrow"))
                    },
                    "maplink" = {  # Map link
                        edge_color <- "#FFD700"  # Gold
                        edge_width <- 1
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.4, type = "arrow"))
                    },
                    {  # Default
                        edge_color <- "#666666"
                        edge_width <- 2
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "arrow"))
                    }
                )
            }
            
                        # Further style based on subtype
            if (!is.null(edges$subtype) && !is.na(edges$subtype[i]) && edges$subtype[i] != "") {
                subtype <- edges$subtype[i]
                
                # Handle multiple subtypes (separated by ;)
                subtypes <- trimws(strsplit(subtype, ";")[[1]])
                primary_subtype <- subtypes[1]
                
                # Check if phosphorylation is present in any subtype for labeling
                has_phosphorylation <- any(grepl("phosphorylation", subtypes))
                has_dephosphorylation <- any(grepl("dephosphorylation", subtypes))
                
                # Apply primary subtype styling
                switch(primary_subtype,
                    "activation" = {
                        edge_color <- "#00AA00"  # Green for activation
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "arrow"))
                        edge_width <- 3
                    },
                    "inhibition" = {
                        edge_color <- "#FF0000"  # Red for inhibition
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.6, type = "bar"))  # T-shaped end for inhibition
                        edge_width <- 3
                    },
                    "expression" = {
                        edge_color <- "#0066CC"  # Blue for expression
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "arrow"))
                        edge_width <- 2
                    },
                    "repression" = {
                        edge_color <- "#CC0066"  # Dark pink for repression
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.6, type = "bar"))  # T-shaped end for repression
                        edge_width <- 3
                    },
                    "indirect effect" = {
                        edge_color <- "#999999"  # Gray for indirect
                        edge_width <- 1
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.4, type = "arrow"))
                    },
                    "state change" = {
                        edge_color <- "#FF8C00"  # Dark orange
                        edge_width <- 2
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "circle"))
                    },
                    "binding/association" = {
                        edge_color <- "#9932CC"  # Dark orchid
                        arrow_style <- list(to = list(enabled = FALSE))  # No arrow for binding
                        edge_width <- 2
                    },
                    "dissociation" = {
                        edge_color <- "#B22222"  # Fire brick
                        arrow_style <- list(to = list(enabled = FALSE))  # No arrow for dissociation
                        edge_width <- 2
                    },
                    "phosphorylation" = {
                        edge_color <- "#1E90FF"  # Dodger blue
                        edge_width <- 2
                        edge_label <- "+p"  # Add phosphorylation label
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "arrow"))
                    },
                    "dephosphorylation" = {
                        edge_color <- "#FF1493"  # Deep pink
                        edge_width <- 2
                        edge_label <- "-p"  # Add dephosphorylation label
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "arrow"))
                    },
                    "glycosylation" = {
                        edge_color <- "#32CD32"  # Lime green
                        edge_width <- 2
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "arrow"))
                    },
                    "ubiquitination" = {
                        edge_color <- "#DAA520"  # Goldenrod
                        edge_width <- 2
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "arrow"))
                    },
                    "methylation" = {
                        edge_color <- "#20B2AA"  # Light sea green
                        edge_width <- 2
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "arrow"))
                    }
                )
                
                # Add phosphorylation label if present (even as secondary subtype)
                if (has_phosphorylation) {
                    edge_label <- "+p"
                } else if (has_dephosphorylation) {
                    edge_label <- "-p"  
                }
            }
            
            # Apply the styling
            vis_edges$color[i] <- edge_color
            vis_edges$width[i] <- edge_width
            vis_edges$arrows[[i]] <- arrow_style
            vis_edges$label[i] <- edge_label
            
            # Add title for hover information
            if (!is.null(edges$relation_type) && !is.null(edges$subtype)) {
                if (!is.na(edges$relation_type[i]) && !is.na(edges$subtype[i])) {
                    vis_edges$title[i] <- paste0(
                        "Type: ", edges$relation_type[i], "\n",
                        "Subtype: ", edges$subtype[i]
                    )
                }
            }
        }
    }
    
    return(vis_edges)
}

# Create edge legend showing different relationship types
create_edge_legend <- function(edges) {
    if (is.null(edges) || nrow(edges) == 0) {
        return(NULL)
    }
    
    legend_data <- data.frame(
        relationship = character(0),
        color = character(0),
        description = character(0),
        stringsAsFactors = FALSE
    )
    
    # Define relationship types and their styling
    relation_types <- list(
        list(type = "activation", color = "#00AA00", desc = "Activation/stimulation (→)"),
        list(type = "inhibition", color = "#FF0000", desc = "Inhibition/repression (T)"),
        list(type = "expression", color = "#0066CC", desc = "Gene expression"),
        list(type = "repression", color = "#CC0066", desc = "Expression repression (T)"),
        list(type = "phosphorylation", color = "#1E90FF", desc = "Phosphorylation (+p)"),
        list(type = "dephosphorylation", color = "#FF1493", desc = "Dephosphorylation (-p)"),
        list(type = "binding/association", color = "#9932CC", desc = "Binding/association"),
        list(type = "dissociation", color = "#B22222", desc = "Dissociation"),
        list(type = "PPrel", color = "#2E8B57", desc = "Protein-protein relation"),
        list(type = "PCrel", color = "#FF6347", desc = "Protein-compound relation"),
        list(type = "ECrel", color = "#4169E1", desc = "Enzyme-compound relation"),
        list(type = "GErel", color = "#8A2BE2", desc = "Gene expression relation")
    )
    
    # Check which types are present in the data
    present_types <- c()
    if (!is.null(edges$relation_type)) {
        present_types <- c(present_types, unique(edges$relation_type[!is.na(edges$relation_type)]))
    }
    if (!is.null(edges$subtype)) {
        subtypes <- unique(edges$subtype[!is.na(edges$subtype) & edges$subtype != ""])
        # Extract individual subtypes (may be semicolon-separated)
        all_subtypes <- c()
        for (st in subtypes) {
            individual <- trimws(strsplit(st, ";")[[1]])
            all_subtypes <- c(all_subtypes, individual)
        }
        present_types <- c(present_types, unique(all_subtypes))
    }
    
    # Build legend for present types
    for (rel_type in relation_types) {
        if (rel_type$type %in% present_types) {
            legend_data <- rbind(legend_data, data.frame(
                relationship = rel_type$type,
                color = rel_type$color,
                description = rel_type$desc,
                stringsAsFactors = FALSE
            ))
        }
    }
    
    if (nrow(legend_data) > 0) {
        cat("Created edge legend with", nrow(legend_data), "relationship types\n")
        return(legend_data)
    } else {
        return(NULL)
    }
}

#' Load phylostratum mapping data
#' @return data.frame with GeneID and Stratum columns
load_phylomap <- function() {
    phylomap_path <- file.path("data", "phylomap_hgnc.tsv")
    if (!file.exists(phylomap_path)) {
        stop("Phylomap file not found: ", phylomap_path)
    }
    
    phylomap <- read.table(phylomap_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    return(phylomap)
}

#' Map gene symbols to phylostrata
#' @param gene_symbols character vector of gene symbols
#' @param phylomap data.frame with GeneID and Stratum columns
#' @return named vector of strata (names are gene symbols)
map_genes_to_phylostrata <- function(gene_symbols, phylomap) {
    # Create a mapping from gene symbols to strata
    strata_map <- setNames(phylomap$Stratum, phylomap$GeneID)
    
    # Map the gene symbols
    gene_strata <- strata_map[gene_symbols]
    names(gene_strata) <- gene_symbols
    
    return(gene_strata)
}

#' Apply phylostratum coloring to nodes
#' @param nodes data.frame with node information
#' @param phylomap data.frame with phylostratum mapping
#' @return nodes data.frame with phylostratum colors applied
apply_phylostratum_coloring <- function(nodes, phylomap = NULL) {
    if (is.null(phylomap)) {
        phylomap <- load_phylomap()
    }
    
    # Check if myTAI is available
    if (!requireNamespace("myTAI", quietly = TRUE)) {
        stop("myTAI package is required for phylostratum coloring. Please install it with: install.packages('myTAI')")
    }
    
    # Extract gene symbols from node labels
    gene_symbols <- nodes$label
    
    # Map genes to phylostrata
    gene_strata <- map_genes_to_phylostrata(gene_symbols, phylomap)
    
    # Get unique strata present in the data
    unique_strata <- unique(gene_strata[!is.na(gene_strata)])
    
    if (length(unique_strata) == 0) {
        warning("No phylostratum matches found for genes")
        return(nodes)
    }
    
    # Get colors using myTAI PS_colours function
    all_colors <- myTAI::PS_colours(max(unique_strata, na.rm = TRUE))
    
    # Create color mapping
    color_map <- setNames(all_colors[unique_strata], unique_strata)
    
    # Apply colors to nodes
    nodes$color.background <- ifelse(
        is.na(gene_strata), 
        "#D3D3D3",  # Gray for unmapped genes
        color_map[as.character(gene_strata)]
    )
    
    # Set border color
    nodes$color.border <- "#666666"
    
    # Set text color based on background brightness
    nodes$font.color <- ifelse(
        is.na(gene_strata),
        "#000000",  # Black text for gray
        sapply(color_map[as.character(gene_strata)], function(color) {
            if (is.na(color)) return("#000000")
            # Convert hex to RGB and calculate brightness
            rgb_vals <- col2rgb(color)
            brightness <- (rgb_vals[1] * 0.299 + rgb_vals[2] * 0.587 + rgb_vals[3] * 0.114) / 255
            if (brightness < 0.5) "#FFFFFF" else "#000000"  # White text for dark colors
        })
    )
    
    # Add stratum information to title for tooltip
    legend_data <- load_phylostratum_legend()
    if (!is.null(legend_data)) {
        # Create a mapping from rank to name
        rank_to_name <- setNames(legend_data$Name, legend_data$Rank)
        
        nodes$title <- ifelse(
            is.na(gene_strata),
            paste0(nodes$label, "\nPhylostratum: Unknown"),
            paste0(nodes$label, "\nPhylostratum ", gene_strata, ": ", 
                   rank_to_name[as.character(gene_strata)])
        )
    } else {
        nodes$title <- ifelse(
            is.na(gene_strata),
            paste0(nodes$label, "\nPhylostratum: Unknown"),
            paste0(nodes$label, "\nPhylostratum: ", gene_strata)
        )
    }
    
    cat("Applied phylostratum coloring to", sum(!is.na(gene_strata)), "out of", length(gene_strata), "genes\n")
    cat("Strata represented:", paste(sort(unique_strata), collapse = ", "), "\n")
    
    return(nodes)
}

# Load phylostratum legend data
load_phylostratum_legend <- function() {
  legend_file <- file.path("data", "strata_legend.tsv")
  if (!file.exists(legend_file)) {
    warning("Phylostratum legend file not found: ", legend_file)
    return(NULL)
  }
  
  legend_data <- read.table(legend_file, header = TRUE, sep = "\t", 
                           stringsAsFactors = FALSE, quote = '"')
  return(legend_data)
}

# Generate phylostratum legend with colors and names
generate_phylostratum_legend <- function() {
  # Load legend data
  legend_data <- load_phylostratum_legend()
  if (is.null(legend_data)) {
    return(NULL)
  }
  
  # Get colors from myTAI
  if (!requireNamespace("myTAI", quietly = TRUE)) {
    warning("myTAI package required for phylostratum coloring")
    return(NULL)
  }
  
  # Get unique strata (1 to max rank)
  max_stratum <- max(legend_data$Rank)
  all_colors <- myTAI::PS_colours(max_stratum)
  
  # Create legend data frame
  legend_df <- data.frame(
    Rank = legend_data$Rank,
    Name = legend_data$Name,
    Color = all_colors[legend_data$Rank],
    stringsAsFactors = FALSE
  )
  
  return(legend_df)
}
