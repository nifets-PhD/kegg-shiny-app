# Helper functions for KEGG data processing

#' Get centralized KEGG edge relationship color mapping
#' This ensures consistent colors across the entire application
#' @return named vector of colors for KEGG relationship types
get_kegg_edge_colors <- function() {
    c(
        "activation" = "#00AA00",           # Green
        "inhibition" = "#FF0000",           # Red
        "expression" = "#0066CC",           # Blue
        "repression" = "#CC0066",           # Pink
        "phosphorylation" = "#1E90FF",      # Dodger blue
        "dephosphorylation" = "#FF1493",    # Deep pink
        "binding/association" = "#9932CC",   # Dark violet
        "dissociation" = "#B22222",         # Fire brick
        "PPrel" = "#2E8B57",               # Sea green
        "PCrel" = "#FF6347",               # Tomato
        "ECrel" = "#4169E1",               # Royal blue
        "GErel" = "#8A2BE2",               # Blue violet
        "unknown" = "#9E9E9E"              # Grey
    )
}

#' Determine primary interaction relationship from edge data
#' This ensures consistent relationship determination across network edges and gene profiles
#' @param subtype character, subtype string (may be semicolon-separated)
#' @param relation_type character, relation type string
#' @return character, primary relationship type
determine_interaction_relationship <- function(subtype = NULL, relation_type = NULL) {
    # Handle subtype first (more specific)
    if (!is.null(subtype) && !is.na(subtype) && subtype != "") {
        # Handle semicolon-separated subtypes (e.g., "activation; phosphorylation")
        subtype_parts <- trimws(strsplit(subtype, ";")[[1]])
        # Return the first meaningful subtype
        return(subtype_parts[1])
    } 
    
    # Fallback to relation_type
    if (!is.null(relation_type) && !is.na(relation_type) && relation_type != "") {
        return(relation_type)
    }
    
    # Default
    return("unknown")
}

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
                    entry_link <- xml2::xml_attr(entry, "link")
                    
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
                        
                        # Process different entry types differently
                        node_label <- label_text
                        hgnc_symbol <- NA
                        kegg_id <- entry_id
                        description <- paste("KEGG", entry_type)
                        
                        if (entry_type == "gene") {
                            # For genes, extract first gene symbol from graphics name
                            if (!is.na(label_text)) {
                                first_alias <- strsplit(label_text, ",", fixed = TRUE)[[1]][1]
                                first_alias <- trimws(first_alias)
                                first_alias <- gsub("\\.\\.\\.$", "", first_alias)
                                node_label <- first_alias
                                hgnc_symbol <- first_alias
                            }
                            # Extract first numeric ID from entry name
                            numeric_ids <- regmatches(entry_name, gregexpr("\\d+", entry_name))[[1]]
                            if (length(numeric_ids) > 0) {
                                kegg_id <- numeric_ids[1]
                            }
                            description <- paste("Gene:", node_label)
                        } else if (entry_type == "compound") {
                            # For compounds, use compound ID and name from graphics
                            if (!is.na(label_text)) {
                                node_label <- label_text
                            }
                            description <- paste("Compound:", node_label)
                        } else if (entry_type == "map") {
                            # For pathway maps, extract pathway info
                            if (!is.na(label_text)) {
                                node_label <- label_text
                            }
                            description <- paste("Pathway:", node_label)
                        } else if (entry_type == "group") {
                            # For groups (complexes), use group name
                            if (!is.na(label_text)) {
                                node_label <- label_text
                            }
                            description <- paste("Complex/Group:", node_label)
                        } else if (entry_type == "ortholog") {
                            # For orthologs
                            if (!is.na(label_text)) {
                                node_label <- label_text
                            }
                            description <- paste("Ortholog:", node_label)
                        }
                        
                        cat("  Graphics - X:", x, "Y:", y, "Type:", entry_type, "Label:", node_label, "\n")
                        
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
                                label = node_label,
                                hgnc_symbol = hgnc_symbol,
                                kegg_id = kegg_id,
                                description = description,
                                link = ifelse(is.na(entry_link), "", entry_link)
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
                            hgnc_symbol = x$hgnc_symbol,
                            kegg_id = x$kegg_id,
                            description = x$description,
                            link = x$link,
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

#' Get genes in a specific KEGG pathway
#' @param pathway_id KEGG pathway ID (e.g., "hsa04110")
#' @return character vector of gene symbols found in the pathway
get_pathway_genes <- function(pathway_id) {
    tryCatch({
        # Clean pathway ID
        pathway_id <- gsub("^path:", "", pathway_id)
        
        # Use KEGGREST to get pathway gene list
        if (requireNamespace("KEGGREST", quietly = TRUE)) {
            # Get genes associated with this pathway
            pathway_genes <- KEGGREST::keggGet(pathway_id)
            
            if (length(pathway_genes) > 0 && !is.null(pathway_genes[[1]]$GENE)) {
                # Extract gene information
                gene_info <- pathway_genes[[1]]$GENE
                
                # Gene info comes as named vector where names are KEGG IDs and values are "Symbol Description"
                gene_symbols <- c()
                
                if (is.character(gene_info)) {
                    # Parse gene symbols from the descriptions
                    for (gene_entry in gene_info) {
                        # Split by semicolon and take first part, then extract first word
                        gene_parts <- unlist(strsplit(gene_entry, ";"))
                        if (length(gene_parts) > 0) {
                            # Extract the gene symbol (first word before any description)
                            symbol_part <- trimws(gene_parts[1])
                            gene_symbol <- unlist(strsplit(symbol_part, "\\s+"))[1]
                            if (!is.na(gene_symbol) && gene_symbol != "") {
                                gene_symbols <- c(gene_symbols, gene_symbol)
                            }
                        }
                    }
                } else if (is.vector(gene_info) && !is.null(names(gene_info))) {
                    # Alternative format: named vector where names are KEGG IDs
                    for (gene_entry in gene_info) {
                        # Extract gene symbol (first word before any description)
                        gene_symbol <- unlist(strsplit(trimws(gene_entry), "\\s+"))[1]
                        if (!is.na(gene_symbol) && gene_symbol != "") {
                            gene_symbols <- c(gene_symbols, gene_symbol)
                        }
                    }
                }
                
                return(unique(toupper(gene_symbols)))
            }
        } else {
            # Fallback: try to extract from pathway names (previous method)
            cat("KEGGREST not available, using basic search for", pathway_id, "\n")
            return(character(0))
        }
        
        return(character(0))
        
    }, error = function(e) {
        cat("Error getting genes for pathway", pathway_id, ":", e$message, "\n")
        return(character(0))
    })
}

filter_pathways_by_category <- function(pathways_list, category) {
    if (is.null(pathways_list)) {
        return(NULL)
    }
    
    # Extract the first 2 digits from category for prefix matching
    category_prefix <- substr(category, 1, 2)
    
    # Filter pathways based on category prefix (e.g., "01" for metabolism)
    category_matches <- grepl(paste0("^hsa", category_prefix), pathways_list$pathway_id)
    
    result <- pathways_list[category_matches, ]
    cat("Category", category, "(prefix", category_prefix, ") matches:", nrow(result), "pathways\n")
    
    return(result)
}

load_kegg_pathway <- function(pathway_id) {
    tryCatch({
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
            cat("Sample node IDs:", paste(head(nodes_data$id, 3), collapse = ", "), "\n")
            
            # Also get edges from KEGG XML (like the old version)
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
create_kegg_network_visualization <- function(nodes, edges, show_labels = TRUE, show_edges = TRUE, coloring_mode = "kegg_default", highlight_genes = NULL, id_type = "symbol") {
    
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
            nodes <- apply_phylostratum_coloring(nodes, id_type = id_type)
        }, error = function(e) {
            warning("Failed to apply phylostratum coloring: ", e$message)
            cat("Using default KEGG coloring instead\n")
        })
    }
    
    # Apply gene highlighting if requested - do this AFTER KEGG styling to override colors
    if (!is.null(highlight_genes) && length(highlight_genes) > 0) {
        nodes <- apply_gene_highlighting(nodes, highlight_genes)
    }
    
    # Prepare nodes with KEGG styling
    vis_nodes <- prepare_kegg_nodes(nodes, show_labels, NULL)  # Don't pass highlight_genes again
    vis_edges <- prepare_kegg_edges(edges, show_edges)
    

    
        # Create network with KEGG layout
    network <- visNetwork(vis_nodes, vis_edges) %>%
        visNodes(
            font = list(
                size = 14,  # Increased base font size
                strokeWidth = 1,
                strokeColor = "white"  # White outline for better contrast
            ),
            borderWidth = 1,
            shadow = FALSE,  # Disable shadows that might affect text rendering
            scaling = list(
                min = 10,
                max = 30,
                label = list(
                    enabled = TRUE,
                    min = 14,
                    max = 24,
                    maxVisible = 30,
                    drawThreshold = 1  # Very low threshold to always show labels
                )
            )
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
        network <- network %>%
            visPhysics(enabled = FALSE) %>%
            visLayout(randomSeed = 123)
    }
    
    return(network)
}

# Prepare nodes with authentic KEGG styling
prepare_kegg_nodes <- function(nodes, show_labels = TRUE, highlight_genes = NULL) {
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
        hgnc_symbol = if(!is.null(nodes$hgnc_symbol)) nodes$hgnc_symbol else nodes$label,
        hidden = FALSE,     # Default to visible
        stringsAsFactors = FALSE
    )
    
    # Hide group nodes
    if (!is.null(nodes$type)) {
        vis_nodes$hidden[nodes$type == "group"] <- TRUE
    }
    
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
        # Apply direct color property (for gene highlighting)
        else if (!is.null(nodes$color) && !is.na(nodes$color[i]) && 
                 is.character(nodes$color) && nchar(nodes$color[i]) > 0) {
            vis_nodes$color[i] <- nodes$color[i]
        }
        # Apply KEGG background color if no other color is set
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
    
    # Apply gene highlighting if requested
    if (!is.null(highlight_genes) && length(highlight_genes) > 0) {
        vis_nodes <- apply_gene_highlighting(vis_nodes, highlight_genes)
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
                # Use centralized function to determine primary relationship
                primary_relationship <- determine_interaction_relationship(
                    subtype = edges$subtype[i],
                    relation_type = edges$relation_type[i]
                )
                
                # Get centralized edge colors
                edge_colors <- get_kegg_edge_colors()
                
                # Apply styling based on primary relationship
                if (primary_relationship %in% names(edge_colors)) {
                    edge_color <- edge_colors[primary_relationship]
                }
                
                # Apply specific styling rules for different relationships
                switch(primary_relationship,
                    "activation" = {
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "arrow"))
                        edge_width <- 3
                    },
                    "inhibition" = {
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.6, type = "bar"))  # T-shaped end for inhibition
                        edge_width <- 3
                    },
                    "expression" = {
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "arrow"))
                        edge_width <- 2
                    },
                    "repression" = {
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.6, type = "bar"))  # T-shaped end for repression
                        edge_width <- 3
                    },
                    "phosphorylation" = {
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "arrow"))
                        edge_width <- 2
                    },
                    "dephosphorylation" = {
                        arrow_style <- list(to = list(enabled = TRUE, scaleFactor = 0.5, type = "arrow"))
                        edge_width <- 2
                    }
                )
                
                # Handle phosphorylation labeling (check if any subtype contains phosphorylation)
                if (!is.null(edges$subtype[i])) {
                    subtypes <- trimws(strsplit(edges$subtype[i], ";")[[1]])
                    has_phosphorylation <- any(grepl("phosphorylation", subtypes))
                    has_dephosphorylation <- any(grepl("dephosphorylation", subtypes))
                    
                    if (has_phosphorylation) {
                        edge_label <- "+p"
                    } else if (has_dephosphorylation) {
                        edge_label <- "-p"
                    }
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
    
    # Define relationship types and their styling using centralized colors
    edge_colors <- get_kegg_edge_colors()
    relation_types <- list(
        list(type = "activation", color = edge_colors["activation"], desc = "Activation/stimulation (→)"),
        list(type = "inhibition", color = edge_colors["inhibition"], desc = "Inhibition/repression (T)"),
        list(type = "expression", color = edge_colors["expression"], desc = "Gene expression"),
        list(type = "repression", color = edge_colors["repression"], desc = "Expression repression (T)"),
        list(type = "phosphorylation", color = edge_colors["phosphorylation"], desc = "Phosphorylation (+p)"),
        list(type = "dephosphorylation", color = edge_colors["dephosphorylation"], desc = "Dephosphorylation (-p)"),
        list(type = "binding/association", color = edge_colors["binding/association"], desc = "Binding/association"),
        list(type = "dissociation", color = edge_colors["dissociation"], desc = "Dissociation"),
        list(type = "PPrel", color = edge_colors["PPrel"], desc = "Protein-protein relation"),
        list(type = "PCrel", color = edge_colors["PCrel"], desc = "Protein-compound relation"),
        list(type = "ECrel", color = edge_colors["ECrel"], desc = "Enzyme-compound relation"),
        list(type = "GErel", color = edge_colors["GErel"], desc = "Gene expression relation")
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
    # Use comprehensive phylomap with Ensembl data
    phylomap_path <- file.path("data", "phylomap.tsv")
    if (!file.exists(phylomap_path)) {
        stop("Comprehensive phylomap file not found: ", phylomap_path)
    }
    
    phylomap <- read.table(phylomap_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    cat("Loaded comprehensive phylomap with", nrow(phylomap), "protein entries\n")
    return(phylomap)
}

#' Map gene symbols to phylostrata
#' @param gene_symbols character vector of gene symbols
#' @param phylomap data.frame with GeneID and Stratum columns
#' @return named vector of strata (names are gene symbols)
map_genes_to_phylostrata <- function(gene_ids, id_type = "symbol", phylomap = NULL) {
    # Enhanced function to support multiple gene ID types
    
    # Load phylomap if not provided
    if (is.null(phylomap)) {
        phylomap <- load_phylomap()
    }
    
    if (id_type == "symbol") {
        # Use HGNC column for gene symbols with duplicate handling
        # Filter out empty HGNC symbols
        valid_hgnc <- phylomap[phylomap$HGNC != "" & !is.na(phylomap$HGNC), ]
        
        # Handle duplicates by taking the minimum stratum (most ancient/conservative)
        hgnc_strata <- aggregate(Stratum ~ HGNC, data = valid_hgnc, FUN = min)
        strata_map <- setNames(hgnc_strata$Stratum, hgnc_strata$HGNC)
        gene_strata <- strata_map[gene_ids]
        
        # Try case-insensitive matching for unmapped genes
        unmapped <- is.na(gene_strata)
        if (any(unmapped)) {
            # Create case-insensitive mapping
            upper_hgnc_strata <- aggregate(Stratum ~ HGNC, 
                                         data = transform(valid_hgnc, HGNC = toupper(HGNC)), 
                                         FUN = min)
            upper_strata_map <- setNames(upper_hgnc_strata$Stratum, upper_hgnc_strata$HGNC)
            gene_strata[unmapped] <- upper_strata_map[toupper(gene_ids[unmapped])]
        }
        
    } else if (id_type == "entrez") {
        # Use Entrez_ID column for Entrez IDs with duplicate handling
        # Filter out empty/NA Entrez IDs
        valid_entrez <- phylomap[!is.na(phylomap$Entrez_ID) & phylomap$Entrez_ID != "", ]
        
        # Handle duplicates by taking the minimum stratum (most ancient/conservative)
        entrez_strata <- aggregate(Stratum ~ Entrez_ID, data = valid_entrez, FUN = min)
        entrez_strata_map <- setNames(entrez_strata$Stratum, as.character(entrez_strata$Entrez_ID))
        gene_strata <- entrez_strata_map[as.character(gene_ids)]
        
    } else if (id_type == "uniprot") {
        # Use UniProt column for UniProt IDs with duplicate handling
        # Filter out empty UniProt IDs
        valid_uniprot <- phylomap[phylomap$UniProt != "" & !is.na(phylomap$UniProt), ]
        
        # Handle duplicates by taking the minimum stratum (most ancient/conservative)
        uniprot_strata <- aggregate(Stratum ~ UniProt, data = valid_uniprot, FUN = min)
        uniprot_strata_map <- setNames(uniprot_strata$Stratum, uniprot_strata$UniProt)
        gene_strata <- uniprot_strata_map[gene_ids]
        
    } else if (id_type == "ensembl") {
        # Use ENSG column for Ensembl Gene IDs with duplicate handling
        # Filter out empty ENSG IDs
        valid_ensembl <- phylomap[phylomap$ENSG != "" & !is.na(phylomap$ENSG), ]
        
        # Handle duplicates by taking the minimum stratum (most ancient/conservative)
        ensembl_strata <- aggregate(Stratum ~ ENSG, data = valid_ensembl, FUN = min)
        ensembl_strata_map <- setNames(ensembl_strata$Stratum, ensembl_strata$ENSG)
        gene_strata <- ensembl_strata_map[gene_ids]
        
    } else {
        # Fallback to HGNC symbols
        warning("Unknown id_type '", id_type, "', falling back to HGNC symbol mapping")
        strata_map <- setNames(phylomap$Stratum, phylomap$HGNC)
        gene_strata <- strata_map[gene_ids]
    }
    
    names(gene_strata) <- gene_ids
    
    # Report mapping statistics
    mapped_count <- sum(!is.na(gene_strata))
    total_count <- length(gene_ids)
    cat("Phylostratum mapping (", id_type, "): ", mapped_count, "/", total_count, 
        " (", round(100 * mapped_count/total_count, 1), "%)\n", sep="")
    
    return(gene_strata)
}

# Phylostratum color generation function (from MyTAI dev version)
# This function uses the pre-computed global color palette for consistency
PS_colours <- function(n) {
    # Use the global color palette defined in global.R
    return(GLOBAL_PHYLOSTRATA_COLORS[1:min(n, length(GLOBAL_PHYLOSTRATA_COLORS))])
}

#' Apply phylostratum coloring to nodes
#' @param nodes data.frame with node information
#' @param phylomap data.frame with phylostratum mapping
#' @return nodes data.frame with phylostratum colors applied
apply_phylostratum_coloring <- function(nodes, phylomap = NULL, id_type = "symbol") {
    if (is.null(phylomap)) {
        phylomap <- load_phylomap()
    }
    
    # Filter nodes to only apply phylostratum coloring to gene nodes
    # Identify gene nodes based on the type column
    is_gene_node <- rep(FALSE, nrow(nodes))
    if (!is.null(nodes$type)) {
        is_gene_node <- nodes$type == "gene" | is.na(nodes$type)
    } else {
        # If no type column, assume all nodes with valid gene IDs are genes
        is_gene_node <- rep(TRUE, nrow(nodes))
    }
    
    # Further filter by checking if nodes have actual gene IDs
    # Exclude nodes that are clearly not genes (like pathway titles)
    if (!is.null(nodes$kegg_id)) {
        # Only consider nodes with valid Entrez IDs as potential genes
        has_valid_entrez <- !is.na(nodes$kegg_id) & nodes$kegg_id != "" & 
                           grepl("^\\d+$", nodes$kegg_id)  # Must be numeric
        is_gene_node <- is_gene_node & has_valid_entrez
    }
    
    # If no gene nodes found, return nodes unchanged
    gene_node_count <- sum(is_gene_node)
    if (gene_node_count == 0) {
        cat("No gene nodes found for phylostratum coloring\n")
        return(nodes)
    }
    
    cat("Applying phylostratum coloring to", gene_node_count, "out of", nrow(nodes), "nodes\n")
    
    # Extract gene IDs only from gene nodes
    gene_nodes <- nodes[is_gene_node, ]
    gene_ids <- NULL
    actual_id_type <- id_type
    
    # For KEGG pathway nodes, prioritize Entrez IDs over other ID types
    if (!is.null(gene_nodes$kegg_id) && any(!is.na(gene_nodes$kegg_id) & gene_nodes$kegg_id != "")) {
        cat("Using Entrez IDs (kegg_id) from gene nodes for phylostratum mapping\n")
        gene_ids <- gene_nodes$kegg_id
        actual_id_type <- "entrez"
    }
    # If no Entrez IDs, fall back to HGNC symbols
    else if (!is.null(gene_nodes$hgnc_symbol) && any(!is.na(gene_nodes$hgnc_symbol) & gene_nodes$hgnc_symbol != "")) {
        cat("Using HGNC symbols from gene nodes for phylostratum mapping\n")
        gene_ids <- gene_nodes$hgnc_symbol
        actual_id_type <- "symbol"
    }
    # If no HGNC symbols, try to use labels as last resort
    else if (!is.null(gene_nodes$label)) {
        cat("Using node labels for phylostratum mapping with id_type:", id_type, "\n")
        gene_ids <- gene_nodes$label
        # Keep the original id_type for labels
    }
    else {
        warning("No suitable gene IDs found in gene nodes for phylostratum mapping")
        return(nodes)
    }
    
    # Map genes to phylostrata using the determined ID type
    gene_strata <- map_genes_to_phylostrata(gene_ids, id_type = actual_id_type, phylomap = phylomap)
    
    # Get unique strata present in the data
    unique_strata <- unique(gene_strata[!is.na(gene_strata)])
    
    if (length(unique_strata) == 0) {
        warning("No phylostratum matches found for genes")
        return(nodes)
    }
    
    # Get colors using our PS_colours function
    all_colors <- PS_colours(max(unique_strata, na.rm = TRUE))
    
    # Create color mapping
    color_map <- setNames(all_colors[unique_strata], unique_strata)
    
    # Initialize color columns for all nodes (only gene nodes will be modified)
    nodes$color.background <- NA
    nodes$color.border <- NA
    nodes$font.color <- NA
    
    # Apply colors only to gene nodes
    gene_indices <- which(is_gene_node)
    
    nodes$color.background[gene_indices] <- ifelse(
        is.na(gene_strata), 
        "#D3D3D3",  # Gray for unmapped genes
        color_map[as.character(gene_strata)]
    )
    
    # Set border color only for gene nodes
    nodes$color.border[gene_indices] <- "#666666"
    
    # Set text color based on background brightness only for gene nodes
    nodes$font.color[gene_indices] <- ifelse(
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
    
    # Add stratum information to title for tooltip - only for gene nodes
    legend_data <- load_phylostratum_legend()
    
    # Initialize title column if it doesn't exist
    if (is.null(nodes$title)) {
        nodes$title <- nodes$label  # Default title is the label
    }
    
    if (!is.null(legend_data)) {
        # Create a mapping from rank to name
        rank_to_name <- setNames(legend_data$Name, legend_data$Rank)
        
        # Only update titles for gene nodes that have phylostratum data
        valid_gene_indices <- gene_indices[!is.na(gene_strata)]
        valid_gene_strata <- gene_strata[!is.na(gene_strata)]
        
        if (length(valid_gene_indices) > 0) {
            nodes$title[valid_gene_indices] <- paste0(
                nodes$label[valid_gene_indices], 
                "\nPhylostratum ", valid_gene_strata, ": ", 
                rank_to_name[as.character(valid_gene_strata)]
            )
        }
        
        # Update titles for gene nodes without phylostratum data
        unknown_gene_indices <- gene_indices[is.na(gene_strata)]
        if (length(unknown_gene_indices) > 0) {
            nodes$title[unknown_gene_indices] <- paste0(
                nodes$label[unknown_gene_indices], 
                "\nPhylostratum: Unknown"
            )
        }
    } else {
        # Only update titles for gene nodes that have phylostratum data
        valid_gene_indices <- gene_indices[!is.na(gene_strata)]
        valid_gene_strata <- gene_strata[!is.na(gene_strata)]
        
        if (length(valid_gene_indices) > 0) {
            nodes$title[valid_gene_indices] <- paste0(
                nodes$label[valid_gene_indices], 
                "\nPhylostratum: ", valid_gene_strata
            )
        }
        
        # Update titles for gene nodes without phylostratum data
        unknown_gene_indices <- gene_indices[is.na(gene_strata)]
        if (length(unknown_gene_indices) > 0) {
            nodes$title[unknown_gene_indices] <- paste0(
                nodes$label[unknown_gene_indices], 
                "\nPhylostratum: Unknown"
            )
        }
    }
    
    cat("Applied phylostratum coloring to", sum(!is.na(gene_strata)), "out of", length(gene_strata), "gene nodes\n")
    cat("Strata represented:", paste(sort(unique_strata), collapse = ", "), "\n")
    cat("Non-gene nodes (titles, compounds, etc.) left unchanged:", nrow(nodes) - gene_node_count, "nodes\n")
    
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
  
  # Get unique strata (1 to max rank)
  max_stratum <- max(legend_data$Rank)
  all_colors <- PS_colours(max_stratum)
  
  # Create legend data frame
  legend_df <- data.frame(
    Rank = legend_data$Rank,
    Name = legend_data$Name,
    Color = all_colors[legend_data$Rank],
    stringsAsFactors = FALSE
  )
  
  return(legend_df)
}

#' Parse gene list from text input
#' @param gene_text character string with gene symbols separated by newlines
#' @details Gene symbols MUST be separated by newlines (one gene per line).
#'   Commas, semicolons, or other separators are NOT supported and will be 
#'   treated as part of the gene name. Each line should contain exactly one gene symbol.
#' @return character vector of cleaned gene symbols
parse_gene_text <- function(gene_text) {
    if (is.null(gene_text) || gene_text == "") {
        return(character(0))
    }
    
    # Check for common separator mistakes and warn user
    if (grepl(",", gene_text) || grepl(";", gene_text) || grepl("\\t", gene_text)) {
        # DEBUG: Show exactly what was found
        found_chars <- c()
        if (grepl(",", gene_text)) found_chars <- c(found_chars, "comma")
        if (grepl(";", gene_text)) found_chars <- c(found_chars, "semicolon") 
        if (grepl("\\t", gene_text)) found_chars <- c(found_chars, "tab")
        
        warning_msg <- paste0("PARSING ERROR: Gene symbols must be separated by newlines only! Found ", 
                             paste(found_chars, collapse=", "), " in input. Please use one gene per line.")
        warning(warning_msg)
    }
    
    # Split by newlines ONLY - this enforces the newline-separated format
    lines <- unlist(strsplit(gene_text, "[\n\r]+"))
    
    # Clean up - remove empty strings and whitespace
    lines <- trimws(lines)
    lines <- lines[lines != ""]
    
    # Convert to uppercase for consistency
    genes <- toupper(lines)
    
    return(unique(genes))
}

#' Load gene list from file
#' @param file_path path to file containing gene symbols
#' @return character vector of gene symbols
load_gene_file <- function(file_path) {
    tryCatch({
        # Check if file exists
        if (!file.exists(file_path)) {
            stop("File does not exist: ", file_path)
        }
        
        # Determine file type
        ext <- tools::file_ext(file_path)
        
        if (ext %in% c("txt", "tsv")) {
            # Read as plain text
            lines <- readLines(file_path, warn = FALSE)
            
            # Clean up lines directly (same logic as parse_gene_text but without comma warning)
            lines <- trimws(lines)
            lines <- lines[lines != ""]
            
            genes <- toupper(lines)
            return(unique(genes))
        } else if (ext == "csv") {
            # Read as CSV and take first column
            data <- read.csv(file_path, stringsAsFactors = FALSE, header = FALSE)
            genes <- as.character(data[, 1])
            # Clean up genes directly (same logic as parse_gene_text but without comma warning)
            genes <- trimws(genes)
            genes <- genes[genes != ""]
            genes <- toupper(genes)
            return(unique(genes))
        } else {
            warning("Unsupported file type: ", ext)
            return(character(0))
        }
    }, error = function(e) {
        warning("Error reading gene file: ", e$message)
        return(character(0))
    })
}

#' Find pathways containing specific genes
#' @param gene_symbols character vector of HGNC gene symbols
#' @param pathways_list data.frame of pathway information
#' @return data.frame of pathways with gene overlap information
find_pathways_with_genes <- function(gene_symbols, pathways_list) {
    if (length(gene_symbols) == 0 || is.null(pathways_list)) {
        return(data.frame())
    }
    
    cat("Searching for pathways containing", length(gene_symbols), "genes...\n")
    
    # Initialize results with gene overlap information
    pathway_results <- pathways_list
    pathway_results$genes_found <- ""
    pathway_results$gene_count <- 0
    pathway_results$priority_score <- 0
    
    # Convert gene symbols to uppercase for matching
    gene_symbols_upper <- toupper(gene_symbols)
    
    # Limit search to prevent server overload (search first 50 pathways)
    max_pathways_to_search <- min(50, nrow(pathway_results))
    cat("Limiting search to first", max_pathways_to_search, "pathways to avoid server overload\n")
    
    # Send initial progress notification
    if (exists("showNotification", mode = "function")) {
        showNotification(
            paste("Starting gene search across", max_pathways_to_search, "pathways..."),
            type = "message",
            duration = 3
        )
    }
    
    # Progress counter
    successful_queries <- 0
    matches_found <- 0
    
    # Query each pathway for its actual gene content
    for (i in seq_len(max_pathways_to_search)) {
        # Progress update every 10 pathways
        if (i %% 10 == 0) {  
            cat("Processing pathway", i, "of", max_pathways_to_search, "...\n")
            
            # Send progress notification to UI
            if (exists("showNotification", mode = "function")) {
                progress_msg <- paste("Processed", i, "of", max_pathways_to_search, "pathways.", 
                                    "Found", matches_found, "matches so far...")
                showNotification(progress_msg, type = "message", duration = 2)
            }
        }
        
        pathway_id <- pathway_results$pathway_id[i]
        pathway_name <- pathway_results$pathway_name[i]
        
        # Get genes in this pathway from KEGG
        pathway_genes <- get_pathway_genes(pathway_id)
        
        if (length(pathway_genes) > 0) {
            successful_queries <- successful_queries + 1
            
            # Find intersection between user genes and pathway genes
            matched_genes <- intersect(gene_symbols_upper, pathway_genes)
            
            if (length(matched_genes) > 0) {
                matches_found <- matches_found + 1
                pathway_results$genes_found[i] <- paste(matched_genes, collapse = ", ")
                pathway_results$gene_count[i] <- length(matched_genes)
                pathway_results$priority_score[i] <- length(matched_genes)
                
                cat("  ✓ Found", length(matched_genes), "matching genes in", pathway_name, "\n")
                
                # Notify about significant matches
                if (length(matched_genes) >= 2) {
                    if (exists("showNotification", mode = "function")) {
                        showNotification(
                            paste("Great match!", length(matched_genes), "genes found in", pathway_name),
                            type = "message",
                            duration = 3
                        )
                    }
                }
            }
        }
        
        # Small delay to avoid overwhelming KEGG servers
        Sys.sleep(0.2)
        
        # Break early if we find some good results to speed up response
        if (successful_queries >= 20 && matches_found >= 5) {
            cat("Found enough results, stopping early to improve speed\n")
            if (exists("showNotification", mode = "function")) {
                showNotification("Found enough good matches, finishing search early!", type = "message", duration = 3)
            }
            break
        }
    }
    
    # Final progress notification
    if (exists("showNotification", mode = "function")) {
        showNotification(
            paste("Search complete! Found", matches_found, "pathways with your genes."),
            type = if (matches_found > 0) "message" else "warning",
            duration = 5
        )
    }
    
    # Sort by priority score (pathways with more matching genes first)
    pathway_results <- pathway_results[order(-pathway_results$priority_score, pathway_results$pathway_name), ]
    
    # Filter to only show pathways with at least one matching gene
    relevant_pathways <- pathway_results[pathway_results$gene_count > 0, ]
    
    if (nrow(relevant_pathways) == 0) {
        cat("No pathways found containing your genes.\n")
        cat("This could mean:\n")
        cat("1. Gene symbols don't match KEGG database (try different gene names)\n")
        cat("2. Genes are not in any human pathways\n") 
        cat("3. Network connectivity issues with KEGG server\n")
        
        # Return a helpful message instead of empty results
        return(data.frame(
            pathway_id = "no_results",
            pathway_name = paste("No pathways found containing:", paste(gene_symbols[1:min(3, length(gene_symbols))], collapse = ", ")),
            description = "Try: 1) Check gene symbols are correct (e.g., TP53, BRCA1), 2) Use official HGNC symbols, 3) Try category search instead",
            genes_found = "",
            gene_count = 0,
            priority_score = 0,
            stringsAsFactors = FALSE
        ))
    }
    
    cat("Found", nrow(relevant_pathways), "pathways containing your genes\n")
    cat("Top matches:", paste(head(relevant_pathways$pathway_name, 3), collapse = ", "), "\n")
    
    return(relevant_pathways)
}

#' Apply gene highlighting to nodes
#' @param nodes data.frame with node information
#' @param highlight_genes character vector of gene identifiers to highlight (can be Entrez IDs, symbols, etc.)
#' @return nodes data.frame with highlighting applied
apply_gene_highlighting <- function(nodes, highlight_genes) {
    if (length(highlight_genes) == 0) {
        return(nodes)
    }
    
    # Convert highlight genes to character for matching
    highlight_genes <- as.character(highlight_genes)
    highlight_genes_upper <- toupper(highlight_genes)
    
    # Check which nodes match the genes of interest
    # Only highlight nodes that are actually genes (not compounds or other types)
    nodes$is_highlighted <- FALSE
    
    # Only consider nodes that are genes (check type column if it exists)
    is_gene_node <- rep(TRUE, nrow(nodes))  # Default: assume all are genes
    if (!is.null(nodes$type)) {
        is_gene_node <- nodes$type == "gene" | is.na(nodes$type)
    }
    
    # Priority 1: Match by Entrez IDs (most reliable for KEGG pathways)
    if (!is.null(nodes$kegg_id)) {
        # Check if highlight genes are numeric (Entrez IDs)
        numeric_highlights <- highlight_genes[grepl("^\\d+$", highlight_genes)]
        if (length(numeric_highlights) > 0) {
            kegg_matches <- is_gene_node & (nodes$kegg_id %in% numeric_highlights)
            nodes$is_highlighted <- nodes$is_highlighted | kegg_matches
            cat("Matched", sum(kegg_matches, na.rm = TRUE), "genes by Entrez ID (kegg_id)\n")
        }
    }
    
    # Also check gene_name field for Entrez IDs (like "hsa:2475 hsa:57521")
    if (!is.null(nodes$gene_name)) {
        numeric_highlights <- highlight_genes[grepl("^\\d+$", highlight_genes)]
        if (length(numeric_highlights) > 0) {
            for (i in seq_len(nrow(nodes))) {
                if (is_gene_node[i] && !is.na(nodes$gene_name[i])) {
                    # Extract numeric IDs from entry names like "hsa:2475 hsa:57521"
                    node_entrez_ids <- regmatches(nodes$gene_name[i], gregexpr("\\d+", nodes$gene_name[i]))[[1]]
                    if (any(numeric_highlights %in% node_entrez_ids)) {
                        nodes$is_highlighted[i] <- TRUE
                    }
                }
            }
        }
    }
    
    # Priority 2: Match by HGNC symbols (for backward compatibility)
    if (!is.null(nodes$hgnc_symbol)) {
        hgnc_matches <- is_gene_node & (toupper(nodes$hgnc_symbol) %in% highlight_genes_upper)
        nodes$is_highlighted <- nodes$is_highlighted | hgnc_matches
        cat("Matched", sum(hgnc_matches, na.rm = TRUE), "genes by HGNC symbol\n")
    }
    
    # Priority 3: Match by labels (fallback)
    if (!is.null(nodes$label)) {
        label_matches <- is_gene_node & (toupper(nodes$label) %in% highlight_genes_upper)
        nodes$is_highlighted <- nodes$is_highlighted | label_matches
        cat("Matched", sum(label_matches, na.rm = TRUE), "genes by node label\n")
    }
    
    # Apply highlighting colors - only modify highlighted nodes
    # For highlighted nodes, set red color
    nodes$color[nodes$is_highlighted] <- "#FF6B6B"  # Bright red for highlighted genes
    
    # Set borders for highlighted nodes
    nodes$color.border[nodes$is_highlighted] <- "#FF0000"  # Red border for highlighted genes
    
    # Make highlighted genes larger
    if (!is.null(nodes$size)) {
        nodes$size[nodes$is_highlighted] <- 35
    }
    
    # Set font colors for highlighted nodes only
    nodes$font.color[nodes$is_highlighted] <- "#FFFFFF"
    
    # Add tooltip information
    highlighted_count <- sum(nodes$is_highlighted)
    nodes$title <- ifelse(
        nodes$is_highlighted,
        paste0(nodes$label, "\n⭐ GENE OF INTEREST ⭐\n", 
               "One of ", highlighted_count, " genes you're tracking"),
        if (!is.null(nodes$title)) nodes$title else nodes$label
    )
    
    cat("Applied highlighting to", highlighted_count, "out of", nrow(nodes), "nodes\n")
    
    return(nodes)
}

#' Parse KEGG pathway (simplified without HSA data)
#' @param pathway_id KEGG pathway ID (e.g., "hsa04152")
#' @param use_cached logical, whether to use cached KEGG data
#' @return list containing nodes, edges, and kegg_data
parse_kegg_pathway_with_hsa <- function(pathway_id, use_cached = TRUE) {
    cat("Loading KEGG pathway", pathway_id, "...\n")
    
    # Load KEGG pathway data
    kegg_result <- load_kegg_pathway(pathway_id)
    
    if (is.null(kegg_result$nodes) || nrow(kegg_result$nodes) == 0) {
        cat("No pathway data loaded\n")
        return(list(nodes = data.frame(), edges = data.frame()))
    }
    
    # The kegg_data should be the parsed XML from load_kegg_pathway
    kegg_data <- kegg_result$graph  # Try using 'graph' instead of 'kegg_data'
    
    return(list(
        nodes = kegg_result$nodes,
        edges = kegg_result$edges,
        kegg_data = kegg_data
    ))
}
