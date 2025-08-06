# Evolutionary Transcriptomics Utility Functions
# Functions for integrating myTAI with the KEGG Shiny app

#' Validate expression data format
#' @param data data.frame with expression data
#' @param gene_id_type character, type of gene IDs ("symbol", "entrez", "ensembl", "uniprot")
#' @return list with validation results
validate_expression_data <- function(data, gene_id_type = "symbol") {
    results <- list(
        is_valid = FALSE,
        errors = character(0),
        warnings = character(0),
        n_genes = 0,
        n_samples = 0,
        gene_column = character(0),
        sample_columns = character(0)
    )
    
    # Check if data exists
    if (is.null(data) || nrow(data) == 0) {
        results$errors <- c(results$errors, "No data provided or data is empty")
        return(results)
    }
    
    # Check minimum dimensions
    if (ncol(data) < 2) {
        results$errors <- c(results$errors, "Data must have at least 2 columns (gene IDs + 1 sample)")
        return(results)
    }
    
    # Identify gene column (should be first column)
    gene_column <- names(data)[1]
    sample_columns <- names(data)[-1]
    
    results$gene_column <- gene_column
    results$sample_columns <- sample_columns
    results$n_genes <- nrow(data)
    results$n_samples <- length(sample_columns)
    
    # Check for missing gene IDs
    missing_genes <- is.na(data[[gene_column]]) | data[[gene_column]] == ""
    if (any(missing_genes)) {
        results$errors <- c(results$errors, 
            paste("Found", sum(missing_genes), "missing gene IDs"))
    }
    
    # Check for duplicate gene IDs
    duplicate_genes <- duplicated(data[[gene_column]])
    if (any(duplicate_genes)) {
        results$warnings <- c(results$warnings, 
            paste("Found", sum(duplicate_genes), "duplicate gene IDs"))
    }
    
    # Validate gene ID format based on type
    gene_ids <- data[[gene_column]][!missing_genes]
    if (length(gene_ids) > 0) {
        format_issues <- validate_gene_id_format(gene_ids, gene_id_type)
        if (length(format_issues) > 0) {
            results$warnings <- c(results$warnings, format_issues)
        }
    }
    
    # Check expression values are numeric
    for (col in sample_columns) {
        if (!is.numeric(data[[col]])) {
            tryCatch({
                data[[col]] <- as.numeric(data[[col]])
            }, warning = function(w) {
                results$errors <- c(results$errors, 
                    paste("Column", col, "contains non-numeric values"))
            })
        }
        
        # Check for negative values (warn but don't fail)
        if (any(data[[col]] < 0, na.rm = TRUE)) {
            results$warnings <- c(results$warnings, 
                paste("Column", col, "contains negative expression values"))
        }
    }
    
    # Overall validation
    results$is_valid <- length(results$errors) == 0
    
    return(results)
}

#' Validate gene ID format
#' @param gene_ids character vector of gene IDs
#' @param gene_id_type character, expected ID type
#' @return character vector of issues found
validate_gene_id_format <- function(gene_ids, gene_id_type) {
    issues <- character(0)
    
    if (gene_id_type == "entrez") {
        # Entrez IDs should be numeric
        non_numeric <- !grepl("^\\d+$", gene_ids)
        if (any(non_numeric)) {
            issues <- c(issues, paste("Found", sum(non_numeric), 
                "non-numeric Entrez IDs"))
        }
    } else if (gene_id_type == "ensembl") {
        # Ensembl IDs should start with ENS
        non_ensembl <- !grepl("^ENS", gene_ids)
        if (any(non_ensembl)) {
            issues <- c(issues, paste("Found", sum(non_ensembl), 
                "IDs not starting with 'ENS'"))
        }
    } else if (gene_id_type == "symbol") {
        # HGNC symbols should be uppercase letters, numbers, and some symbols
        non_symbol <- !grepl("^[A-Z0-9][A-Z0-9@-]*$", toupper(gene_ids))
        if (any(non_symbol)) {
            issues <- c(issues, paste("Found", sum(non_symbol), 
                "IDs with unusual gene symbol format"))
        }
    }
    # UniProt validation could be added here
    
    return(issues)
}

#' Create myTAI BulkPhyloExpressionSet using myTAI's proper mapping
#' @param expression_data data.frame with gene expression data (first column = gene IDs)
#' @param gene_id_type character, type of gene IDs
#' @param groups character/factor vector indicating sample groups (optional)
#' @param phylomap data.frame, phylostratum mapping (optional, will load default if NULL)
#' @param strata_legend data.frame, phylostratum legend (optional, will load default if NULL)
#' @param name character, dataset name (optional)
#' @return myTAI BulkPhyloExpressionSet object
create_bulk_phyloexpression_set <- function(expression_data, 
                                            gene_id_type = "symbol", 
                                            groups = NULL,
                                            phylomap = NULL,
                                            strata_legend = NULL,
                                            name = "Expression Dataset") {
    
    # Check if myTAI is available
    if (!requireNamespace("myTAI", quietly = TRUE)) {
        stop("myTAI package is required. Please install it first.")
    }
    
    # Validate expression data first
    validation <- validate_expression_data(expression_data, gene_id_type)
    if (!validation$is_valid) {
        stop("Expression data validation failed: ", paste(validation$errors, collapse = "; "))
    }
    
    # Load and prepare phylomap if not provided
    if (is.null(phylomap)) {
        phylomap <- load_phylomap()
    }
    
    # Load legend if not provided
    if (is.null(strata_legend)) {
        legend_file <- file.path("data", "strata_legend.tsv")
        if (file.exists(legend_file)) {
            strata_legend <- read.table(legend_file, header = TRUE, sep = "\t", 
                               stringsAsFactors = FALSE, quote = '"')
        }
    }
    
    # Convert phylomap to myTAI format (Stratum, GeneID)
    # Map user gene IDs to phylostrata using existing mapping function
    gene_column <- names(expression_data)[1]
    gene_ids <- expression_data[[gene_column]]
    
    # Use existing mapping function to convert gene IDs to appropriate type for phylomap matching
    gene_strata <- map_genes_to_phylostrata(gene_ids, id_type = gene_id_type, phylomap = phylomap)
    
    # Create the phylomap in myTAI format with only successfully mapped genes
    mapped_genes <- !is.na(gene_strata)
    mytai_phylomap <- data.frame(
        Stratum = gene_strata[mapped_genes],
        GeneID = gene_ids[mapped_genes],
        stringsAsFactors = FALSE
    )
    
    # Filter expression data to only mapped genes
    filtered_expression_data <- expression_data[mapped_genes, ]
    
    if (nrow(filtered_expression_data) == 0) {
        stop("No genes could be mapped to phylostrata")
    }
    
    # Set up groups if not provided
    sample_columns <- names(filtered_expression_data)[-1]
    if (is.null(groups)) {
        groups <- sample_columns  # Use sample names as groups
    }
    
    # Ensure groups is the right length
    if (length(groups) != length(sample_columns)) {
        stop("Length of groups (", length(groups), ") must match number of samples (", length(sample_columns), ")")
    }
    
    # Use myTAI's match_map function to create BulkPhyloExpressionSet
    tryCatch({
        # Add strata_legend to match_map call if available
        if (!is.null(strata_legend)) {
            bulk_set <- myTAI::match_map(
                data = filtered_expression_data,
                phylomap = mytai_phylomap,
                groups = groups,
                strata_legend = strata_legend,
                name = name
            )
        } else {
            bulk_set <- myTAI::match_map(
                data = filtered_expression_data,
                phylomap = mytai_phylomap,
                groups = groups,
                name = name
            )
        }
        
        # Calculate mapping statistics
        mapping_stats <- list(
            n_total = nrow(expression_data),
            n_mapped = nrow(filtered_expression_data),
            mapping_rate = round(100 * nrow(filtered_expression_data) / nrow(expression_data), 1),
            strata_distribution = table(mytai_phylomap$Stratum)
        )
        
        cat("Successfully created BulkPhyloExpressionSet:\n")
        cat("  Total genes input:", mapping_stats$n_total, "\n")
        cat("  Genes mapped to phylostrata:", mapping_stats$n_mapped, "\n")
        cat("  Mapping rate:", mapping_stats$mapping_rate, "%\n")
        
        return(list(
            phyloset = bulk_set,
            mapping_stats = mapping_stats,
            is_myTAI_object = TRUE
        ))
        
    }, error = function(e) {
        stop("Failed to create BulkPhyloExpressionSet: ", e$message)
    })
}

#' Load example expression data
#' @return data.frame with example expression data
load_example_expression_data <- function() {
    example_file <- file.path("data", "example_expression.tsv")
    if (!file.exists(example_file)) {
        stop("Example expression data file not found: ", example_file)
    }
    
    # Read the data
    data <- read.table(example_file, header = TRUE, sep = "\t", 
                      stringsAsFactors = FALSE, quote = '"')
    
    cat("Loaded example expression data:\n")
    cat("  Genes:", nrow(data), "\n")
    cat("  Samples:", ncol(data) - 1, "\n")
    cat("  Sample names:", paste(names(data)[-1], collapse = ", "), "\n")
    
    return(data)
}

#' Load phylostratum legend
#' @return data.frame with phylostratum names and descriptions
load_strata_legend <- function() {
    legend_file <- file.path("data", "strata_legend.tsv")
    if (!file.exists(legend_file)) {
        warning("Strata legend file not found: ", legend_file)
        return(NULL)
    }
    
    # Read the legend data
    legend <- read.table(legend_file, header = TRUE, sep = "\t", 
                        stringsAsFactors = FALSE, quote = '"')
    
    cat("Loaded strata legend with", nrow(legend), "phylostrata\n")
    return(legend)
}

#' Wrapper functions for myTAI plotting functions with better error handling
#' These functions provide consistent interfaces and error handling for the Shiny app

#' Create TAI signature plot using myTAI
#' @param bulk_phyloset myTAI BulkPhyloExpressionSet object
#' @param title character, plot title
#' @param show_stats logical, whether to show statistical test results
#' @return ggplot object
create_tai_signature_plot <- function(bulk_phyloset, title = "TAI Signature", show_stats = TRUE) {
    if (!requireNamespace("myTAI", quietly = TRUE)) {
        stop("myTAI package is required")
    }
    
    tryCatch({
        p <- myTAI::plot_signature(bulk_phyloset, show_p_val = show_stats)
        p <- p + ggplot2::labs(title = title)
        return(p)
    }, error = function(e) {
        stop("Error creating TAI signature plot: ", e$message)
    })
}

#' Create distribution of phylostrata plot using myTAI
#' @param bulk_phyloset myTAI BulkPhyloExpressionSet object
#' @param selected_genes character vector of genes to include (optional)
#' @param title character, plot title
#' @return ggplot object
create_distribution_strata_plot <- function(bulk_phyloset, selected_genes = NULL, title = "Distribution of Phylostrata") {
    if (!requireNamespace("myTAI", quietly = TRUE)) {
        stop("myTAI package is required")
    }
    
    tryCatch({
        if (!is.null(selected_genes) && length(selected_genes) > 0) {
            # Filter selected genes to only those present in the phyloset
            # Get the gene IDs from the phyloset using the correct slot
            phyloset_genes <- bulk_phyloset@gene_ids
            valid_selected_genes <- intersect(selected_genes, phyloset_genes)
            
            if (length(valid_selected_genes) > 0) {
                cat("Using", length(valid_selected_genes), "out of", length(selected_genes), "selected genes for distribution plot\n")
                # Use selected genes in myTAI plot_distribution_strata
                p <- myTAI::plot_distribution_strata(bulk_phyloset@strata, selected_gene_ids = valid_selected_genes)
            } else {
                cat("None of the selected genes are present in the phyloset. Using all genes for distribution plot.\n")
                p <- myTAI::plot_distribution_strata(bulk_phyloset@strata)
            }
        } else {
            p <- myTAI::plot_distribution_strata(bulk_phyloset@strata)
        }
        p <- p + ggplot2::labs(title = title)
        return(p)
    }, error = function(e) {
        stop("Error creating distribution strata plot: ", e$message)
    })
}

#' Create gene heatmap plot using myTAI
#' @param bulk_phyloset myTAI BulkPhyloExpressionSet object
#' @param selected_genes character vector of genes to plot (optional)
#' @param title character, plot title
#' @return ggplot object
create_gene_heatmap_plot <- function(bulk_phyloset, selected_genes = NULL, title = "Gene Expression Heatmap") {
    if (!requireNamespace("myTAI", quietly = TRUE)) {
        stop("myTAI package is required")
    }
    
    tryCatch({
        if (!is.null(selected_genes) && length(selected_genes) > 0) {
            # Filter selected genes to only those present in the phyloset
            # Get the gene IDs from the phyloset using the correct slot
            phyloset_genes <- bulk_phyloset@gene_ids
            valid_selected_genes <- intersect(selected_genes, phyloset_genes)
            
            if (length(valid_selected_genes) > 0) {
                cat("Using", length(valid_selected_genes), "out of", length(selected_genes), "selected genes for heatmap\n")
                # Use selected genes parameter in myTAI plot_gene_heatmap
                p <- myTAI::plot_gene_heatmap(bulk_phyloset, genes = valid_selected_genes)
            } else {
                cat("None of the selected genes are present in the phyloset. Using all genes for heatmap.\n")
                p <- myTAI::plot_gene_heatmap(bulk_phyloset)
            }
        } else {
            p <- myTAI::plot_gene_heatmap(bulk_phyloset)
        }
        p <- p + ggplot2::labs(title = title)
        return(p)
    }, error = function(e) {
        stop("Error creating gene heatmap plot: ", e$message)
    })
}

#' Create contribution plot using myTAI
#' @param bulk_phyloset myTAI BulkPhyloExpressionSet object
#' @param title character, plot title
#' @return ggplot object
create_contribution_plot <- function(bulk_phyloset, title = "Phylostratum Contribution to TAI") {
    if (!requireNamespace("myTAI", quietly = TRUE)) {
        stop("myTAI package is required")
    }
    
    tryCatch({
        p <- myTAI::plot_contribution(bulk_phyloset)
        p <- p + ggplot2::labs(title = title)
        return(p)
    }, error = function(e) {
        stop("Error creating contribution plot: ", e$message)
    })
}

#' Create gene space plot using myTAI
#' @param bulk_phyloset myTAI BulkPhyloExpressionSet object
#' @param title character, plot title
#' @return ggplot object
create_gene_space_plot <- function(bulk_phyloset, title = "Gene Space Analysis") {
    if (!requireNamespace("myTAI", quietly = TRUE)) {
        stop("myTAI package is required")
    }
    
    tryCatch({
        p <- myTAI::plot_gene_space(bulk_phyloset)
        p <- p + ggplot2::labs(title = title)
        return(p)
    }, error = function(e) {
        stop("Error creating gene space plot: ", e$message)
    })
}

#' Create sample space plot using myTAI  
#' @param bulk_phyloset myTAI BulkPhyloExpressionSet object
#' @param title character, plot title
#' @return ggplot object
create_sample_space_plot <- function(bulk_phyloset, title = "Sample Space Analysis") {
    if (!requireNamespace("myTAI", quietly = TRUE)) {
        stop("myTAI package is required")
    }
    
    tryCatch({
        p <- myTAI::plot_sample_space(bulk_phyloset)
        p <- p + ggplot2::labs(title = title)
        return(p)
    }, error = function(e) {
        stop("Error creating sample space plot: ", e$message)
    })
}

#' Create gene profiles plot using myTAI
#' @param bulk_phyloset myTAI BulkPhyloExpressionSet object
#' @param selected_genes character vector of genes to plot (optional)
#' @param interaction_colors named vector of colors for genes based on interactions (optional)
#' @param title character, plot title
#' @return ggplot object
create_gene_profiles_plot <- function(bulk_phyloset, selected_genes = NULL, interaction_colors = NULL, title = "Gene Expression Profiles") {
    if (!requireNamespace("myTAI", quietly = TRUE)) {
        stop("myTAI package is required")
    }
    
    tryCatch({
        if (!is.null(selected_genes) && length(selected_genes) > 0) {
            # Filter selected genes to only those present in the phyloset
            # Get the gene IDs from the phyloset using the correct slot
            phyloset_genes <- bulk_phyloset@gene_ids
            valid_selected_genes <- intersect(selected_genes, phyloset_genes)
            
            if (length(valid_selected_genes) > 0) {
                cat("Using", length(valid_selected_genes), "out of", length(selected_genes), "selected genes that are present in phyloset\n")
                
                if (!is.null(interaction_colors) && length(interaction_colors) > 0) {
                    # Filter colors to match valid genes
                    valid_colors <- interaction_colors[names(interaction_colors) %in% valid_selected_genes]
                    if (length(valid_colors) > 0) {
                        # Use manual coloring with interaction colors
                        p <- myTAI::plot_gene_profiles(bulk_phyloset, genes = valid_selected_genes, 
                                                     colour_by = "manual", colours = valid_colors)
                    } else {
                        # Use default gene-based coloring
                        p <- myTAI::plot_gene_profiles(bulk_phyloset, genes = valid_selected_genes)
                    }
                } else {
                    # Use default gene-based coloring
                    p <- myTAI::plot_gene_profiles(bulk_phyloset, genes = valid_selected_genes)
                }
            } else {
                cat("None of the selected genes are present in the phyloset. Using all genes.\n")
                # Use all genes with default settings
                p <- myTAI::plot_gene_profiles(bulk_phyloset)
            }
        } else {
            # Use all genes with default settings
            p <- myTAI::plot_gene_profiles(bulk_phyloset)
        }
        p <- p + ggplot2::labs(title = title)
        return(p)
    }, error = function(e) {
        stop("Error creating gene profiles plot: ", e$message)
    })
}

#' Get genes interacting with a selected gene from pathway data
#' @param selected_gene character, gene identifier
#' @param pathway_edges data.frame with pathway edge information (from server.R values$edges)
#' @param pathway_nodes data.frame with pathway node information (from server.R values$nodes)
#' @param gene_id_type character, type of gene ID for matching
#' @return list with interacting genes, interaction types, and colors
get_interacting_genes <- function(selected_gene, pathway_edges = NULL, pathway_nodes = NULL, gene_id_type = "symbol") {
    
    if (is.null(pathway_edges) || nrow(pathway_edges) == 0 || 
        is.null(pathway_nodes) || nrow(pathway_nodes) == 0) {
        return(list(
            interacting_genes = character(0),
            interaction_types = character(0),
            interaction_colors = NULL,
            message = "No pathway interaction data available"
        ))
    }
    
    # Find the node corresponding to the selected gene
    selected_node_id <- NULL
    
    # Try different matching strategies based on gene_id_type
    if (gene_id_type == "symbol") {
        # Match by HGNC symbol
        matches <- which(toupper(pathway_nodes$hgnc_symbol) == toupper(selected_gene) | 
                        toupper(pathway_nodes$label) == toupper(selected_gene))
    } else if (gene_id_type == "entrez") {
        # Match by Entrez/KEGG ID
        matches <- which(pathway_nodes$kegg_id == selected_gene)
    } else {
        # Generic matching
        matches <- which(toupper(pathway_nodes$label) == toupper(selected_gene) |
                        toupper(pathway_nodes$hgnc_symbol) == toupper(selected_gene))
    }
    
    if (length(matches) == 0) {
        return(list(
            interacting_genes = character(0),
            interaction_types = character(0),
            interaction_colors = NULL,
            message = paste("Gene", selected_gene, "not found in current pathway")
        ))
    }
    
    selected_node_id <- pathway_nodes$id[matches[1]]
    
    # Find all edges involving this node
    incoming_edges <- pathway_edges[pathway_edges$to == selected_node_id, ]
    outgoing_edges <- pathway_edges[pathway_edges$from == selected_node_id, ]
    all_edges <- rbind(incoming_edges, outgoing_edges)
    
    if (nrow(all_edges) == 0) {
        return(list(
            interacting_genes = character(0),
            interaction_types = character(0),
            interaction_colors = NULL,
            message = paste("No interactions found for gene", selected_gene)
        ))
    }
    
    # Get all interacting node IDs
    interacting_node_ids <- unique(c(
        all_edges$from[all_edges$to == selected_node_id],
        all_edges$to[all_edges$from == selected_node_id]
    ))
    
    # Remove the selected node itself
    interacting_node_ids <- interacting_node_ids[interacting_node_ids != selected_node_id]
    
    # Get gene symbols for interacting nodes
    interacting_nodes <- pathway_nodes[pathway_nodes$id %in% interacting_node_ids, ]
    
    # Prioritize HGNC symbols, fallback to labels
    interacting_genes <- ifelse(
        !is.na(interacting_nodes$hgnc_symbol) & interacting_nodes$hgnc_symbol != "",
        interacting_nodes$hgnc_symbol,
        interacting_nodes$label
    )
    
    # Get interaction types and create color mapping
    interaction_types <- character(0)
    interaction_colors <- character(0)
    
    # Create color mapping based on edge relationships
    edge_colors <- c(
        "activation" = "#4CAF50",       # Green
        "inhibition" = "#F44336",       # Red  
        "binding" = "#2196F3",          # Blue
        "expression" = "#FF9800",       # Orange
        "catalysis" = "#9C27B0",        # Purple
        "reaction" = "#607D8B",         # Blue-grey
        "unknown" = "#9E9E9E"           # Grey
    )
    
    for (i in seq_along(interacting_genes)) {
        node_id <- interacting_nodes$id[i]
        gene_name <- interacting_genes[i]
        
        # Find edges involving this gene
        gene_edges <- all_edges[all_edges$from == node_id | all_edges$to == node_id, ]
        
        if (nrow(gene_edges) > 0) {
            # Use the most common relationship type or first one
            relationship <- gene_edges$relationship[1]
            interaction_types <- c(interaction_types, relationship)
            
            # Assign color based on relationship
            color <- edge_colors[relationship]
            if (is.na(color)) {
                color <- edge_colors["unknown"]
            }
            interaction_colors <- c(interaction_colors, color)
        } else {
            interaction_types <- c(interaction_types, "unknown")
            interaction_colors <- c(interaction_colors, edge_colors["unknown"])
        }
    }
    
    # Create named color vector for myTAI plotting
    names(interaction_colors) <- interacting_genes
    
    # Include the selected gene itself with a special color
    all_genes <- c(selected_gene, interacting_genes)
    all_colors <- c("#000000", interaction_colors)  # Black for selected gene
    names(all_colors) <- all_genes
    
    return(list(
        interacting_genes = interacting_genes,
        interaction_types = interaction_types,
        interaction_colors = all_colors,
        selected_gene = selected_gene,
        message = paste("Found", length(interacting_genes), "interacting genes for", selected_gene)
    ))
}

#' Get genes from pathway nodes
#' @param pathway_nodes data.frame with pathway node information
#' @param gene_id_type character, preferred gene ID type to return
#' @return character vector of gene IDs from the pathway
get_pathway_genes <- function(pathway_nodes = NULL, gene_id_type = "symbol") {
    
    if (is.null(pathway_nodes) || nrow(pathway_nodes) == 0) {
        return(character(0))
    }
    
    # Extract genes based on preferred ID type
    if (gene_id_type == "symbol") {
        # Prioritize HGNC symbols, fallback to labels
        genes <- ifelse(
            !is.na(pathway_nodes$hgnc_symbol) & pathway_nodes$hgnc_symbol != "",
            pathway_nodes$hgnc_symbol,
            pathway_nodes$label
        )
    } else if (gene_id_type == "entrez") {
        # Use KEGG IDs (which are often Entrez-based)
        genes <- pathway_nodes$kegg_id
    } else {
        # Default to labels
        genes <- pathway_nodes$label
    }
    
    # Remove empty or NA values
    genes <- genes[!is.na(genes) & genes != ""]
    
    return(unique(genes))
}

#' Generate summary statistics for phylostratum mapping
#' @param mapping_stats list with mapping statistics
#' @return character vector with formatted summary
format_mapping_summary <- function(mapping_stats) {
    if (is.null(mapping_stats)) {
        return("No mapping statistics available")
    }
    
    summary_lines <- c(
        paste("Total genes:", mapping_stats$n_total),
        paste("Successfully mapped:", mapping_stats$n_mapped),
        paste("Mapping rate:", mapping_stats$mapping_rate, "%")
    )
    
    # Don't include phylostrata distribution in summary
    return(summary_lines)
}
