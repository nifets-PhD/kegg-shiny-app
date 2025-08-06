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
        cat("  Phylostrata represented:", paste(names(mapping_stats$strata_distribution), collapse = ", "), "\n")
        
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

#' Create distribution of strata plot using myTAI
#' @param bulk_phyloset myTAI BulkPhyloExpressionSet object
#' @param title character, plot title
#' @return ggplot object
create_distribution_strata_plot <- function(bulk_phyloset, title = "Distribution of Phylostrata") {
    if (!requireNamespace("myTAI", quietly = TRUE)) {
        stop("myTAI package is required")
    }
    
    tryCatch({
        # Try different approaches for plotting distribution
        # First try with the phyloset directly
        myTAI::plot_distribution_strata(bulk_phyloset@strata) + ggplot2::labs(title = title)
    }, error = function(e) {
        stop("Error creating distribution strata plot: ", e$message)
    })
}

#' Create gene heatmap plot using myTAI
#' @param bulk_phyloset myTAI BulkPhyloExpressionSet object
#' @param title character, plot title
#' @return ggplot object
create_gene_heatmap_plot <- function(bulk_phyloset, title = "Gene Expression Heatmap") {
    if (!requireNamespace("myTAI", quietly = TRUE)) {
        stop("myTAI package is required")
    }
    
    tryCatch({
        p <- myTAI::plot_gene_heatmap(bulk_phyloset)
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
#' @param title character, plot title
#' @return ggplot object
create_gene_profiles_plot <- function(bulk_phyloset, selected_genes = NULL, title = "Gene Expression Profiles") {
    if (!requireNamespace("myTAI", quietly = TRUE)) {
        stop("myTAI package is required")
    }
    
    tryCatch({
        if (!is.null(selected_genes)) {
            # Filter to selected genes using myTAI's select_genes method
            filtered_set <- myTAI::select_genes(bulk_phyloset, selected_genes)
            p <- myTAI::plot_gene_profiles(filtered_set)
        } else {
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
#' @param pathway_edges data.frame with pathway edge information
#' @param gene_id_type character, type of gene ID for matching
#' @return list with interacting genes and interaction types
get_interacting_genes <- function(selected_gene, pathway_edges = NULL, gene_id_type = "symbol") {
    
    if (is.null(pathway_edges) || nrow(pathway_edges) == 0) {
        return(list(
            interacting_genes = character(0),
            interaction_types = character(0),
            message = "No pathway interaction data available"
        ))
    }
    
    # Find edges involving the selected gene
    # This depends on the structure of pathway_edges
    # For now, return empty result with informative message
    
    return(list(
        interacting_genes = character(0),
        interaction_types = character(0),
        message = paste("Interaction analysis for", selected_gene, "not yet implemented")
    ))
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
        paste("Mapping rate:", mapping_stats$mapping_rate, "%"),
        "",
        "Phylostrata distribution:"
    )
    
    if (!is.null(mapping_stats$strata_distribution)) {
        strata_lines <- paste("  Stratum", names(mapping_stats$strata_distribution), ":",
                             mapping_stats$strata_distribution, "genes")
        summary_lines <- c(summary_lines, strata_lines)
    }
    
    return(summary_lines)
}
