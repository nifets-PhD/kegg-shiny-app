# KEGG Enrichment Analysis Utility Functions
# Functions for performing KEGG pathway enrichment analysis

#' Convert gene symbols to Entrez IDs
#' @param gene_symbols Character vector of HGNC gene symbols
#' @return Character vector of Entrez IDs
convert_symbols_to_entrez <- function(gene_symbols) {
    tryCatch({
        # Remove any empty or NA symbols
        gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]
        
        if (length(gene_symbols) == 0) {
            return(character(0))
        }
        
        # Convert symbols to Entrez IDs using org.Hs.eg.db
        entrez_ids <- clusterProfiler::bitr(
            gene_symbols, 
            fromType = "SYMBOL", 
            toType = "ENTREZID", 
            OrgDb = org.Hs.eg.db::org.Hs.eg.db,
            drop = TRUE
        )
        
        return(entrez_ids$ENTREZID)
        
    }, error = function(e) {
        cat("Error converting gene symbols to Entrez IDs:", e$message, "\n")
        return(character(0))
    })
}

#' Perform KEGG pathway enrichment analysis
#' @param gene_symbols Character vector of HGNC gene symbols
#' @param organism KEGG organism code (default: "hsa" for human)
#' @param pvalue_cutoff P-value cutoff for significance (default: 0.05)
#' @param qvalue_cutoff Q-value cutoff for FDR correction (default: 0.2)
#' @param progress Optional progress object for showing progress updates
#' @return enrichResult object from clusterProfiler
perform_kegg_enrichment <- function(gene_symbols, organism = "hsa", 
                                   pvalue_cutoff = 0.05, qvalue_cutoff = 0.2, progress = NULL) {
    tryCatch({
        cat("Starting KEGG enrichment analysis for", length(gene_symbols), "genes\n")
        
        # Update progress if provided
        if (!is.null(progress)) {
            progress$set(message = "Converting gene symbols to Entrez IDs...", value = 0.3)
        }
        
        # Convert gene symbols to Entrez IDs
        entrez_ids <- convert_symbols_to_entrez(gene_symbols)
        
        if (length(entrez_ids) == 0) {
            stop("No valid Entrez IDs found for the provided gene symbols")
        }
        
        cat("Converted", length(gene_symbols), "symbols to", length(entrez_ids), "Entrez IDs\n")
        
        # Update progress if provided
        if (!is.null(progress)) {
            progress$set(message = "Running KEGG pathway enrichment...", value = 0.6)
        }
        
        # Perform KEGG enrichment
        kegg_result <- clusterProfiler::enrichKEGG(
            gene = entrez_ids,
            organism = organism,
            keyType = "ncbi-geneid",
            pvalueCutoff = pvalue_cutoff,
            pAdjustMethod = "BH",
            qvalueCutoff = qvalue_cutoff,
            use_internal_data = FALSE
        )
        
        # Update progress if provided
        if (!is.null(progress)) {
            progress$set(message = "Analyzing results...", value = 0.8)
        }
        
        if (is.null(kegg_result) || nrow(kegg_result@result) == 0) {
            cat("No significant pathways found\n")
            return(NULL)
        }
        
        cat("Found", nrow(kegg_result@result), "significant pathways\n")
        return(kegg_result)
        
    }, error = function(e) {
        cat("Error in KEGG enrichment analysis:", e$message, "\n")
        return(NULL)
    })
}

#' Format enrichment results for display
#' @param enrichment_result enrichResult object from clusterProfiler
#' @param gene_symbols Original gene symbols for mapping back
#' @return Data frame with formatted results
format_enrichment_results <- function(enrichment_result, gene_symbols = NULL) {
    if (is.null(enrichment_result) || nrow(enrichment_result@result) == 0) {
        return(data.frame(
            Pathway = character(0),
            Description = character(0),
            GeneRatio = character(0),
            BgRatio = character(0),
            pvalue = numeric(0),
            p.adjust = numeric(0),
            qvalue = numeric(0),
            Count = integer(0),
            FoldEnrichment = numeric(0),
            Genes = character(0)
        ))
    }
    
    tryCatch({
        result_df <- enrichment_result@result
        
        # Calculate fold enrichment
        gene_ratios <- sapply(strsplit(result_df$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
        bg_ratios <- sapply(strsplit(result_df$BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
        fold_enrichment <- gene_ratios / bg_ratios
        
        # Convert Entrez IDs back to gene symbols if provided
        genes_formatted <- result_df$geneID
        if (!is.null(gene_symbols)) {
            # This is a simplified mapping - in practice you might want more sophisticated mapping
            genes_formatted <- sapply(strsplit(result_df$geneID, "/"), function(entrez_list) {
                # For now, just return the Entrez IDs formatted nicely
                paste(entrez_list, collapse = ", ")
            })
        }
        
        formatted_results <- data.frame(
            Pathway = result_df$ID,
            Description = result_df$Description,
            GeneRatio = result_df$GeneRatio,
            BgRatio = result_df$BgRatio,
            pvalue = result_df$pvalue,  # Keep full precision
            p.adjust = result_df$p.adjust,  # Keep full precision
            qvalue = result_df$qvalue,  # Keep full precision
            Count = result_df$Count,
            FoldEnrichment = round(fold_enrichment, 2),
            Genes = genes_formatted,
            stringsAsFactors = FALSE
        )
        
        # Sort by p.adjust
        formatted_results <- formatted_results[order(formatted_results$p.adjust), ]
        
        return(formatted_results)
        
    }, error = function(e) {
        cat("Error formatting enrichment results:", e$message, "\n")
        return(data.frame())
    })
}

#' Create enrichment plot data for visualization
#' @param enrichment_result enrichResult object from clusterProfiler
#' @param top_n Number of top pathways to include (default: 20)
#' @return Data frame suitable for plotting
create_enrichment_plot_data <- function(enrichment_result, top_n = 20) {
    if (is.null(enrichment_result) || nrow(enrichment_result@result) == 0) {
        return(NULL)
    }
    
    tryCatch({
        result_df <- enrichment_result@result
        
        # Take top N most significant pathways
        top_pathways <- head(result_df[order(result_df$p.adjust), ], top_n)
        
        # Calculate fold enrichment
        gene_ratios <- sapply(strsplit(top_pathways$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
        bg_ratios <- sapply(strsplit(top_pathways$BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
        fold_enrichment <- gene_ratios / bg_ratios
        
        plot_data <- data.frame(
            Pathway = top_pathways$Description,
            FoldEnrichment = fold_enrichment,
            pvalue = top_pathways$pvalue,
            p.adjust = top_pathways$p.adjust,
            Count = top_pathways$Count,
            # Truncate long pathway names for better display
            PathwayShort = ifelse(nchar(top_pathways$Description) > 50, 
                                paste0(substr(top_pathways$Description, 1, 47), "..."), 
                                top_pathways$Description),
            stringsAsFactors = FALSE
        )
        
        # Order by fold enrichment for plotting
        plot_data <- plot_data[order(plot_data$FoldEnrichment, decreasing = TRUE), ]
        plot_data$PathwayShort <- factor(plot_data$PathwayShort, levels = plot_data$PathwayShort)
        
        return(plot_data)
        
    }, error = function(e) {
        cat("Error creating enrichment plot data:", e$message, "\n")
        return(NULL)
    })
}

#' Validate gene symbols
#' @param gene_symbols Character vector of gene symbols
#' @return List with valid genes, invalid genes, and conversion statistics
validate_gene_symbols <- function(gene_symbols) {
    tryCatch({
        # Remove empty and NA values
        clean_genes <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]
        clean_genes <- unique(trimws(clean_genes))
        
        if (length(clean_genes) == 0) {
            return(list(
                valid_genes = character(0),
                invalid_genes = character(0),
                entrez_mapping = data.frame(),
                conversion_rate = 0
            ))
        }
        
        # Try to convert to Entrez IDs
        conversion_result <- clusterProfiler::bitr(
            clean_genes, 
            fromType = "SYMBOL", 
            toType = "ENTREZID", 
            OrgDb = org.Hs.eg.db::org.Hs.eg.db,
            drop = FALSE
        )
        
        valid_genes <- conversion_result$SYMBOL[!is.na(conversion_result$ENTREZID)]
        invalid_genes <- setdiff(clean_genes, valid_genes)
        conversion_rate <- length(valid_genes) / length(clean_genes)
        
        return(list(
            valid_genes = valid_genes,
            invalid_genes = invalid_genes,
            entrez_mapping = conversion_result,
            conversion_rate = conversion_rate
        ))
        
    }, error = function(e) {
        cat("Error validating gene symbols:", e$message, "\n")
        return(list(
            valid_genes = character(0),
            invalid_genes = gene_symbols,
            entrez_mapping = data.frame(),
            conversion_rate = 0
        ))
    })
}
