# Test phylostratum coloring filtering for non-gene nodes
# This should verify that pathway titles are not colored

# Source the functions
source("R/kegg_utils.R")

# Load a cached pathway
cache_files <- list.files("kegg_cache", pattern = "\\.rds$", full.names = TRUE)

if (length(cache_files) > 0) {
    test_file <- cache_files[1]
    cat("Loading cached pathway from:", test_file, "\n")
    
    pathway_data <- readRDS(test_file)
    nodes <- pathway_data$nodes
    
    cat("\n=== NODE TYPE ANALYSIS ===\n")
    
    # Show node types
    if (!is.null(nodes$type)) {
        type_summary <- table(nodes$type, useNA = "always")
        cat("Node types:\n")
        print(type_summary)
    }
    
    # Look for pathway titles
    title_nodes <- nodes[grepl("^TITLE:", nodes$label, ignore.case = TRUE), ]
    cat("\nPathway title nodes found:", nrow(title_nodes), "\n")
    if (nrow(title_nodes) > 0) {
        cat("Title labels:\n")
        for (i in 1:nrow(title_nodes)) {
            cat(" -", title_nodes$label[i], "\n")
            cat("   Type:", title_nodes$type[i], "\n") 
            cat("   KEGG_ID:", title_nodes$kegg_id[i], "\n")
            cat("   HGNC:", title_nodes$hgnc_symbol[i], "\n")
        }
    }
    
    cat("\n=== TESTING PHYLOSTRATUM FILTERING ===\n")
    
    # Test the corrected apply_phylostratum_coloring function
    cat("Before phylostratum coloring:\n")
    cat("Nodes with color.background:", sum(!is.na(nodes$color.background) & nodes$color.background != ""), "\n")
    
    result_nodes <- apply_phylostratum_coloring(nodes, id_type = "entrez")
    
    cat("\nAfter phylostratum coloring:\n")
    cat("Nodes with color.background:", sum(!is.na(result_nodes$color.background) & result_nodes$color.background != ""), "\n")
    
    # Check if title nodes were colored
    if (nrow(title_nodes) > 0) {
        title_indices <- which(grepl("^TITLE:", nodes$label, ignore.case = TRUE))
        colored_titles <- sum(!is.na(result_nodes$color.background[title_indices]) & 
                            result_nodes$color.background[title_indices] != "")
        cat("Title nodes colored:", colored_titles, "out of", length(title_indices), "\n")
        
        if (colored_titles > 0) {
            cat("❌ ERROR: Title nodes should NOT be colored!\n")
        } else {
            cat("✅ SUCCESS: Title nodes were not colored\n")
        }
    }
    
    # Check gene nodes were colored
    gene_indices <- which(nodes$type == "gene" | is.na(nodes$type))
    colored_genes <- sum(!is.na(result_nodes$color.background[gene_indices]) & 
                        result_nodes$color.background[gene_indices] != "")
    cat("Gene nodes colored:", colored_genes, "out of", length(gene_indices), "\n")
    
} else {
    cat("No cached pathway files found. Please load a pathway first.\n")
}
