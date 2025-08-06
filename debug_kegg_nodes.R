# Debug script to check KEGG pathway node structure
# This will help us understand what IDs are actually available

# Source the functions
source("R/kegg_utils.R")

# Test with a simple pathway (let's use a cached one)
cache_files <- list.files("kegg_cache", pattern = "\\.rds$", full.names = TRUE)

if (length(cache_files) > 0) {
    # Load the first cached pathway
    test_file <- cache_files[1]
    cat("Loading cached pathway from:", test_file, "\n")
    
    pathway_data <- readRDS(test_file)
    nodes <- pathway_data$nodes
    
    cat("\n=== NODE STRUCTURE ANALYSIS ===\n")
    cat("Number of nodes:", nrow(nodes), "\n")
    cat("Column names:", paste(names(nodes), collapse = ", "), "\n\n")
    
    if (!is.null(nodes$hgnc_symbol)) {
        hgnc_available <- sum(!is.na(nodes$hgnc_symbol) & nodes$hgnc_symbol != "")
        cat("HGNC symbols available:", hgnc_available, "out of", nrow(nodes), "\n")
        if (hgnc_available > 0) {
            cat("Sample HGNC symbols:", paste(head(nodes$hgnc_symbol[!is.na(nodes$hgnc_symbol) & nodes$hgnc_symbol != ""], 5), collapse = ", "), "\n")
        }
    } else {
        cat("HGNC symbol column: NOT PRESENT\n")
    }
    
    if (!is.null(nodes$label)) {
        label_available <- sum(!is.na(nodes$label) & nodes$label != "")
        cat("Labels available:", label_available, "out of", nrow(nodes), "\n")
        if (label_available > 0) {
            cat("Sample labels:", paste(head(nodes$label[!is.na(nodes$label) & nodes$label != ""], 5), collapse = ", "), "\n")
        }
    } else {
        cat("Label column: NOT PRESENT\n")
    }
    
    if (!is.null(nodes$kegg_id)) {
        kegg_id_available <- sum(!is.na(nodes$kegg_id) & nodes$kegg_id != "")
        cat("KEGG IDs available:", kegg_id_available, "out of", nrow(nodes), "\n")
        if (kegg_id_available > 0) {
            cat("Sample KEGG IDs:", paste(head(nodes$kegg_id[!is.na(nodes$kegg_id) & nodes$kegg_id != ""], 5), collapse = ", "), "\n")
        }
    } else {
        cat("KEGG ID column: NOT PRESENT\n")
    }
    
    # Test the phylostratum coloring logic
    cat("\n=== TESTING PHYLOSTRATUM LOGIC ===\n")
    
    # Test what the function would choose
    if (!is.null(nodes$hgnc_symbol) && any(!is.na(nodes$hgnc_symbol) & nodes$hgnc_symbol != "")) {
        cat("✅ Function SHOULD use HGNC symbols\n")
        test_genes <- nodes$hgnc_symbol[!is.na(nodes$hgnc_symbol) & nodes$hgnc_symbol != ""][1:min(5, sum(!is.na(nodes$hgnc_symbol) & nodes$hgnc_symbol != ""))]
        cat("Testing with genes:", paste(test_genes, collapse = ", "), "\n")
        
        # Test the mapping
        phylo_result <- map_genes_to_phylostrata(test_genes, id_type = "symbol")
        cat("Mapping result:", sum(!is.na(phylo_result)), "out of", length(test_genes), "found\n")
        
    } else {
        cat("❌ Function will NOT find HGNC symbols\n")
        if (!is.null(nodes$label)) {
            cat("Will fall back to labels\n")
        }
    }
    
} else {
    cat("No cached pathway files found. Please load a pathway first.\n")
}
