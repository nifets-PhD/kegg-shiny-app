# Test the corrected phylostratum logic for KEGG pathways
# This should now use Entrez IDs instead of HGNC symbols

# Source the functions
source("R/kegg_utils.R")

# Load a cached pathway
cache_files <- list.files("kegg_cache", pattern = "\\.rds$", full.names = TRUE)

if (length(cache_files) > 0) {
    # Load the first cached pathway
    test_file <- cache_files[1]
    cat("Loading cached pathway from:", test_file, "\n")
    
    pathway_data <- readRDS(test_file)
    nodes <- pathway_data$nodes
    
    cat("\n=== TESTING CORRECTED PHYLOSTRATUM LOGIC ===\n")
    
    # Show what IDs we have
    if (!is.null(nodes$kegg_id)) {
        kegg_id_available <- sum(!is.na(nodes$kegg_id) & nodes$kegg_id != "")
        cat("Entrez IDs (kegg_id) available:", kegg_id_available, "out of", nrow(nodes), "\n")
        if (kegg_id_available > 0) {
            test_entrez_ids <- nodes$kegg_id[!is.na(nodes$kegg_id) & nodes$kegg_id != ""][1:min(5, kegg_id_available)]
            cat("Sample Entrez IDs:", paste(test_entrez_ids, collapse = ", "), "\n")
        }
    }
    
    # Test the apply_phylostratum_coloring function directly
    cat("\n=== TESTING apply_phylostratum_coloring() ===\n")
    
    # This should now prioritize Entrez IDs
    test_nodes <- nodes[1:10, ]  # Test with first 10 nodes
    result_nodes <- apply_phylostratum_coloring(test_nodes, id_type = "uniprot")  # Even with uniprot, should use Entrez
    
    cat("\nResult: Function should have used Entrez IDs regardless of id_type parameter\n")
    
} else {
    cat("No cached pathway files found. Please load a pathway first.\n")
}
