#!/usr/bin/env Rscript

# Test script for evolutionary transcriptomics functionality
source('R/kegg_utils.R')
source('R/evolutionary_utils.R')

# Load example data and create phyloset
cat("Loading example data...\n")
example_data <- load_example_expression_data()

cat("\nLoading strata legend...\n")
legend <- load_strata_legend()

cat("\nCreating BulkPhyloExpressionSet with strata_legend...\n")
result <- create_bulk_phyloexpression_set(
    expression_data = example_data,
    gene_id_type = 'symbol',
    strata_legend = legend,
    name = 'Example Dataset'
)

# Test plotting functions
cat("\nTesting TAI signature plot...\n")
tryCatch({
    p1 <- create_tai_signature_plot(result$phyloset, title = 'Test TAI Signature')
    cat("✓ TAI signature plot created successfully!\n")
    cat("  Plot class:", class(p1), "\n")
}, error = function(e) {
    cat("✗ Error creating TAI signature plot:", e$message, "\n")
})

cat("\nTesting distribution strata plot...\n")
tryCatch({
    p2 <- create_distribution_strata_plot(result$phyloset)
    cat("✓ Distribution strata plot created successfully!\n")
}, error = function(e) {
    cat("✗ Error creating distribution strata plot:", e$message, "\n")
})

cat("\nTesting gene heatmap plot...\n")
tryCatch({
    p3 <- create_gene_heatmap_plot(result$phyloset)
    cat("✓ Gene heatmap plot created successfully!\n")
}, error = function(e) {
    cat("✗ Error creating gene heatmap plot:", e$message, "\n")
})

cat("\nTesting sample space plot...\n")
tryCatch({
    p4 <- create_sample_space_plot(result$phyloset)
    cat("✓ Sample space plot created successfully!\n")
}, error = function(e) {
    cat("✗ Error creating sample space plot:", e$message, "\n")
})

# Print summary
cat("\n=== SUMMARY ===\n")
cat("BulkPhyloExpressionSet object class:", class(result$phyloset), "\n")
cat("Mapping statistics:\n")
stats <- result$mapping_stats
cat("  Total genes:", stats$n_total, "\n")
cat("  Mapped genes:", stats$n_mapped, "\n") 
cat("  Mapping rate:", stats$mapping_rate, "%\n")
cat("  Phylostrata count:", length(stats$strata_distribution), "\n")

cat("\nAll tests completed!\n")
