#!/usr/bin/env Rscript

# Test script to check if the Shiny app can start without errors
cat("Testing Shiny app startup...\n")

# Source the global.R to check for any loading errors
tryCatch({
    source("global.R")
    cat("✓ global.R loaded successfully\n")
}, error = function(e) {
    cat("✗ Error loading global.R:", e$message, "\n")
    quit(status = 1)
})

# Check if UI can be created
tryCatch({
    source("ui.R")
    cat("✓ ui.R loaded successfully\n")
    cat("✓ UI structure created\n")
}, error = function(e) {
    cat("✗ Error loading ui.R:", e$message, "\n")
    quit(status = 1)
})

# Check if server can be created (without running it)
tryCatch({
    source("server.R")
    cat("✓ server.R loaded successfully\n")
}, error = function(e) {
    cat("✗ Error loading server.R:", e$message, "\n")
    quit(status = 1)
})

# Test evolutionary functions are accessible
tryCatch({
    # Test function exists
    if (exists("create_bulk_phyloexpression_set")) {
        cat("✓ Evolutionary transcriptomics functions are available\n")
    } else {
        cat("✗ Evolutionary transcriptomics functions not found\n")
    }
    
    # Test example data loading
    example_data <- load_example_expression_data()
    cat("✓ Example expression data loading works\n")
    
    # Test strata legend loading
    legend <- load_strata_legend()
    if (!is.null(legend)) {
        cat("✓ Strata legend loading works\n")
    } else {
        cat("! Strata legend not found (optional)\n")
    }
    
}, error = function(e) {
    cat("✗ Error testing evolutionary functions:", e$message, "\n")
})

cat("\n=== STARTUP TEST SUMMARY ===\n")
cat("All core components loaded successfully!\n")
cat("The Shiny app should be ready to run.\n")
cat("\nTo start the app, run:\n")
cat("  shiny::runApp()\n")
cat("Or in RStudio: click the 'Run App' button\n")
