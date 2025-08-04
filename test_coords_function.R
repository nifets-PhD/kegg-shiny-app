#!/usr/bin/env Rscript

# Test get_kegg_node_coordinates function directly
cat("Testing get_kegg_node_coordinates function...\n")

# Source the function
source('R/kegg_utils.R')

# Test the function
tryCatch({
    cat("Calling get_kegg_node_coordinates('hsa04152')...\n")
    coords <- get_kegg_node_coordinates('hsa04152')
    
    if (is.null(coords)) {
        cat("Function returned NULL\n")
    } else {
        cat("Function returned data with", nrow(coords), "entries\n")
        
        # Look for MTOR
        mtor_entries <- coords[coords$entry_id == '92', ]
        if (nrow(mtor_entries) > 0) {
            cat("MTOR entry found:\n")
            print(mtor_entries[c('entry_id', 'kegg_name', 'label')])
        } else {
            cat("MTOR entry not found\n")
            cat("Available entry IDs:", head(coords$entry_id), "\n")
        }
    }
}, error = function(e) {
    cat("Error occurred:", e$message, "\n")
    traceback()
})
