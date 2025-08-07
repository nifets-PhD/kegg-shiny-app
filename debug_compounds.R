#!/usr/bin/env Rscript

# Debug script for compound labeling issue
cat("===     
    # Check for any nodes that might be compounds based on kegg_name  
    kegg_compounds <- nodes[grepl("^cpd:", nodes$kegg_name), ]
    cat('\nNodes with compound-like kegg_name (cpd:):', nrow(kegg_compounds), "\n")
    
    if (nrow(kegg_compounds) > 0) {
        cat("Sample compound-like node:\n")
        print(kegg_compounds[1, c('entry_id', 'label', 'entry_type', 'kegg_name', 'description')])
    } COMPOUND LABELING ===\n")

# Load required functions
source('R/kegg_utils.R')
source('R/gene_consolidation_utils.R')

# Test with mTOR pathway which has compounds
pathway_id <- 'hsa04152'
cat("Testing pathway:", pathway_id, "\n\n")

# Step 1: Check raw XML parsing
cat("1. CHECKING RAW XML PARSING:\n")
xml_data <- xml2::read_xml('kegg_cache/hsa04152.xml')
nodes_ns <- xml2::xml_find_all(xml_data, './/entry')

# Find first few compound entries
compound_entries <- c()
for(i in 1:min(20, length(nodes_ns))) {
    entry <- nodes_ns[[i]]
    entry_type <- xml2::xml_attr(entry, 'type')
    if (entry_type == 'compound') {
        entry_id <- xml2::xml_attr(entry, 'id')
        entry_name <- xml2::xml_attr(entry, 'name')
        
        # Get graphics info
        graphics <- xml2::xml_find_first(entry, './/graphics')
        label_text <- if (!is.null(graphics)) xml2::xml_attr(graphics, 'name') else NA
        
        compound_entries <- rbind(compound_entries, data.frame(
            id = entry_id,
            name = entry_name,
            label = label_text,
            stringsAsFactors = FALSE
        ))
        
        cat("Entry", entry_id, "- Type:", entry_type, "Name:", entry_name, "Label:", label_text, "\n")
        
        if (nrow(compound_entries) >= 3) break
    }
}

cat("\nFound", nrow(compound_entries), "compound entries in XML\n\n")

# Step 2: Test full coordinate extraction
cat("2. TESTING FULL COORDINATE EXTRACTION:\n")
nodes <- get_kegg_node_coordinates(pathway_id)
cat("Total nodes returned:", nrow(nodes), "\n")

if (nrow(nodes) > 0) {
    cat("Columns in nodes data:\n")
    print(colnames(nodes))
    
    cat("\nNode types:\n")
    if ("type" %in% colnames(nodes)) {
        print(table(nodes$type, useNA = "always"))
    } else {
        cat("No 'type' column found!\n")
    }
    
    # Check for compound nodes using entry_type instead of type
    if ("entry_type" %in% colnames(nodes)) {
        print(table(nodes$entry_type, useNA = "always"))
        compound_nodes <- nodes[!is.na(nodes$entry_type) & nodes$entry_type == 'compound', ]
        cat("\nCompound nodes found:", nrow(compound_nodes), "\n")
        
        if (nrow(compound_nodes) > 0) {
            cat("Sample compound node info:\n")
            sample_compound <- compound_nodes[1, ]
            print(sample_compound[c('entry_id', 'label', 'entry_type', 'kegg_name', 'description')])
        }
    } else {
        cat("No 'entry_type' column found either!\n")
    }
    
    # Check for any nodes that might be compounds based on gene_name (use kegg_name instead)
    kegg_compounds <- nodes[grepl("^cpd:", nodes$kegg_name), ]
    cat("\nNodes with compound-like kegg_name (cpd:):", nrow(kegg_compounds), "\n")
    
    if (nrow(kegg_compounds) > 0) {
        cat("Sample compound-like node:\n")
        print(kegg_compounds[1, c('id', 'label', 'type', 'gene_name', 'description')])
    }
    
} else {
    cat("ERROR: No nodes returned from get_kegg_node_coordinates!\n")
}

cat("\n=== DEBUG COMPLETE ===\n")
