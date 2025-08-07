#!/usr/bin/env Rscript

# Test compound name mapping fix
cat("=== TESTING COMPOUND NAME MAPPING FIX ===\n")

# Load global settings and functions
source('global.R')
source('R/gene_consolidation_utils.R')

# Test compound info function
cat("1. Testing compound info function:\n")
compound_test <- get_compound_info('cpd:C00031')
cat('Compound C00031 info:\n')
print(compound_test)

# Test another compound
compound_test2 <- get_compound_info('C00162')
cat('\nCompound C00162 info:\n')
print(compound_test2)

# Test consolidation with mTOR pathway
cat("\n2. Testing compound consolidation in mTOR pathway:\n")
nodes <- get_kegg_node_coordinates('hsa04152')
compound_nodes <- nodes[nodes$type == 'compound' & !is.na(nodes$type), ]
cat('Found', nrow(compound_nodes), 'compound nodes\n')

if (nrow(compound_nodes) > 0) {
    cat('Sample compound node:\n')
    sample <- compound_nodes[1, ]
    print(sample[c('id', 'label', 'type', 'gene_name', 'description')])
}

cat("\n=== TEST COMPLETE ===\n")
