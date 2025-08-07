# Debug script to check why all nodes are showing as "gene" type
# Load the kegg_utils functions
source("R/kegg_utils.R")

cat("=== DEBUGGING NODE TYPES ===\n")

# Test the coordinate extraction function directly
cat("\n1. Testing coordinate extraction...\n")
coords <- get_kegg_node_coordinates('hsa05200')

if (!is.null(coords)) {
  cat('Node types found in coordinates:', paste(unique(coords$entry_type), collapse=', '), '\n')
  cat('Total nodes:', nrow(coords), '\n')
  cat('Gene nodes:', sum(coords$entry_type == 'gene'), '\n')
  cat('Compound nodes:', sum(coords$entry_type == 'compound'), '\n')  
  cat('Map nodes:', sum(coords$entry_type == 'map'), '\n')
  
  cat('\nFirst 10 nodes:\n')
  print(coords[1:10, c('entry_id', 'entry_type', 'label', 'kegg_name')])
  
  # Check specific node types
  cat('\nCompound nodes:\n')
  compound_nodes <- coords[coords$entry_type == 'compound', ]
  if (nrow(compound_nodes) > 0) {
    print(compound_nodes[, c('entry_id', 'entry_type', 'label', 'kegg_name')])
  } else {
    cat("No compound nodes found in coordinates\n")
  }
  
  cat('\nMap nodes:\n')
  map_nodes <- coords[coords$entry_type == 'map', ]
  if (nrow(map_nodes) > 0) {
    print(map_nodes[, c('entry_id', 'entry_type', 'label', 'kegg_name')])
  } else {
    cat("No map nodes found in coordinates\n")
  }
} else {
  cat('No coordinates returned\n')
}

# Test the full pathway parsing function
cat("\n2. Testing full pathway parsing...\n")
pathway_data <- parse_kegg_pathway_with_hsa('hsa05200', use_cached = TRUE)

if (!is.null(pathway_data) && !is.null(pathway_data$nodes)) {
  nodes <- pathway_data$nodes
  cat('Node types found in parsed data:', paste(unique(nodes$type), collapse=', '), '\n')
  cat('Total nodes:', nrow(nodes), '\n')
  cat('Gene nodes:', sum(nodes$type == 'gene'), '\n')
  cat('Compound nodes:', sum(nodes$type == 'compound'), '\n')  
  cat('Map nodes:', sum(nodes$type == 'map'), '\n')
  
  cat('\nFirst 10 nodes from parsed data:\n')
  print(nodes[1:10, c('id', 'type', 'label', 'gene_name')])
  
  # Check specific node types
  cat('\nCompound nodes from parsed data:\n')
  compound_nodes <- nodes[nodes$type == 'compound', ]
  if (nrow(compound_nodes) > 0) {
    print(compound_nodes[, c('id', 'type', 'label', 'gene_name')])
  } else {
    cat("No compound nodes found in parsed data\n")
  }
  
  cat('\nMap nodes from parsed data:\n')
  map_nodes <- nodes[nodes$type == 'map', ]
  if (nrow(map_nodes) > 0) {
    print(map_nodes[, c('id', 'type', 'label', 'gene_name')])
  } else {
    cat("No map nodes found in parsed data\n")
  }
} else {
  cat('No pathway data returned\n')
}

cat("\n=== DEBUG COMPLETE ===\n")
