# Gene Mapping Utilities
# Simple utilities for loading cached gene mapping data

#' Load comprehensive gene mapping from cached file
#' @return data.frame with gene mapping information or NULL if not available
load_comprehensive_gene_mapping <- function() {
    mapping_file <- "data/gene_id_mapping.rds"
    
    if (file.exists(mapping_file)) {
        tryCatch({
            mapping_data <- readRDS(mapping_file)
            cat("Loaded comprehensive gene mapping:", nrow(mapping_data), "entries\n")
            return(mapping_data)
        }, error = function(e) {
            cat("Error loading gene mapping:", e$message, "\n")
            return(NULL)
        })
    } else {
        cat("Warning: Gene mapping file not found at", mapping_file, "\n")
        return(NULL)
    }
}

#' Get UniProt ID for HGNC gene symbol using mapping data
#' @param hgnc_symbol HGNC gene symbol
#' @param mapping_data cached mapping data
#' @return UniProt ID or NULL if not found
get_uniprot_id <- function(hgnc_symbol, mapping_data = NULL) {
    if (is.null(hgnc_symbol) || is.na(hgnc_symbol) || hgnc_symbol == "") {
        return(NULL)
    }
    
    if (is.null(mapping_data)) {
        mapping_data <- load_comprehensive_gene_mapping()
        if (is.null(mapping_data)) return(NULL)
    }
    
    # Find the UniProt ID for this gene
    tryCatch({
        match_idx <- which(mapping_data$hgnc_symbol == hgnc_symbol)
        
        if (length(match_idx) > 0) {
            # Look for UniProt ID in uniprotswissprot column
            uniprot_id <- mapping_data$uniprotswissprot[match_idx[1]]
            if (!is.null(uniprot_id) && !is.na(uniprot_id) && uniprot_id != "") {
                return(uniprot_id)
            }
        }
        
        return(NULL)
    }, error = function(e) {
        return(NULL)
    })
}

#' Create HTML for UniProt structure link
#' @param uniprot_id UniProt ID
#' @return HTML string for AlphaFold structure link
create_uniprot_structure_link <- function(uniprot_id) {
    if (is.null(uniprot_id) || is.na(uniprot_id) || uniprot_id == "") {
        return("")
    }
    
    uniprot_url <- paste0("https://www.uniprot.org/uniprotkb/", uniprot_id, "/entry#structure")
    return(paste0('<a href="', uniprot_url, '" target="_blank" style="color: #007bff;">View 3D Structure ↗</a>'))
}
