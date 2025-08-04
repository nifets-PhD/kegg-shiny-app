#!/usr/bin/env Rscript

# Test the complete MTOR gene processing pipeline
cat("Testing MTOR gene processing pipeline...\n")

library(xml2)

# Test the coordinate extraction first
xml_file <- "kegg_cache/hsa04152.xml"
if (file.exists(xml_file)) {
    xml_doc <- read_xml(xml_file)
    entries <- xml_find_all(xml_doc, "//entry[@id='92']")
    
    if (length(entries) > 0) {
        entry <- entries[[1]]
        entry_id <- xml_attr(entry, "id")
        entry_name <- xml_attr(entry, "name")
        entry_type <- xml_attr(entry, "type")
        
        graphics <- xml_find_first(entry, "graphics")
        if (!is.null(graphics)) {
            name_attr <- xml_attr(graphics, "name")
            cat("Raw graphics name:", name_attr, "\n")
            
            # Apply the extraction logic
            if (!is.na(name_attr) && entry_type == "gene") {
                first_alias <- strsplit(name_attr, ",", fixed = TRUE)[[1]][1]
                first_alias <- trimws(first_alias)
                first_alias <- gsub("\\.\\.\\.$", "", first_alias)
                cat("Extracted gene symbol:", first_alias, "\n")
                
                # Show what would be stored
                cat("Entry ID:", entry_id, "\n")
                cat("Entry name (KEGG IDs):", entry_name, "\n")
                cat("Entry type:", entry_type, "\n")
                cat("Final label:", first_alias, "\n")
            }
        }
    }
}
