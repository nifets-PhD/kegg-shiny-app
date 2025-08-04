#!/usr/bin/env Rscript

# Test script to debug MTOR gene alias extraction
cat("Testing MTOR gene alias extraction...\n")

# Test XML parsing directly
library(xml2)

xml_file <- "kegg_cache/hsa04152.xml"
if (file.exists(xml_file)) {
    cat("XML file exists\n")
    
    xml_doc <- read_xml(xml_file)
    entries <- xml_find_all(xml_doc, "//entry[@id='92']")
    
    if (length(entries) > 0) {
        cat("Found entry 92\n")
        graphics <- xml_find_first(entries[[1]], "graphics")
        if (!is.null(graphics)) {
            name_attr <- xml_attr(graphics, "name")
            cat("Graphics name attribute:", name_attr, "\n")
            
            # Test the extraction logic
            if (!is.na(name_attr)) {
                first_alias <- strsplit(name_attr, ",", fixed = TRUE)[[1]][1]
                first_alias <- trimws(first_alias)
                first_alias <- gsub("\\.\\.\\.$", "", first_alias)
                cat("Extracted first alias:", first_alias, "\n")
            }
        } else {
            cat("No graphics element found\n")
        }
    } else {
        cat("Entry 92 not found\n")
    }
} else {
    cat("XML file does not exist:", xml_file, "\n")
}
