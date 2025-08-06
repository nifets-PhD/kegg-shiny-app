#!/usr/bin/env Rscript
# Test script to check gene ID mapping functionality

cat("=== Gene ID Mapping Test ===\n\n")

# Source the required functions
source("R/network_utils.R")

# Test 1: Check if the mapping file loads correctly
cat("1. Testing comprehensive gene mapping load...\n")
mapping_data <- load_comprehensive_gene_mapping()

if (is.null(mapping_data)) {
    cat("❌ Failed to load mapping data\n")
    quit(status = 1)
} else {
    cat("✅ Mapping data loaded successfully\n")
    cat("   - Rows:", nrow(mapping_data), "\n")
    cat("   - Columns:", paste(names(mapping_data), collapse = ", "), "\n")
}

# Test 2: Check column names and structure
cat("\n2. Checking mapping data structure...\n")
required_cols <- c("entrezgene_id", "hgnc_symbol")
missing_cols <- required_cols[!required_cols %in% names(mapping_data)]

if (length(missing_cols) > 0) {
    cat("❌ Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
    cat("Available columns:", paste(names(mapping_data), collapse = ", "), "\n")
} else {
    cat("✅ Required columns found\n")
}

# Test 3: Test conversion from Entrez to Symbol
cat("\n3. Testing Entrez to Symbol conversion...\n")
test_entrez_ids <- c("1", "2", "3", "4", "5", "10", "100")  # Common Entrez IDs

# Test the conversion logic used in server.R
if (!is.null(mapping_data) && "entrezgene_id" %in% names(mapping_data) && "hgnc_symbol" %in% names(mapping_data)) {
    mapping_df <- mapping_data[mapping_data$entrezgene_id %in% test_entrez_ids, ]
    converted_symbols <- unique(mapping_df$hgnc_symbol[!is.na(mapping_df$hgnc_symbol)])
    
    cat("   Input Entrez IDs:", length(test_entrez_ids), "\n")
    cat("   Found mappings:", nrow(mapping_df), "\n")
    cat("   Converted symbols:", length(converted_symbols), "\n")
    
    if (length(converted_symbols) > 0) {
        cat("   Examples:", paste(head(converted_symbols, 5), collapse = ", "), "\n")
        cat("✅ Entrez to Symbol conversion working\n")
    } else {
        cat("❌ No symbols converted from Entrez IDs\n")
    }
} else {
    cat("❌ Cannot test conversion - missing columns\n")
}

# Test 4: Test conversion from Symbol to Entrez
cat("\n4. Testing Symbol to Entrez conversion...\n")
test_symbols <- c("A1BG", "A2M", "ACADM", "ACE", "ACTN4")  # Common gene symbols

if (!is.null(mapping_data) && "entrezgene_id" %in% names(mapping_data) && "hgnc_symbol" %in% names(mapping_data)) {
    mapping_df <- mapping_data[toupper(mapping_data$hgnc_symbol) %in% toupper(test_symbols), ]
    converted_entrez <- unique(mapping_df$entrezgene_id[!is.na(mapping_df$entrezgene_id)])
    
    cat("   Input symbols:", length(test_symbols), "\n")
    cat("   Found mappings:", nrow(mapping_df), "\n")
    cat("   Converted Entrez IDs:", length(converted_entrez), "\n")
    
    if (length(converted_entrez) > 0) {
        cat("   Examples:", paste(head(converted_entrez, 5), collapse = ", "), "\n")
        cat("✅ Symbol to Entrez conversion working\n")
    } else {
        cat("❌ No Entrez IDs converted from symbols\n")
    }
} else {
    cat("❌ Cannot test conversion - missing columns\n")
}

# Test 5: Examine actual data structure
cat("\n5. Sample of mapping data:\n")
if (!is.null(mapping_data) && nrow(mapping_data) > 0) {
    sample_rows <- head(mapping_data, 5)
    print(sample_rows)
    
    # Check for empty or NA values
    if ("entrezgene_id" %in% names(mapping_data)) {
        na_entrez <- sum(is.na(mapping_data$entrezgene_id))
        empty_entrez <- sum(mapping_data$entrezgene_id == "", na.rm = TRUE)
        cat("\n   entrezgene_id column: ", na_entrez, "NAs,", empty_entrez, "empty strings\n")
    }
    
    if ("hgnc_symbol" %in% names(mapping_data)) {
        na_symbol <- sum(is.na(mapping_data$hgnc_symbol))
        empty_symbol <- sum(mapping_data$hgnc_symbol == "", na.rm = TRUE)
        cat("   hgnc_symbol column: ", na_symbol, "NAs,", empty_symbol, "empty strings\n")
    }
}

# Test 6: Test with user's specific Entrez IDs (if available)
cat("\n6. Testing with realistic gene set...\n")
# These are common genes that should have mappings
realistic_entrez <- c("7157", "4193", "5728", "2068", "1956", "6790", "1499", "5925")  # p53, MDM2, PTEN, ERCC2, EGFR, AURKA, CTSK, RB1

if (!is.null(mapping_data) && "entrezgene_id" %in% names(mapping_data) && "hgnc_symbol" %in% names(mapping_data)) {
    mapping_df <- mapping_data[mapping_data$entrezgene_id %in% realistic_entrez, ]
    converted_symbols <- unique(mapping_df$hgnc_symbol[!is.na(mapping_df$hgnc_symbol)])
    
    cat("   Realistic Entrez IDs tested:", length(realistic_entrez), "\n")
    cat("   Successfully converted:", length(converted_symbols), "\n")
    
    if (length(converted_symbols) > 0) {
        cat("   Converted symbols:", paste(converted_symbols, collapse = ", "), "\n")
        cat("✅ Realistic gene conversion working\n")
    } else {
        cat("❌ No realistic genes could be converted\n")
        
        # Debug: Show what's in the entrezgene_id column
        cat("\nDEBUG: Sample entrezgene_id values:\n")
        sample_entrez <- head(unique(mapping_data$entrezgene_id[!is.na(mapping_data$entrezgene_id)]), 10)
        cat("   ", paste(sample_entrez, collapse = ", "), "\n")
        cat("   Data types - entrezgene_id:", class(mapping_data$entrezgene_id), "\n")
        cat("   Data types - Test IDs:", class(realistic_entrez), "\n")
    }
}

# Test 7: Test UniProt to Symbol conversion
cat("\n7. Testing UniProt to Symbol conversion...\n")
test_uniprot_ids <- c("P04637", "P38398", "P00533", "P01106", "P04406")  # p53, BRCA1, EGFR, MYC, GAPDH

if (!is.null(mapping_data) && "uniprotswissprot" %in% names(mapping_data) && "hgnc_symbol" %in% names(mapping_data)) {
    mapping_df <- mapping_data[mapping_data$uniprotswissprot %in% test_uniprot_ids, ]
    converted_symbols <- unique(mapping_df$hgnc_symbol[!is.na(mapping_df$hgnc_symbol)])
    
    cat("   Input UniProt IDs:", length(test_uniprot_ids), "\n")
    cat("   Found mappings:", nrow(mapping_df), "\n")
    cat("   Converted symbols:", length(converted_symbols), "\n")
    
    if (length(converted_symbols) > 0) {
        cat("   Examples:", paste(head(converted_symbols, 5), collapse = ", "), "\n")
        cat("✅ UniProt to Symbol conversion working\n")
    } else {
        cat("❌ No symbols converted from UniProt IDs\n")
    }
} else {
    cat("❌ Cannot test UniProt conversion - missing columns\n")
}

cat("\n=== Test Complete ===\n")
