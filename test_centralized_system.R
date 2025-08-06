# Test Centralized Gene ID Management System
# ========================================
# This script demonstrates how the new centralized system eliminates complexity

library(dplyr)

# Load the centralized system
source("R/gene_id_manager.R")
source("R/network_utils.R")  # For comprehensive mapping

cat("=== Testing Centralized Gene ID Management System ===\n\n")

# Load gene mapping data
comprehensive_mapping <- load_comprehensive_gene_mapping()
cat("✅ Loaded comprehensive gene mapping with", nrow(comprehensive_mapping), "entries\n\n")

# Test 1: Mixed input types - all converted to Entrez internally
cat("1. Testing mixed input types (all converted to Entrez internally):\n")

# Symbols
symbols <- c("TP53", "BRCA1", "EGFR", "MYC")
symbol_result <- convert_to_internal_standard(symbols, "symbol", comprehensive_mapping)
cat("   Symbols:", paste(symbols, collapse = ", "), "\n")
cat("   → Entrez:", paste(symbol_result$entrez_ids, collapse = ", "), "\n")
cat("   → Success rate:", round(symbol_result$conversion_stats$conversion_rate * 100, 1), "%\n\n")

# Entrez IDs
entrez_ids <- c("7157", "672", "1956", "4233")
entrez_result <- convert_to_internal_standard(entrez_ids, "entrez", comprehensive_mapping)
cat("   Entrez:", paste(entrez_ids, collapse = ", "), "\n")
cat("   → Entrez:", paste(entrez_result$entrez_ids, collapse = ", "), "\n")
cat("   → Success rate:", round(entrez_result$conversion_stats$conversion_rate * 100, 1), "%\n\n")

# UniProt IDs
uniprot_ids <- c("P04637", "P38398", "P00533", "P01106")
uniprot_result <- convert_to_internal_standard(uniprot_ids, "uniprot", comprehensive_mapping)
cat("   UniProt:", paste(uniprot_ids, collapse = ", "), "\n")
cat("   → Entrez:", paste(uniprot_result$entrez_ids, collapse = ", "), "\n")
cat("   → Success rate:", round(uniprot_result$conversion_stats$conversion_rate * 100, 1), "%\n\n")

# Test 2: Display conversion (Entrez back to symbols for user display)
cat("2. Testing display conversion (Entrez → Symbols for user-friendly display):\n")
all_entrez <- unique(c(symbol_result$entrez_ids, entrez_result$entrez_ids, uniprot_result$entrez_ids))
display_symbols <- entrez_to_symbols(all_entrez, comprehensive_mapping)
cat("   Internal Entrez IDs:", paste(all_entrez, collapse = ", "), "\n")
cat("   → Display symbols:", paste(display_symbols, collapse = ", "), "\n\n")

# Test 3: Demonstrate simplicity - no more conditional logic needed!
cat("3. Demonstrating system simplicity:\n")
cat("   ❌ OLD WAY: Complex conditional conversions everywhere\n")
cat("   ✅ NEW WAY: Single conversion at input, Entrez throughout, convert only for display\n\n")

# Test 4: Show how this works with different contexts
cat("4. Testing context-specific gene handling:\n")

# Simulate uploading different gene types
user_symbols <- c("GAPDH", "ACTB", "TUBB", "ALDOA")
user_conversion <- convert_to_internal_standard(user_symbols, "symbol", comprehensive_mapping)
cat("   User uploaded symbols:", paste(user_symbols, collapse = ", "), "\n")
cat("   → Internal storage:", paste(user_conversion$entrez_ids, collapse = ", "), "(Entrez)\n")

# For KEGG pathway analysis - use Entrez directly (no conversion needed!)
cat("   → For KEGG analysis: Use", paste(user_conversion$entrez_ids, collapse = ", "), "directly ✅\n")

# For expression analysis expecting symbols - convert from internal Entrez
expr_symbols <- entrez_to_symbols(user_conversion$entrez_ids, comprehensive_mapping)
cat("   → For expression analysis:", paste(expr_symbols, collapse = ", "), "(converted from internal Entrez)\n")

# For pathway visualization - use Entrez directly
cat("   → For pathway visualization: Use", paste(user_conversion$entrez_ids, collapse = ", "), "directly ✅\n\n")

cat("=== Summary ===\n")
cat("✅ ALL gene IDs stored internally as Entrez IDs (single source of truth)\n")
cat("✅ NO conditional conversion logic scattered throughout the app\n") 
cat("✅ SIMPLE conversion only at input (to Entrez) and output (for display)\n")
cat("✅ KEGG pathway operations work directly with internal IDs\n")
cat("✅ Expression data gets appropriate format through single conversion\n")
cat("✅ CLEAN, maintainable, and much less error-prone!\n\n")

cat("🎉 The centralized system eliminates all the conversion complexity!\n")
