# KEGG Shiny App Improvement Plan

## Overview
This document outlines the priority improvements made to address gene ID handling, phylomap integration, and pathway visualization issues in the KEGG Shiny app.

## Priority 1: Enhanced Gene ID Type Support ✅ COMPLETED

### Problem
- App only supported HGNC gene symbols as input
- KEGG pathways are protein-centric and work better with other ID types
- Limited gene conversion capabilities

### Solutions Implemented
1. **Multi-ID Type Support**: Added support for:
   - **Entrez Gene ID** (recommended for KEGG) - numeric IDs like 1956, 7157
   - **Gene Symbol (HGNC)** - traditional symbols like TP53, BRCA1  
   - **Ensembl Gene ID** - stable IDs like ENSG00000141510
   - **UniProt ID** - protein IDs like P04637

2. **UI Updates**:
   - Added gene ID type selector dropdown
   - Dynamic instruction panels that update based on selected ID type
   - Default examples appropriate for each ID type
   - Better user guidance with explanatory text

3. **Backend Updates**:
   - Updated `convert_ids_to_entrez()` function to handle all ID types
   - Enhanced `perform_kegg_enrichment()` to accept `id_type` parameter
   - Improved `validate_gene_ids()` for comprehensive validation
   - Server-side tracking of current gene ID type

### Benefits
- Better KEGG pathway mapping (Entrez IDs are optimal)
- More flexible input options for users
- Improved validation and error handling
- Clearer user feedback on ID conversion rates

## Priority 2: Comprehensive Phylomap Integration ✅ COMPLETED

### Problem  
- App used limited `phylomap_hgnc.tsv` with only gene symbols
- New comprehensive `phylomap.tsv` contains richer data with Ensembl IDs and protein information
- Missing integration of evolutionary data

### Solutions Implemented
1. **Data Source Update**:
   - Switched from `data/phylomap_hgnc.tsv` to `data/phylomap.tsv`
   - Updated `load_phylomap()` function to use comprehensive data
   - Enhanced data structure with columns: ENSP, HGNC, ENSG, Stratum, Name, etc.

2. **Phylostratum Mapping Enhancement**:
   - Updated `map_genes_to_phylostrata()` to use HGNC column from new data
   - Added case-insensitive matching for better gene symbol mapping
   - Improved phylostratum coloring with comprehensive data support

3. **Server Integration**:
   - Updated server to load comprehensive phylomap on startup
   - Enhanced gene table display to show phylomap matches
   - Better integration with gene validation workflows

### Benefits
- More genes covered in phylostratum analysis
- Richer evolutionary data integration
- Better gene-to-stratum mapping accuracy
- Future-ready for Ensembl ID integration

## Priority 3: Enhanced XML Pathway Parsing ✅ IN PROGRESS

### Problem
- XML parsing wasn't capturing all node types (missing metabolites, compounds)
- Some edges not being displayed properly
- Limited selection capabilities for non-gene nodes

### Solutions Implemented
1. **Comprehensive Node Type Support**:
   - Enhanced XML parsing to capture ALL entry types:
     - **Genes**: Traditional gene entries with symbols
     - **Compounds**: Metabolites and chemical compounds  
     - **Maps**: Links to other KEGG pathways
     - **Groups**: Protein complexes and functional groups
     - **Orthologs**: Orthologous gene groups
   
2. **Improved Node Processing**:
   - Type-specific label generation and formatting
   - Better description generation for each node type
   - Enhanced graphics attribute extraction
   - Proper handling of compound and metabolite names

3. **Better Edge Extraction**:
   - Enhanced relationship parsing from XML `<relation>` elements
   - Capture of subtypes and relationship details
   - Improved edge styling based on relationship types

### Benefits
- Complete pathway visualization (not just genes)
- Better representation of metabolic pathways
- Improved user interaction with all node types
- More accurate biological pathway representation

## Implementation Status

### ✅ Completed
- [x] Multi-gene ID type support (UI + backend)
- [x] Comprehensive phylomap integration  
- [x] Enhanced XML node type parsing
- [x] Updated validation and error handling
- [x] Improved user interface with dynamic instructions
- [x] Server-side gene ID type tracking

### 🔄 In Progress  
- [ ] Enhanced edge relationship visualization
- [ ] Improved node selection for non-gene types
- [ ] Better metabolite/compound information display

### 📋 Future Enhancements
- [ ] Direct Ensembl ID pathway mapping
- [ ] UniProt protein information integration
- [ ] Advanced pathway filtering options
- [ ] Export capabilities for comprehensive results

## Usage Recommendations

### For Users
1. **Prefer Entrez Gene IDs**: Best compatibility with KEGG pathways
2. **Use comprehensive data**: New phylomap provides better coverage
3. **Explore all node types**: Pathways now show genes, compounds, and complexes

### For Developers  
1. **Test with different ID types**: Ensure all conversion paths work
2. **Validate phylostratum coloring**: Check with comprehensive data
3. **Verify XML parsing**: Test with various pathway types

## Files Modified

### Core Functions
- `R/enrichment_utils.R` - Enhanced gene ID handling
- `R/kegg_utils.R` - Improved XML parsing and phylomap integration
- `ui.R` - Dynamic gene ID type selection
- `server.R` - Comprehensive backend updates

### Data Files  
- Using `data/phylomap.tsv` (comprehensive) instead of `data/phylomap_hgnc.tsv`
- `data/strata_legend.tsv` for phylostratum legends

## Testing Checklist

- [x] Test all 4 gene ID types (Entrez, Symbol, Ensembl, UniProt)
- [x] **FIXED: Gene highlighting now works with all ID types via Entrez conversion**
- [x] Verify phylostratum coloring with comprehensive data
- [x] Test pathway loading with different pathway types
- [x] **IMPROVED: Gene conversion rates displayed in UI**
- [x] Check node selection for all node types (genes, compounds, maps)
- [x] **ENHANCED: Gene matching in pathways uses Entrez IDs for reliability**
- [x] Verify edge display and styling

## Latest Fix: Universal Gene ID Support for Network Visualization ✅ COMPLETED

### Problem Identified

- Users inputting UniProt IDs (or other non-Entrez ID types) could not see their genes highlighted in pathway networks
- Gene matching was failing because KEGG pathways internally use Entrez IDs, but highlighting was attempted using original input IDs
- Network visualization showed "none of the genes were found" even when genes were actually present in the pathway

### Solution Implemented

1. **Early ID Conversion System**:
   - All input gene IDs are now converted to Entrez IDs immediately after loading
   - Both original IDs and converted Entrez IDs are stored in reactive values
   - `converted_entrez_genes()` reactive value stores Entrez IDs for KEGG operations

2. **Enhanced Gene Highlighting**:
   - Updated `apply_gene_highlighting()` function to prioritize Entrez ID matching
   - Gene highlighting now uses converted Entrez IDs instead of original input IDs
   - Added multi-level matching: Entrez IDs → HGNC symbols → node labels (fallback)

3. **Improved Gene Pathway Matching**:
   - `genes_in_pathway` output now uses Entrez ID matching for reliable gene detection
   - Shows which genes were found via Entrez ID conversion vs traditional matching
   - Better error messages when genes aren't found, with conversion rate information

4. **Enhanced User Feedback**:
   - Gene statistics display shows conversion rates for non-Entrez ID inputs
   - Validation results clearly indicate how many genes are "KEGG-ready"
   - Instructions updated to mention automatic conversion to Entrez IDs

5. **Optimized Enrichment Analysis**:
   - KEGG enrichment now directly uses converted Entrez IDs for better performance
   - Eliminates redundant ID conversion during enrichment analysis

## Latest Enhancement: Comprehensive Selected Gene Card Display ✅ COMPLETED

### Issues Addressed

1. **Gene Symbol Display Issue**: Selected gene cards were showing Entrez IDs (numbers) instead of readable gene symbols, even though network nodes displayed symbols
2. **Reduced Functionality**: The comprehensive gene information display from the original implementation was simplified and missing key features

### New Features Implemented

1. **Comprehensive Gene ID Display**:
   - Shows all available gene ID types (Symbol, Entrez, Ensembl, UniProt, KEGG)
   - Dynamically pulls from the comprehensive gene mapping database
   - Fallback handling when comprehensive mapping is not available
   - Clear section organization with "GENE IDENTIFIERS" header

2. **Enhanced UniProt Information**:
   - Displays protein names and aliases when available
   - Shows comprehensive mapping information beyond just IDs
   - Maintains 3D structure links for protein visualization
   - Better error handling for missing UniProt data

3. **Improved Relationship Display**:
   - Relationship connections now prefer gene symbols over IDs for readability
   - Smart fallback: Symbol → Label → ID for connection display
   - Maintains all relationship types (incoming/outgoing connections)
   - Shows relationship subtypes and detailed connection information

4. **User-Friendly Display Priority**:
   - Primary display uses gene symbols instead of Entrez IDs when available
   - Comprehensive ID mapping shown in organized sections
   - Clear visual hierarchy with distinct sections
   - Maintains backward compatibility with existing functionality

### Technical Implementation

- **Enhanced Node Selection Handler**: Complete rewrite of the `observe` block for pathway network selection
- **Comprehensive Mapping Integration**: Uses the `comprehensive_mapping` reactive value loaded at startup
- **Smart ID Resolution**: Automatically resolves the best display names for each gene
- **Error Handling**: Graceful fallbacks when mapping data is unavailable
- **Performance Optimized**: Efficient lookups using pre-loaded mapping data

### User Experience Improvements

- ✅ **Readable Gene Names**: Selected gene cards now show symbols (e.g., "TP53") instead of numbers (e.g., "7157")
- ✅ **Complete ID Coverage**: Display all relevant gene ID types in one place
- ✅ **Comprehensive Information**: Restored full functionality with enhanced data sources
- ✅ **Better Navigation**: Clear sections and organized information hierarchy
- ✅ **Consistent Experience**: Network labels and selected gene cards now use consistent naming

### Technical Files Updated

- `server.R`: Enhanced node selection observer with comprehensive gene ID display and improved relationship handling

## Bug Fix: Node Selection Error Resolution ✅ FIXED

### Issue Details

- Error when clicking on gene nodes in pathway networks: "missing value where TRUE/FALSE needed"
- Error occurred in `renderUI` function at line 817 in server.R
- Issue was caused by `!is.na()` checks on potentially NULL values

### Root Cause

- When accessing columns that don't exist in data frames (e.g., `uniprot_info$name` or `gene_record$alias_symbol`), R returns `NULL`
- Using `!is.na(NULL)` causes the error because `is.na()` expects a vector, not NULL
- The comprehensive gene mapping data structure varies, and some columns might not exist

### Fix Applied

- **Enhanced NULL Checking**: Added `!is.null()` checks before all `!is.na()` conditions
- **Safe Column Access**: All data frame column access now follows the pattern: `!is.null(df$col) && !is.na(df$col) && df$col != ""`
- **Comprehensive Error Prevention**: Fixed all instances in gene ID display, UniProt information, and relationship display sections

### Fixed Code Sections

1. **Gene ID Information Display**: Safe access to `hgnc_symbol`, `entrezgene_id`, `ensembl_gene_id`, `uniprotswissprot`
2. **UniProt Information Section**: Safe access to `name`, `alias_symbol`, and other protein attributes
3. **Relationship Display**: Safe access to `hgnc_symbol` and `label` in source/target nodes
4. **Fallback Information**: Safe handling of gene symbol display

### Technical Details

**Before (causing error):**

```r
if (!is.na(gene_record$name) && gene_record$name != "")
```

**After (error-free):**

```r
if (!is.null(gene_record$name) && !is.na(gene_record$name) && gene_record$name != "")
```

### Benefits for Users

Users can now click on any gene node in pathway networks and see:

- **All Gene ID Types**: Symbol, Entrez, Ensembl, UniProt, and KEGG IDs in one place
- **Readable Names**: Gene symbols prominently displayed instead of numeric IDs
- **Complete Information**: Protein information, aliases, and comprehensive mapping data
- **Clear Relationships**: Connection information using readable gene names
- **Consistent Interface**: Gene names match between network display and information cards

### New Benefits

- ✅ **Universal ID Support**: Users can input any supported ID type and see genes highlighted
- ✅ **Reliable Matching**: Uses Entrez IDs for consistent gene detection in KEGG pathways  
- ✅ **Better UX**: Clear feedback about conversion success rates and gene availability
- ✅ **Performance**: Eliminates redundant ID conversions during analysis
- ✅ **Backward Compatibility**: Still supports traditional gene symbol matching as fallback

### Additional Files Modified

- `server.R`: Added `converted_entrez_genes` reactive, enhanced gene loading, improved UI feedback
- `R/kegg_utils.R`: Updated `apply_gene_highlighting()` for Entrez ID priority matching
- `R/enrichment_utils.R`: Already had robust ID conversion functions (no changes needed)

### User Impact

Users can now confidently upload **any supported gene ID type** (Entrez, Symbol, Ensembl, UniProt) and:

- See immediate feedback about conversion success rates
- Have their genes properly highlighted in pathway networks
- Get reliable gene matching regardless of input ID type
- Receive clear explanations when genes cannot be converted or found
