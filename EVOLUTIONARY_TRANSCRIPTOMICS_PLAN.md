# Evolutionary Transcriptomics Tab Implementation Plan

## Overview
Add a 5th tab called "Evolutionary Transcriptomics" to the KEGG Shiny app that integrates myTAI functionality for analyzing gene expression datasets with evolutionary information.

## Project Status: 🚧 PLANNING PHASE

### ✅ Completed Tasks
- [x] **Project Planning**: Created comprehensive implementation plan
- [x] **myTAI Analysis**: Reviewed myTAI package structure and main functions
- [x] **Current App Analysis**: Analyzed existing KEGG Shiny app structure
- [x] **Integration Point Identification**: Identified how to integrate with existing phylomap data
- [x] **Function Scope Refinement**: Removed `plot_strata_expression`, `plot_distribution_expression`, and separate statistical tests (included in `plot_signature`)
- [x] **Example Dataset**: Added `data/example_expression.tsv` with HGNC gene symbols and 15 cell type columns
- [x] **Function List Finalization**: Updated to include all 7 myTAI functions including `plot_sample_space`

### 🔄 Current Phase: Implementation - Phase 1 (Backend Functions)

### 📋 TODO: Implementation Tasks

## Feature Requirements

### 1. Data Upload & Processing
- **File Upload**: Support CSV/TSV format gene expression data
  - First column: Gene IDs (support for 4 types: entrez, symbol, ensembl, uniprot)
  - Remaining columns: Sample expression values (bulk transcriptomics)
  - Validation for proper format and gene ID types
  
- **Gene ID Mapping**: 
  - Reuse existing `map_genes_to_phylostrata()` function from kegg_utils.R
  - Support the same 4 ID types as the current app
  - Map user genes to phylostratum information using existing phylomap.tsv

### 2. Core myTAI Integration
- **PhyloExpressionSet Creation**: Convert uploaded data to myTAI format
- **TAI Calculation**: Compute Transcriptome Age Index for samples
- **Statistical Tests**: Implement myTAI statistical functions

### 3. Visualization Functions (from myTAI)

- **`plot_signature()`**: Main TAI signature plot across samples (includes statistical tests)
- **`plot_distribution_strata()`**: Distribution of genes across phylostrata
- **`plot_gene_heatmap()`**: Heatmap of selected genes with phylostratum coloring
- **`plot_contribution()`**: Contribution to overall TAI by phylostratum
- **`plot_gene_space()`**: Gene space analysis
- **`plot_sample_space()`**: Sample space analysis
- **`plot_gene_profiles()`**: Individual gene expression profiles with interaction coloring

### 4. Gene Selection & Highlighting
- **Dropdown for Gene Selection**: 
  - Option 1: "Selected Gene Set" (from uploaded genes)
  - Option 2: "Genes in Current Pathway" (from KEGG pathway browser)
  
- **Interactive Gene Selection**:
  - Select individual genes from pathway network
  - Show interacting genes with color-coding by interaction type
  - Use `plot_gene_profiles()` with interaction-based coloring

### 5. Integration with Existing Features
- **Phylomap Integration**: Use existing `data/phylomap.tsv` 
- **Gene Set Integration**: Connect with gene sets from other tabs
- **Pathway Integration**: Use genes from selected KEGG pathways
- **ID Type Consistency**: Maintain same gene ID type support as other tabs

## Technical Implementation Plan

### Phase 1: Backend Functions (R/evolutionary_utils.R)
```r
# New utility functions to create:
- create_phyloexpression_set()    # Convert uploaded data to myTAI format
- validate_expression_data()      # Validate uploaded expression data
- map_expression_to_phylostrata() # Map expression data to evolutionary info
- plot_tai_signature_custom()     # Custom TAI plotting with highlighting
- plot_gene_profiles_interactions() # Gene profiles colored by interactions
- get_interacting_genes()         # Extract genes that interact with selected gene
```

### Phase 2: UI Components (ui.R)
```r
# Add new tabItem for "evolutionary":
- File upload widget for expression data
- Gene ID type selector (reuse existing)
- Sample selection controls
- Plot type selector (signature, distribution, heatmap, etc.)
- Gene highlighting controls
- Download buttons for plots and results
```

### Phase 3: Server Logic (server.R)  
```r
# New server components:
- Expression data upload handling
- myTAI PhyloExpressionSet creation
- Plot generation with myTAI functions
- Gene selection and highlighting logic
- Integration with existing pathway data
```

### Phase 4: Testing & Refinement
- Test with various expression datasets
- Validate phylostratum mapping accuracy
- Ensure proper integration with existing tabs
- Performance optimization for large datasets

## Data Structure Requirements

### Input Data Format (Example: data/example_expression.tsv)
```tsv
GeneID          Unannotated  CP-like  Cilia-like  IPC-like  Neurons  ...
OR4F29          0            0        0           0         0        ...
SAMD11          134          32       5           21        1        ...
NOC2L           296          226      50          42        26       ...
KLHL17          54           81       17          12        29       ...
...
```

### myTAI PhyloExpressionSet Format
```tsv
Phylostratum    GeneID      Unannotated  CP-like  Cilia-like  IPC-like  Neurons  ...
1               OR4F29      0            0        0           0         0        ...
3               SAMD11      134          32       5           21        1        ...
2               NOC2L       296          226      50          42        26       ...
1               KLHL17      54           81       17          12        29       ...
...
```

## User Interface Mockup

### Layout Structure
```
Tab: "Evolutionary Transcriptomics"
├── Row 1: Data Upload & Configuration
│   ├── Box 1: Expression Data Upload
│   │   ├── Gene ID type selector
│   │   ├── File upload widget
│   │   └── Data validation results
│   └── Box 2: Analysis Configuration  
│       ├── Plot type selector
│       ├── Gene highlighting options
│       └── Statistical test options
├── Row 2: Main Visualization
│   └── Box 3: TAI Analysis Plot (full width)
│       ├── plot_signature() output
│       └── Plot controls and options
├── Row 3: Secondary Analyses
│   ├── Box 4: Phylostratum Distribution
│   │   └── plot_distribution_strata()
│   └── Box 5: Gene Selection & Profiles
│       ├── Interactive gene selector
│       └── plot_gene_profiles() with interaction coloring
└── Row 4: Results & Export
    ├── Box 6: Statistical Results
    └── Box 7: Data Export Options
```

## Integration Points with Existing App

### 1. Phylomap Data
- **File**: `data/phylomap.tsv` (already exists)
- **Function**: `load_phylomap()` and `map_genes_to_phylostrata()` (already implemented)
- **Usage**: Direct integration for evolutionary information

### 2. Gene ID Mapping
- **Functions**: `convert_kegg_to_hgnc()` and related mapping functions
- **Usage**: Consistent gene ID handling across all tabs

### 3. Gene Set Integration  
- **Variables**: `uploaded_genes` reactive values
- **Usage**: Allow highlighting of user-uploaded gene sets in evolutionary plots

### 4. Pathway Integration
- **Variables**: `values$nodes` from pathway visualization
- **Usage**: Extract genes from currently loaded pathway for analysis

## Key myTAI Functions to Implement

### myTAI Visualization Functions
1. **`plot_signature()`** - Main TAI signature across samples (includes statistical tests)
2. **`plot_distribution_strata()`** - Gene distribution by phylostrata  
3. **`plot_gene_heatmap()`** - Expression heatmap with evolutionary coloring
4. **`plot_contribution()`** - Phylostratum contribution to TAI
5. **`plot_gene_profiles()`** - Individual gene expression profiles
6. **`plot_gene_space()`** - Gene space visualization
7. **`plot_sample_space()`** - Sample space visualization

## Dependencies & Requirements

### R Packages
- **myTAI**: Main evolutionary transcriptomics package
- **ggplot2**: Already used, for plot customization
- **DT**: Already used, for data tables
- **shiny/shinydashboard**: Already used
- **plotly**: Already used, for interactive plots

### Data Requirements  
- **phylomap.tsv**: Already available
- **Gene ID mapping**: Already implemented
- **Expression data**: User-provided (CSV/TSV)
- **Example dataset**: `data/example_expression.tsv` with HGNC gene symbols and 15 cell type/condition columns

## Risk Assessment & Mitigation

### Technical Risks
1. **myTAI Integration Complexity**
   - *Risk*: myTAI might have complex data structure requirements
   - *Mitigation*: Start with simple examples, gradually add complexity

2. **Performance with Large Datasets**  
   - *Risk*: Large expression matrices might be slow to process
   - *Mitigation*: Implement data filtering and sampling options

3. **Gene ID Mapping Issues**
   - *Risk*: Some genes might not map to phylostrata
   - *Mitigation*: Provide clear feedback on mapping success rates

### User Experience Risks
1. **Complex Interface**
   - *Risk*: Too many options might confuse users
   - *Mitigation*: Use progressive disclosure, start with basic options

2. **Data Format Confusion**
   - *Risk*: Users might upload incorrectly formatted data
   - *Mitigation*: Provide clear examples and validation feedback

## Success Criteria

### Functional Requirements
- [ ] Successfully upload and process gene expression data
- [ ] Generate main myTAI plots (signature, distribution, heatmap)
- [ ] Integrate with existing gene sets and pathway data  
- [ ] Provide gene highlighting and interaction analysis
- [ ] Export results and plots

### Quality Requirements
- [ ] Handle datasets with 10k+ genes and 10+ samples
- [ ] Maintain consistent UI/UX with existing tabs
- [ ] Provide helpful error messages and validation
- [ ] Generate publication-quality plots

### Integration Requirements  
- [ ] Seamless integration with existing phylomap data
- [ ] Compatible with all 4 supported gene ID types
- [ ] Connect with gene sets from other tabs
- [ ] Use genes from selected KEGG pathways

## Timeline Estimate

- **Phase 1 (Backend)**: 3-4 days
- **Phase 2 (UI)**: 2-3 days  
- **Phase 3 (Server)**: 3-4 days
- **Phase 4 (Testing)**: 2-3 days
- **Total**: ~10-14 days

## Next Steps

1. **Set up myTAI development environment**
2. **Create evolutionary_utils.R with basic functions**  
3. **Implement simple PhyloExpressionSet creation**
4. **Add basic UI components to ui.R**
5. **Implement file upload and validation logic**
6. **Add first myTAI plot integration**
7. **Iterate and expand functionality**

---

*Last Updated: August 6, 2025*
*Status: Planning Phase - Ready to begin implementation*
