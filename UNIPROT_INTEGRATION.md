# UniProt Structure Integration

## Overview
The KEGG Shiny app now includes direct integration with UniProt for accessing protein structure information. When you click on a gene in the network visualization, the gene information panel will show a link to the UniProt structure page for that protein.

## Features
- **Automatic HGNC to UniProt Mapping**: Pre-cached mapping of ~19,800 human genes to their UniProt IDs
- **Direct Structure Links**: Click "Go to UniProt Structure" to access detailed protein structure information
- **Swiss-Prot Preference**: Prefers high-quality Swiss-Prot entries over TrEMBL when available
- **Cached Performance**: Fast lookup using local cache, no network requests during app usage

## Implementation Details

### Files Created
- `data/hgnc_uniprot_mapping.rds` - Cached mapping data (binary format)
- `data/hgnc_uniprot_mapping.tsv` - Human-readable mapping file
- `download_uniprot_mapping.R` - Initial download script
- `refresh_uniprot_mapping.R` - Refresh script for periodic updates

### Functions Added
- `download_hgnc_uniprot_mapping()` - Downloads complete mapping from Ensembl
- `get_uniprot_id()` - Fast lookup using cached data
- `create_uniprot_structure_link()` - Creates HTML link to UniProt structure page

### Usage
1. Select any gene node in the network visualization
2. The gene information panel now shows comprehensive details including:
   - **Basic Gene Info**: KEGG ID, HGNC symbol, node type, display label
   - **Phylostratum**: Evolutionary age information (when available)
   - **UniProt Details**: UniProt ID, SwissProt/TrEMBL IDs, Ensembl ID
   - **Description**: Clean gene description from UniProt/Ensembl
   - **Structure Link**: Direct link to UniProt structure page
   - **Relationships**: Network connections and interaction types
3. Click "Go to UniProt Structure" to open UniProt in a new tab
4. The UniProt page will show detailed protein information including:
   - 3D structure data (if available)
   - Structural annotations
   - Domain information
   - Experimental structures from PDB

## Maintenance
- The mapping cache is valid for 30 days
- Run `Rscript refresh_uniprot_mapping.R` to update the cache
- The initial download contains mappings for 19,803 genes with UniProt IDs

## Example Links
- TP53 → https://www.uniprot.org/uniprotkb/E7EQX7/entry#structure
- EGFR → https://www.uniprot.org/uniprotkb/P00533/entry#structure
- GAPDH → https://www.uniprot.org/uniprotkb/P04406/entry#structure
