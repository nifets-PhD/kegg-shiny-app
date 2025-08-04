# KEGG Pathway Visualizer

An interactive Shiny application for exploring KEGG biological pathways with gene highlighting and phylostratum analysis.

## 🌟 Features

- **Interactive Pathway Visualization**: Explore KEGG pathways with authentic layouts and styling
- **Gene-Based Pathway Search**: Find pathways containing your genes using real KEGG API queries
- **Smart Gene Highlighting**: Upload your gene lists and see them highlighted in pathway networks
- **Phylostratum Coloring**: Color genes by evolutionary age (requires myTAI package)
- **Progress Notifications**: Real-time feedback during pathway searches
- **Comprehensive Caching**: Optimized performance with local pathway caching

## 🚀 Live Demo

**Deployed App**: https://[your-account].shinyapps.io/kegg-pathway-visualizer/

## 🛠️ Local Installation

### Prerequisites
```r
# Install required packages
install.packages(c("shiny", "shinydashboard", "visNetwork", "DT", 
                  "xml2", "dplyr", "readr", "rsconnect"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager")
BiocManager::install(c("KEGGREST", "KEGGgraph", "biomaRt"))

# Install myTAI for phylostratum analysis
install.packages("myTAI")
```

### Run Locally
```r
# Clone the repository
git clone https://github.com/nifets-PhD/kegg-shiny-app.git
cd kegg-shiny-app

# Start the app
R -e "shiny::runApp()"
```

## 🌐 Deployment to shinyapps.io

### Step 1: Setup shinyapps.io Account
1. Create account at https://www.shinyapps.io/
2. Get your tokens from https://www.shinyapps.io/admin/#/tokens
3. Configure in R:
```r
rsconnect::setAccountInfo(
  name = 'your-account-name',
  token = 'your-token', 
  secret = 'your-secret'
)
```

### Step 2: Deploy with Enhanced Caching
```r
# Source deployment script
source("deploy.R")

# Check requirements
check_deployment_requirements()

# Deploy with pre-cached pathways (recommended)
deploy_to_shinyapps("kegg-pathway-visualizer", pre_cache = TRUE)
```

### Step 3: Automatic Deployment (Optional)
Set up GitHub Actions for automatic deployment:
```r
create_github_actions_workflow()
```

Then add these secrets to your GitHub repository settings:
- `SHINYAPPS_NAME`: Your shinyapps.io account name
- `SHINYAPPS_TOKEN`: Your token
- `SHINYAPPS_SECRET`: Your secret

## 📊 Global Caching Strategy

The app uses a smart caching strategy to provide fast user experience:

### Pre-Cached Pathways
Essential pathways are pre-loaded during deployment:
- Cell cycle (hsa04110)
- PI3K-Akt signaling (hsa04151) 
- AMPK signaling (hsa04152)
- MAPK signaling (hsa04010)
- HIF-1 signaling (hsa04066)
- Glycolysis (hsa00010)
- TCA cycle (hsa00020)
- And more...

### Dynamic Caching
- New pathways are cached on first load
- Cache persists across user sessions
- Respectful API usage with rate limiting

## 🧬 Usage Guide

### 1. Gene Input
- **Text Input**: One gene symbol per line (e.g., TP53, BRCA1, EGFR)
- **File Upload**: TXT/CSV files with newline-separated gene symbols
- **Format**: HGNC gene symbols (uppercase recommended)

### 2. Pathway Search
- **By Category**: Browse metabolism, signaling, cellular processes
- **By Genes**: Find pathways containing your specific genes
- **By Name**: Search pathway names and descriptions

### 3. Visualization Options
- **KEGG Default**: Authentic KEGG pathway colors and layout
- **Phylostratum**: Color by evolutionary age
- **Gene Highlighting**: Highlight your uploaded genes in red

## 🔧 Technical Details

### Architecture
- **Frontend**: Shiny with shinydashboard
- **Visualization**: visNetwork for interactive networks
- **Data Source**: KEGG API via KEGGREST
- **Gene Mapping**: biomaRt for ID conversion
- **Caching**: Local RDS files for performance

### Performance Optimizations
- Intelligent caching system
- Rate-limited API calls
- Progress notifications
- Early termination for large searches
- Pre-populated essential pathways

## 📁 File Structure
```
kegg-shiny-app/
├── app.R              # Main application entry
├── ui.R               # User interface definition
├── server.R           # Server logic
├── global.R           # Global variables and setup
├── deploy.R           # Deployment helper script
├── R/                 # Utility functions
│   ├── kegg_utils.R   # KEGG pathway processing
│   ├── hsa_utils.R    # HSA gene data utilities
│   └── ...
├── data/              # Essential data files
│   ├── phylomap_hgnc.tsv
│   └── strata_legend.tsv
├── kegg_cache/        # Local KEGG pathway cache
├── hsa_cache/         # Local HSA gene cache
└── .github/workflows/ # GitHub Actions (optional)
```

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Test locally
5. Submit a pull request

## 📄 License

This project is licensed under the MIT License.

## 🙏 Acknowledgments

- **KEGG**: Kyoto Encyclopedia of Genes and Genomes
- **Bioconductor**: For biological data packages
- **RStudio/Posit**: For Shiny framework and shinyapps.io hosting
- **myTAI**: For phylostratum analysis

## 📞 Support

- **Issues**: GitHub Issues page
- **Documentation**: This README and in-app help text
- **KEGG Data**: https://www.genome.jp/kegg/

---

**Made with ❤️ for the biological research community**
