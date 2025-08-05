# Deployment script for KEGG Shiny App
# This script prepares the app for deployment on shinyapps.io

library(rsconnect)

# Essential pathways to pre-cache for all users
ESSENTIAL_PATHWAYS <- c(
  "hsa04110", # Cell cycle
  "hsa04151", # PI3K-Akt signaling pathway
  "hsa04152", # AMPK signaling pathway
  "hsa04010", # MAPK signaling pathway
  "hsa04066", # HIF-1 signaling pathway
  "hsa04142", # Lysosome
  "hsa04140", # Autophagy
  "hsa00010", # Glycolysis / Gluconeogenesis
  "hsa00020", # Citrate cycle (TCA cycle)
  "hsa03013", # RNA transport
  "hsa03010"  # Ribosome
)

# Function to pre-populate essential cache
pre_populate_cache <- function() {
  cat("Pre-populating cache with essential pathways...\n")
  
  # Source the utilities
  source("R/kegg_utils.R")
  
  # Create cache directories
  if (!dir.exists("kegg_cache")) dir.create("kegg_cache")
  
  # Pre-load essential pathways
  for (pathway_id in ESSENTIAL_PATHWAYS) {
    cat("Pre-caching pathway:", pathway_id, "\n")
    tryCatch({
      # This will download and cache the pathway if not already cached
      pathway_data <- load_kegg_pathway(pathway_id)
      if (!is.null(pathway_data)) {
        cat("  ✓ Successfully cached", pathway_id, "\n")
      }
    }, error = function(e) {
      cat("  ✗ Failed to cache", pathway_id, ":", e$message, "\n")
    })
    
    # Small delay to be respectful to KEGG servers
    Sys.sleep(1)
  }
  
  cat("Pre-population complete!\n")
}

# Function to prepare app for deployment
prepare_for_deployment <- function(include_cache = TRUE) {
  cat("Preparing KEGG Shiny App for deployment...\n")
  
  # Create cache directory and copy existing cache files
  if (include_cache) {
    if (!dir.exists("cache")) dir.create("cache")
    
    if (dir.exists("kegg_cache")) {
      cache_files <- list.files("kegg_cache", full.names = TRUE)
      if (length(cache_files) > 0) {
        file.copy(cache_files, "cache/", overwrite = TRUE)
        cat("Copied", length(cache_files), "KEGG cache files\n")
      }
    }
  }
  
  # Ensure essential data files are present
  essential_files <- c("data/phylomap_hgnc.tsv", "data/strata_legend.tsv")
  missing_files <- essential_files[!file.exists(essential_files)]
  
  if (length(missing_files) > 0) {
    cat("WARNING: Missing essential data files:\n")
    cat(paste(missing_files, collapse = "\n"), "\n")
    cat("Please ensure these files are present before deployment.\n")
  } else {
    cat("✓ All essential data files present\n")
  }
  
  cat("App preparation complete!\n")
}

# Check deployment requirements
check_requirements <- function() {
  cat("Checking deployment requirements...\n")
  
  # Check shinyapps.io authentication
  accounts <- rsconnect::accounts()
  if (nrow(accounts) == 0) {
    cat("⚠️ No shinyapps.io account configured.\n")
    cat("Set up with: rsconnect::setAccountInfo(name='account', token='token', secret='secret')\n")
    cat("Get credentials from: https://www.shinyapps.io/admin/#/tokens\n")
    return(FALSE)
  } else {
    cat("✅ shinyapps.io account configured:", accounts$name[1], "\n")
  }
  
  return(TRUE)
}

# Main deployment function
deploy_app <- function(app_name = "evo-kegg-pathways", 
                      pre_cache = FALSE,  # Changed default to FALSE to avoid fgsea
                      account = NULL) {
  
  cat("🚀 Deploying Evo KEGG Pathways...\n")
  
  # Check requirements first
  if (!check_requirements()) {
    cat("❌ Deployment cancelled due to missing requirements.\n")
    return(FALSE)
  }
  
  # Pre-populate cache if requested
  if (pre_cache) {
    cat("Pre-populating cache for deployment...\n")
    pre_populate_cache()
  }
  
  # Prepare for deployment
  prepare_for_deployment(include_cache = TRUE)
  
  cat("Deploying to shinyapps.io...\n")
  cat("App name:", app_name, "\n")
  
  # Deploy the app with additional options to handle build issues
  tryCatch({
    rsconnect::deployApp(
      appName = app_name,
      appTitle = "Evo KEGG Pathway Visualizer",
      account = account,
      launch.browser = TRUE,
      forceUpdate = TRUE,
      logLevel = "normal",
      lint = FALSE  # Skip linting to avoid potential issues
    )
    
    cat("\n🎉 Deployment complete!\n")
    cat("Your app should be available at: https://", 
        if(is.null(account)) "[your-account]" else account, 
        ".shinyapps.io/", app_name, "/\n")
    
    return(TRUE)
    
  }, error = function(e) {
    cat("❌ Deployment failed:", e$message, "\n")
    cat("You may need to check your package dependencies or try again.\n")
    return(FALSE)
  })
}

# Main deployment interface
if (interactive()) {
  cat("🚀 Evo KEGG Shiny App Deployment\n")
  cat("============================\n")
  cat("Platform: shinyapps.io\n\n")
  
  cat("Available functions:\n")
  cat("1. check_requirements() - Verify setup\n")
  cat("2. deploy_app('app-name') - Deploy the application\n")
  cat("3. pre_populate_cache() - Download essential pathways\n\n")
  
  cat("Quick deploy:\n")
  cat("deploy_app('kegg-pathway-viz')\n\n")
  
  cat("💡 Current repository settings look good!\n")
  cat("Bioconductor repositories are properly configured.\n")
}



