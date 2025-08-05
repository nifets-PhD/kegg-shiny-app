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
  cat("Pre-populating cache with essential pathways...
")
  
  # Source the utilities
  source("R/kegg_utils.R")
  
  # Create cache directories
  if (!dir.exists("kegg_cache")) dir.create("kegg_cache")
  if (!dir.exists("hsa_cache")) dir.create("hsa_cache")
  
  # Pre-load essential pathways
  for (pathway_id in ESSENTIAL_PATHWAYS) {
    cat("Pre-caching pathway:", pathway_id, "
")
    tryCatch({
      # This will download and cache the pathway if not already cached
      pathway_data <- load_kegg_pathway(pathway_id)
      if (!is.null(pathway_data)) {
        cat("  ✓ Successfully cached", pathway_id, "
")
      }
    }, error = function(e) {
      cat("  ✗ Failed to cache", pathway_id, ":", e$message, "
")
    })
    
    # Small delay to be respectful to KEGG servers
    Sys.sleep(1)
  }
  
  cat("Pre-population complete!
")
}

# Function to prepare app for deployment
prepare_for_deployment <- function(include_cache = TRUE) {
  cat("Preparing KEGG Shiny App for deployment...
")
  
  # Create cache directory and copy existing cache files
  if (include_cache) {
    if (!dir.exists("cache")) dir.create("cache")
    
    if (dir.exists("kegg_cache")) {
      cache_files <- list.files("kegg_cache", full.names = TRUE)
      if (length(cache_files) > 0) {
        file.copy(cache_files, "cache/", overwrite = TRUE)
        cat("Copied", length(cache_files), "KEGG cache files
")
      }
    }
    
    if (dir.exists("hsa_cache")) {
      hsa_files <- list.files("hsa_cache", full.names = TRUE)
      if (length(hsa_files) > 0) {
        file.copy(hsa_files, "cache/", overwrite = TRUE)
        cat("Copied", length(hsa_files), "HSA cache files
")
      }
    }
  }
  
  # Ensure essential data files are present
  essential_files <- c("data/phylomap_hgnc.tsv", "data/strata_legend.tsv")
  missing_files <- essential_files[!file.exists(essential_files)]
  
  if (length(missing_files) > 0) {
    cat("WARNING: Missing essential data files:
")
    cat(paste(missing_files, collapse = "
"), "
")
  } else {
    cat("✓ All essential data files present
")
  }
  
  cat("App preparation complete!
")
}

# Simple deployment function
deploy_app <- function(app_name = "evo-kegg-pathway-viz", pre_cache = FALSE, account = NULL) {
  cat("🚀 Deploying Evo KEGG Pathway Visualizer...
")
  
  # Pre-populate cache if requested
  if (pre_cache) {
    pre_populate_cache()
  }
  
  # Prepare for deployment
  prepare_for_deployment(include_cache = TRUE)
  
  cat("Deploying to shinyapps.io...
")
  
  # Deploy the app
  rsconnect::deployApp(
    appName = app_name,
    appTitle = "Evo KEGG Pathway Visualizer - Biological Network Explorer",
    account = account,
    launch.browser = TRUE,
    forceUpdate = TRUE,
    logLevel = "normal"
  )
  
  cat("🎉 Deployment complete!
")
  cat("Your app should be available at: https://", 
      if(is.null(account)) "[your-account]" else account, 
      ".shinyapps.io/", app_name, "/
")
}

# Check deployment requirements
check_requirements <- function() {
  cat("Checking deployment requirements...
")
  
  # Check shinyapps.io authentication
  accounts <- rsconnect::accounts()
  if (nrow(accounts) == 0) {
    cat("⚠️ No shinyapps.io account configured.
")
    cat("Set up with: rsconnect::setAccountInfo(name='account', token='token', secret='secret')
")
    cat("Get credentials from: https://www.shinyapps.io/admin/#/tokens
")
    return(FALSE)
  } else {
    cat("✅ shinyapps.io account configured:", accounts$name[1], "
")
  }
  
  return(TRUE)
}

# Main deployment interface
if (interactive()) {
  cat("🚀 KEGG Shiny App Deployment
")
  cat("============================
")
  cat("Platform: shinyapps.io

")
  
  cat("Available functions:
")
  cat("1. check_requirements() - Verify setup
")
  cat("2. deploy_app('app-name') - Deploy the application
")
  cat("3. pre_populate_cache() - Download essential pathways

")
  
  cat("Quick deploy:
")
  cat("deploy_app('kegg-pathway-viz')

")
  
  cat("💡 Note: Bioconductor packages may cause deployment issues.
")
  cat("If deployment fails, try deploying without pre-caching or
")
  cat("consider using a different deployment platform.
")
}

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
  if (!dir.exists("hsa_cache")) dir.create("hsa_cache")
  
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
  
  # Create necessary directories
  if (!dir.exists("cache")) {
    dir.create("cache", recursive = TRUE)
    cat("Created cache directory\n")
  }
  
  if (!dir.exists("data")) {
    dir.create("data", recursive = TRUE)
    cat("Created data directory\n")
  }
  
  # Copy existing cache files to deployment cache
  if (include_cache) {
    if (dir.exists("kegg_cache")) {
      cache_files <- list.files("kegg_cache", full.names = TRUE)
      if (length(cache_files) > 0) {
        file.copy(cache_files, "cache/", overwrite = TRUE)
        cat("Copied", length(cache_files), "KEGG cache files\n")
      }
    }
    
    if (dir.exists("hsa_cache")) {
      hsa_files <- list.files("hsa_cache", full.names = TRUE)
      if (length(hsa_files) > 0) {
        file.copy(hsa_files, "cache/", overwrite = TRUE)
        cat("Copied", length(hsa_files), "HSA cache files\n")
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

# Function to create a simplified manifest to avoid Bioconductor issues
create_simple_manifest <- function() {
  cat("Creating simplified manifest.json...\n")
  
  # Create a basic manifest with core packages only
  manifest <- list(
    version = 1,
    metadata = list(
      appmode = "shiny",
      entrypoint = "app"
    ),
    packages = list(
      list(name = "shiny", version = as.character(packageVersion("shiny"))),
      list(name = "shinydashboard", version = as.character(packageVersion("shinydashboard"))),
      list(name = "visNetwork", version = as.character(packageVersion("visNetwork"))),
      list(name = "DT", version = as.character(packageVersion("DT"))),
      list(name = "xml2", version = as.character(packageVersion("xml2"))),
      list(name = "dplyr", version = as.character(packageVersion("dplyr"))),
      list(name = "readr", version = as.character(packageVersion("readr")))
    ),
    files = list(
      list(checksum = "", relativePath = "app.R"),
      list(checksum = "", relativePath = "ui.R"),
      list(checksum = "", relativePath = "server.R"),
      list(checksum = "", relativePath = "global.R")
    ),
    R = list(
      Version = paste(R.Version()$major, R.Version()$minor, sep = "."),
      Repository = "https://cran.rstudio.com/"
    )
  )
  
  # Write manifest
  jsonlite::write_json(manifest, "manifest.json", pretty = TRUE, auto_unbox = TRUE)
  cat("Created simplified manifest.json\n")
}

# Function to deploy to shinyapps.io with properly fixed packages
deploy_to_shinyapps <- function(app_name = "kegg-pathway-visualizer", 
                               pre_cache = TRUE, 
                               account = NULL, 
                               ignore_bioc_version = TRUE) {
  
  cat("🚀 Deploying with renv validation bypass...
")
  
  # Pre-populate cache if requested
  if (pre_cache) {
    cat("Pre-populating cache for deployment...
")
    pre_populate_cache()
  }
  
  # Prepare for deployment
  prepare_for_deployment(include_cache = TRUE)
  
  cat("Deploying to shinyapps.io...
")
  cat("App name:", app_name, "
")
  
  # Set environment variables to bypass renv validation
  old_options <- options()
  on.exit(options(old_options))
  
  # Configure options to bypass problematic validation
  options(
    renv.config.snapshot.validate = FALSE,
    renv.config.bioconductor.repos = character(),
    repos = c(CRAN = "https://cran.rstudio.com/")
  )
  
  # Create .Renviron to persist these settings during deployment
  renviron_content <- paste(
    "RENV_CONFIG_SNAPSHOT_VALIDATE=FALSE",
    "RENV_CONFIG_BIOCONDUCTOR_REPOS=",
    "R_REPOS=https://cran.rstudio.com/",
    sep = "
"
  )
  writeLines(renviron_content, ".Renviron")
  cat("Created .Renviron with validation bypass settings
")
  
  # Deploy with bypassed validation
  tryCatch({
    rsconnect::deployApp(
      appName = app_name,
      appTitle = "KEGG Pathway Visualizer - Interactive Biological Network Explorer",
      account = account,
      launch.browser = TRUE,
      forceUpdate = TRUE,
      logLevel = "normal"
    )
  }, error = function(e) {
    cat("Standard deployment failed:", e$message, "
")
    cat("Trying with additional bypass options...
")
    
    # Try with even more aggressive bypassing
    rsconnect::deployApp(
      appName = app_name,
      appTitle = "KEGG Pathway Visualizer",
      account = account,
      launch.browser = TRUE,
      forceUpdate = TRUE,
      logLevel = "quiet",
      lint = FALSE
    )
  })
  
  # Clean up
  if (file.exists(".Renviron")) unlink(".Renviron")
  
  cat("
🎉 Deployment complete!
")
  cat("Your app should be available at: https://", 
      if(is.null(account)) "[your-account]" else account, 
      ".shinyapps.io/", app_name, "/
")
}


