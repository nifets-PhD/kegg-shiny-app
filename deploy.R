# Deployment script for KEGG Shiny App
# This script prepares the app for deployment on shinyapps.io with renv

library(rsconnect)

# Load required packages for deployment
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
library(jsonlite)

# Ensure renv is activated for deployment
if (file.exists("renv/activate.R")) {
  cat("Activating renv environment...\n")
  source("renv/activate.R")
} else {
  cat("Warning: renv/activate.R not found. Ensure renv is properly set up.\n")
}

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

# Function to sync renv environment before deployment
sync_renv <- function() {
  cat("Syncing renv environment...\n")
  
  if (!requireNamespace("renv", quietly = TRUE)) {
    cat("❌ renv package not found. Install with: install.packages('renv')\n")
    return(FALSE)
  }
  
  tryCatch({
    # Restore packages from lockfile
    cat("Restoring packages from renv.lock...\n")
    renv::restore(prompt = FALSE)
    
    # Update snapshot if needed (since you're in explicit mode)
    cat("Updating renv snapshot from DESCRIPTION...\n")
    renv::snapshot(prompt = FALSE)
    
    cat("✅ renv environment synced successfully\n")
    return(TRUE)
  }, error = function(e) {
    cat("❌ Failed to sync renv environment:", e$message, "\n")
    return(FALSE)
  })
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
  
  # Check renv setup
  if (!file.exists("renv.lock")) {
    cat("⚠️ renv.lock file not found. Run renv::snapshot() first.\n")
    return(FALSE)
  } else {
    cat("✅ renv.lock file found\n")
  }
  
  # Check if renv is in explicit mode
  if (file.exists("renv/settings.json")) {
    settings <- jsonlite::fromJSON("renv/settings.json")
    if (settings$snapshot.type == "explicit") {
      cat("✅ renv configured in explicit mode (using DESCRIPTION)\n")
    } else {
      cat("ℹ️ renv in", settings$snapshot.type, "mode\n")
    }
  }
  
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

# Function to terminate any stuck deployments
terminate_deployment <- function(app_name = "evo-kegg-pathways", account = NULL) {
  cat("Attempting to terminate deployment for:", app_name, "\n")
  
  tryCatch({
    rsconnect::terminateApp(appName = app_name, account = account)
    cat("✅ Successfully terminated deployment\n")
    cat("Wait a few minutes before attempting to deploy again\n")
    return(TRUE)
  }, error = function(e) {
    cat("❌ Failed to terminate deployment:", e$message, "\n")
    return(FALSE)
  })
}

# Main deployment function
deploy_app <- function(app_name = "evo-kegg-pathways", 
                      pre_cache = FALSE,  # Changed default to FALSE to avoid fgsea
                      sync_environment = TRUE,  # New parameter for renv sync
                      account = NULL) {
  
  cat("🚀 Deploying Evo KEGG Pathways with renv...\n")
  
  # Check requirements first
  if (!check_requirements()) {
    cat("❌ Deployment cancelled due to missing requirements.\n")
    return(FALSE)
  }
  
  # Sync renv environment if requested
  if (sync_environment) {
    if (!sync_renv()) {
      cat("❌ Failed to sync renv environment. Deployment cancelled.\n")
      return(FALSE)
    }
  }
  
  # Pre-populate cache if requested
  if (pre_cache) {
    cat("Pre-populating cache for deployment...\n")
    pre_populate_cache()
  }
  
  # Prepare for deployment
  prepare_for_deployment(include_cache = TRUE)
  
  cat("Deploying to shinyapps.io with renv...\n")
  cat("App name:", app_name, "\n")
  
  # Deploy the app with additional options to handle build issues
  tryCatch({
    # First, try to terminate any existing deployments
    cat("Checking for existing deployments...\n")
    tryCatch({
      apps <- rsconnect::applications(account = account)
      existing_app <- apps[apps$name == app_name, ]
      if (nrow(existing_app) > 0 && !is.na(existing_app$status) && existing_app$status == "deploying") {
        cat("Found deployment in progress, attempting to terminate...\n")
        rsconnect::terminateApp(appName = app_name, account = account)
        Sys.sleep(5)  # Wait a moment for termination
      }
    }, error = function(e) {
      cat("Note: Could not check/terminate existing deployment:", e$message, "\n")
    })
    
    # Now attempt deployment with retry logic
    max_attempts <- 3
    for (attempt in 1:max_attempts) {
      cat("Deployment attempt", attempt, "of", max_attempts, "...\n")
      
      deploy_result <- tryCatch({
        # Deploy with renv support
        rsconnect::deployApp(
          appName = app_name,
          appTitle = "Evo KEGG Pathway Visualizer",
          account = account,
          launch.browser = (attempt == max_attempts), # Only launch browser on final attempt
          forceUpdate = TRUE,
          logLevel = "normal",
          lint = FALSE,  # Skip linting to avoid potential issues
          # renv will be automatically detected and used
          appFiles = c(
            "app.R", "ui.R", "server.R", "global.R",
            "DESCRIPTION", "renv.lock", "renv/",
            "R/", "data/", "cache/", "kegg_cache/"
          )
        )
        TRUE
      }, error = function(e) {
        if (grepl("409", e$message) && attempt < max_attempts) {
          cat("HTTP 409 conflict detected, waiting before retry...\n")
          Sys.sleep(10 * attempt)  # Exponential backoff
          return(FALSE)
        } else {
          stop(e)
        }
      })
      
      if (deploy_result) break
    }
    
    cat("\n🎉 Deployment with renv complete!\n")
    cat("Your app should be available at: https://", 
        if(is.null(account)) "[your-account]" else account, 
        ".shinyapps.io/", app_name, "/\n")
    cat("📦 All dependencies managed by renv from DESCRIPTION file\n")
    
    return(TRUE)
    
  }, error = function(e) {
    cat("❌ Deployment failed:", e$message, "\n")
    cat("You may need to check your renv.lock file or try again.\n")
    return(FALSE)
  })
}

# Main deployment interface
if (interactive()) {
  cat("🚀 Evo KEGG Shiny App Deployment with renv\n")
  cat("==========================================\n")
  cat("Platform: shinyapps.io\n")
  cat("Environment Management: renv (explicit mode)\n\n")
  
  cat("Available functions:\n")
  cat("1. check_requirements() - Verify setup and renv status\n")
  cat("2. sync_renv() - Sync renv environment with DESCRIPTION\n")
  cat("3. deploy_app('app-name') - Deploy the application with renv\n")
  cat("4. pre_populate_cache() - Download essential pathways\n")
  cat("5. terminate_deployment('app-name') - Stop stuck deployments\n\n")
  
  cat("Quick deploy workflow:\n")
  cat("1. sync_renv()  # Ensure environment is up to date\n")
  cat("2. deploy_app('kegg-pathway-viz')  # Deploy with renv\n\n")
  
  cat("💡 renv Status:\n")
  if (file.exists("renv.lock")) {
    cat("✅ renv.lock found - dependencies locked\n")
  } else {
    cat("⚠️ renv.lock missing - run renv::snapshot()\n")
  }
  
  if (file.exists("DESCRIPTION")) {
    cat("✅ DESCRIPTION file found - explicit mode enabled\n")
  }
  
  cat("\n📦 Dependencies managed by renv from DESCRIPTION file\n")
  cat("🔄 Bioconductor repositories are properly configured.\n")
}



