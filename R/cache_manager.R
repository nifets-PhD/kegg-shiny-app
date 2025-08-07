# Enhanced Server-Side Caching System for KEGG Pathways
# This implements a multi-level caching strategy that works across user sessions

#' Enhanced KEGG Pathway Cache Manager
#' Provides server-side caching that persists across user sessions
#' and implements intelligent cache warming strategies

# Global cache configuration
CACHE_CONFIG <- list(
    kegg_cache_dir = "kegg_cache",
    shared_cache_dir = "cache",
    max_cache_age_days = 30,
    warm_cache_pathways = c(
        "hsa04110", "hsa04151", "hsa04010", "hsa04012", "hsa04014",  # Common signaling
        "hsa04060", "hsa04150", "hsa04210", "hsa04310", "hsa04330",  # Immune/Development
        "hsa00010", "hsa00020", "hsa00030", "hsa00040", "hsa00051",  # Metabolism
        "hsa05200", "hsa05210", "hsa05212", "hsa05213", "hsa05214"   # Cancer pathways
    ),
    preload_gene_mappings = TRUE
)

#' Initialize shared cache system
#' Creates necessary cache directories and sets up permissions
init_shared_cache <- function() {
    tryCatch({
        # Create cache directories
        for (cache_dir in c(CACHE_CONFIG$kegg_cache_dir, CACHE_CONFIG$shared_cache_dir)) {
            if (!dir.exists(cache_dir)) {
                dir.create(cache_dir, recursive = TRUE)
                cat("Created cache directory:", cache_dir, "\n")
            }
        }
        
        # Set up cache monitoring
        cache_info <- get_cache_stats()
        cat("Cache initialized:\n")
        cat("  KEGG cache:", cache_info$kegg_files, "files\n")
        cat("  Shared cache:", cache_info$shared_files, "files\n")
        cat("  Total size:", round(cache_info$total_size_mb, 2), "MB\n")
        
        return(TRUE)
    }, error = function(e) {
        cat("Error initializing cache:", e$message, "\n")
        return(FALSE)
    })
}

#' Get cache statistics
get_cache_stats <- function() {
    kegg_files <- 0
    shared_files <- 0
    total_size <- 0
    
    if (dir.exists(CACHE_CONFIG$kegg_cache_dir)) {
        kegg_list <- list.files(CACHE_CONFIG$kegg_cache_dir, full.names = TRUE)
        kegg_files <- length(kegg_list)
        if (kegg_files > 0) {
            total_size <- total_size + sum(file.size(kegg_list), na.rm = TRUE)
        }
    }
    
    if (dir.exists(CACHE_CONFIG$shared_cache_dir)) {
        shared_list <- list.files(CACHE_CONFIG$shared_cache_dir, full.names = TRUE)
        shared_files <- length(shared_list)
        if (shared_files > 0) {
            total_size <- total_size + sum(file.size(shared_list), na.rm = TRUE)
        }
    }
    
    list(
        kegg_files = kegg_files,
        shared_files = shared_files,
        total_size_mb = total_size / (1024^2)
    )
}

#' Warm cache with commonly used pathways
#' Pre-loads popular pathways when the server starts
warm_pathway_cache <- function() {
    cat("Warming pathway cache...\n")
    
    for (pathway_id in CACHE_CONFIG$warm_cache_pathways) {
        tryCatch({
            cat("Pre-loading pathway:", pathway_id, "\n")
            
            # Check if already cached
            cache_file <- file.path(CACHE_CONFIG$kegg_cache_dir, paste0(pathway_id, ".rds"))
            
            if (!file.exists(cache_file)) {
                # Load and cache the pathway
                pathway_data <- parse_kegg_pathway_with_hsa(pathway_id, use_cached = TRUE)
                
                if (!is.null(pathway_data$nodes) && nrow(pathway_data$nodes) > 0) {
                    cat("  Successfully cached:", pathway_id, "\n")
                } else {
                    cat("  Failed to load:", pathway_id, "\n")
                }
            } else {
                cat("  Already cached:", pathway_id, "\n")
            }
            
            # Small delay to avoid overwhelming KEGG servers
            Sys.sleep(0.5)
            
        }, error = function(e) {
            cat("  Error caching", pathway_id, ":", e$message, "\n")
        })
    }
    
    cat("Cache warming completed.\n")
}

#' Clean old cache files
#' Removes cache files older than specified days
clean_old_cache <- function(max_age_days = CACHE_CONFIG$max_cache_age_days) {
    cat("Cleaning cache files older than", max_age_days, "days...\n")
    
    cleaned_count <- 0
    total_saved_mb <- 0
    
    for (cache_dir in c(CACHE_CONFIG$kegg_cache_dir, CACHE_CONFIG$shared_cache_dir)) {
        if (dir.exists(cache_dir)) {
            files <- list.files(cache_dir, full.names = TRUE)
            
            for (file in files) {
                if (file.exists(file)) {
                    file_age_days <- as.numeric(Sys.Date() - as.Date(file.mtime(file)))
                    
                    if (file_age_days > max_age_days) {
                        file_size_mb <- file.size(file) / (1024^2)
                        unlink(file)
                        cleaned_count <- cleaned_count + 1
                        total_saved_mb <- total_saved_mb + file_size_mb
                        cat("  Removed old file:", basename(file), "\n")
                    }
                }
            }
        }
    }
    
    cat("Cache cleaning completed:", cleaned_count, "files removed,", 
        round(total_saved_mb, 2), "MB freed\n")
}

#' Enhanced pathway loading with intelligent caching
#' Extends the existing load_kegg_pathway function with better caching logic
load_kegg_pathway_cached <- function(pathway_id, force_refresh = FALSE) {
    # Clean pathway ID
    pathway_id <- gsub("^path:", "", pathway_id)
    
    if (is.na(pathway_id) || pathway_id == "" || pathway_id == "hsaNA") {
        stop(paste("Invalid pathway ID:", pathway_id))
    }
    
    if (!grepl("^hsa\\d+", pathway_id)) {
        stop(paste("Invalid pathway ID format:", pathway_id))
    }
    
    # Check shared cache first
    shared_cache_file <- file.path(CACHE_CONFIG$shared_cache_dir, paste0(pathway_id, "_processed.rds"))
    
    if (!force_refresh && file.exists(shared_cache_file)) {
        tryCatch({
            cat("Loading from shared cache:", pathway_id, "\n")
            cached_data <- readRDS(shared_cache_file)
            
            # Validate cache integrity
            if (is.list(cached_data) && 
                !is.null(cached_data$nodes) && 
                !is.null(cached_data$edges) &&
                nrow(cached_data$nodes) > 0) {
                
                # Update access time for cache management
                file.create(shared_cache_file)  # Touch file
                return(cached_data)
            } else {
                cat("Shared cache corrupted, removing...\n")
                unlink(shared_cache_file)
            }
        }, error = function(e) {
            cat("Error reading shared cache:", e$message, "\n")
            unlink(shared_cache_file)
        })
    }
    
    # Load using existing function
    cat("Loading pathway fresh:", pathway_id, "\n")
    pathway_data <- load_kegg_pathway(pathway_id)
    
    # Save to shared cache if successful
    if (!is.null(pathway_data) && 
        !is.null(pathway_data$nodes) && 
        nrow(pathway_data$nodes) > 0) {
        
        tryCatch({
            saveRDS(pathway_data, shared_cache_file)
            cat("Saved to shared cache:", pathway_id, "\n")
        }, error = function(e) {
            cat("Error saving to shared cache:", e$message, "\n")
        })
    }
    
    return(pathway_data)
}

#' Cache management endpoint for monitoring
get_cache_management_info <- function() {
    stats <- get_cache_stats()
    
    # Get pathway list from cache
    cached_pathways <- character(0)
    if (dir.exists(CACHE_CONFIG$kegg_cache_dir)) {
        rds_files <- list.files(CACHE_CONFIG$kegg_cache_dir, pattern = "\\.rds$")
        cached_pathways <- gsub("\\.rds$", "", rds_files)
    }
    
    # Get recently accessed pathways
    if (dir.exists(CACHE_CONFIG$shared_cache_dir)) {
        shared_files <- list.files(CACHE_CONFIG$shared_cache_dir, pattern = "_processed\\.rds$", full.names = TRUE)
        recent_access <- data.frame(
            pathway = gsub("_processed\\.rds$", "", basename(shared_files)),
            last_accessed = file.mtime(shared_files),
            stringsAsFactors = FALSE
        )
        recent_access <- recent_access[order(recent_access$last_accessed, decreasing = TRUE), ]
    } else {
        recent_access <- data.frame(pathway = character(0), last_accessed = character(0))
    }
    
    list(
        statistics = stats,
        cached_pathways = cached_pathways,
        warm_cache_pathways = CACHE_CONFIG$warm_cache_pathways,
        recent_access = head(recent_access, 20)
    )
}

cat("Enhanced caching system loaded. Use init_shared_cache() to initialize.\n")
