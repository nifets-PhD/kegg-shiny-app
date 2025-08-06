# Color Configuration for KEGG Shiny App
# Centralized color palette for easy customization



# Custom CSS Colors for Fine Control
CSS_COLORS <- list(
  # Purple variations
  primary_purple = "#332660ff",      # Medium slate blue
  dark_purple3 = "#332660ff",        # Plum
  light_purple = "#332660ff",        # Plum
  dark_purple = "#1c172fff",         # Rebecca purple
  dark_purple2 = "#2a243fff",
  lavender = "#E6E6FA",            # Light lavender
  
  # Green variations  
  primary_green = "#3f5b4c",       # Bootstrap success green
  light_green = "#D4EDDA",         # Light green background
  dark_green = "#769789",          # Dark green text
  mint_green = "#010201",          # Pale green

  primary_red= "#ffba7e",
  dark_red = "#f38a2e",
  
  # Neutral colors
  white = "#FFFFFF",
  light_gray = "#F8F9FA",
  gray = "#6C757D",
  dark_gray = "#343A40"
)

# Phylostratum Color Palette (using existing GLOBAL_PHYLOSTRATA_COLORS)
# This uses the system defined in global.R with PS_colours function from kegg_utils.R
# No need to redefine - just reference the existing system



# Function to get CSS color
get_css_color <- function(color_name) {
  return(CSS_COLORS[[color_name]] %||% CSS_COLORS$gray)
}

# Function to get phylostratum color (uses existing GLOBAL_PHYLOSTRATA_COLORS)
get_phylostratum_color <- function(stratum) {
  stratum_key <- as.character(stratum)
  return(GLOBAL_PHYLOSTRATA_COLORS[stratum_key] %||% "#808080")
}

# Function to generate custom CSS
generate_custom_css <- function() {
  paste0(
    "
    /* Custom Purple and Green Theme */
    .main-header .navbar {
      background-color: ", CSS_COLORS$dark_purple2, " !important;
    }
    
    .main-header .logo {
      background-color: ", CSS_COLORS$dark_purple, " !important;
      color: ", CSS_COLORS$white, " !important;
    }
    
    .main-sidebar {
      background-color: ", CSS_COLORS$dark_purple, " !important;
    }
    
    .sidebar-menu > li > a {
      color: ", CSS_COLORS$white, " !important;
    }
    
    .sidebar-menu > li.active > a {
      background-color: ", CSS_COLORS$dark_purple3, " !important;
      border-left-color: ", CSS_COLORS$light_purple, " !important;
    }
    
    /* Card customizations - ALL primary cards to purple */
    .box.box-solid.box-primary > .box-header {
      background: ", CSS_COLORS$primary_purple, " !important;
      background-color: ", CSS_COLORS$primary_purple, " !important;
    }
    
    .box.box-solid.box-primary {
      border: 1px solid ", CSS_COLORS$primary_purple, " !important;
    }
    
    .box.box-primary > .box-header {
      background: ", CSS_COLORS$primary_purple, " !important;
      background-color: ", CSS_COLORS$primary_purple, " !important;
      color: ", CSS_COLORS$white, " !important;
    }
    
    .box.box-primary {
      border-top-color: ", CSS_COLORS$primary_purple, " !important;
    }
    
    /* Success cards to green */
    .box.box-solid.box-success > .box-header {
      background: ", CSS_COLORS$light_purple, " !important;
      background-color: ", CSS_COLORS$light_purple, " !important;
    }
    
    .box.box-solid.box-success {
      border: 1px solid ", CSS_COLORS$light_purple, " !important;
    }
    
    .box.box-success > .box-header {
      background: ", CSS_COLORS$light_purple, " !important;
      background-color: ", CSS_COLORS$light_purple, " !important;
      color: ", CSS_COLORS$white, " !important;
    }
    
    .box.box-success {
      border-top-color: ", CSS_COLORS$light_purple, " !important;
    }
    
    /* Info cards to light purple instead of cyan */
    .box.box-solid.box-info > .box-header {
      background: ", CSS_COLORS$light_purple, " !important;
      background-color: ", CSS_COLORS$light_purple, " !important;
    }
    
    .box.box-solid.box-info {
      border: 1px solid ", CSS_COLORS$light_purple, " !important;
    }
    
    .box.box-info > .box-header {
      background: ", CSS_COLORS$light_purple, " !important;
      background-color: ", CSS_COLORS$light_purple, " !important;
      color: ", CSS_COLORS$white, " !important;
    }
    
    .box.box-info {
      border-top-color: ", CSS_COLORS$light_purple, " !important;
    }
    
    /* Button customizations */
    .btn-primary {
      background-color: ", CSS_COLORS$light_green, " !important;
      border-color: ", CSS_COLORS$light_green, " !important;
    }
    
    .btn-primary:hover, .btn-primary:focus, .btn-primary:active {
      background-color: ", CSS_COLORS$dark_green, " !important;
      border-color: ", CSS_COLORS$dark_green, " !important;
    }
    
    .btn-success {
      background-color: ", CSS_COLORS$primary_red, " !important;
      border-color: ", CSS_COLORS$primary_red, " !important;
    }
    
    .btn-success:hover, .btn-success:focus, .btn-success:active {
      background-color: ", CSS_COLORS$dark_red, " !important;
      border-color: ", CSS_COLORS$primary_red, " !important;
    }
    
    .btn-info {
      background-color: ", CSS_COLORS$light_green, " !important;
      border-color: ", CSS_COLORS$light_green, " !important;
    }
    
    .btn-info:hover, .btn-info:focus, .btn-info:active {
      background-color: ", CSS_COLORS$dark_green, " !important;
      border-color: ", CSS_COLORS$dark_green, " !important;
    }
    
    /* Alert/notification backgrounds */
    .alert-info {
      background-color: ", CSS_COLORS$lavender, " !important;
      border-color: ", CSS_COLORS$light_purple, " !important;
      color: ", CSS_COLORS$dark_purple, " !important;
    }
    
    .alert-success {
      background-color: ", CSS_COLORS$light_green, " !important;
      border-color: ", CSS_COLORS$primary_green, " !important;
      color: ", CSS_COLORS$dark_green, " !important;
    }
    "
  )
}
