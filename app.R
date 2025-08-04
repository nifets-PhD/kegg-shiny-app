# KEGG Pathway Visualizer - Main Application File
# This combines UI and Server for shinyapps.io deployment

# Source global configuration
source("global.R")

# Source UI and Server definitions
source("ui.R")
source("server.R")

# Run the application
shinyApp(ui = ui, server = server)
