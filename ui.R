# Define UI
ui <- dashboardPage(
    dashboardHeader(title = "KEGG Pathway Visualizer"),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem("Pathway Explorer", tabName = "pathways", icon = icon("project-diagram")),
            menuItem("Network Visualization", tabName = "network", icon = icon("share-alt")),
            menuItem("Gene Annotations", tabName = "annotations", icon = icon("tags"))
        )
    ),
    
    dashboardBody(
        tags$head(
            tags$style(HTML("
                .content-wrapper, .right-side {
                    background-color: #f4f4f4;
                }
                .box {
                    margin-bottom: 20px;
                }
            "))
        ),
        
        tabItems(
            # Pathway Explorer Tab
            tabItem(tabName = "pathways",
                fluidRow(
                    box(
                        title = "KEGG Pathway Selection", status = "primary", solidHeader = TRUE,
                        width = 12, collapsible = TRUE,
                        
                        fluidRow(
                            column(6,
                                h4("Search Pathways"),
                                textInput("pathway_search", "Search pathway name:", 
                                         placeholder = "e.g., glycolysis, metabolism"),
                                actionButton("search_pathways", "Search", class = "btn-primary"),
                                br(), br(),
                                h5("Or select from categories:"),
                                selectInput("pathway_category", "Pathway Category:",
                                           choices = list(
                                               "Metabolism" = "01100",
                                               "Genetic Information Processing" = "03000", 
                                               "Environmental Information Processing" = "04000",
                                               "Cellular Processes" = "05000",
                                               "Organismal Systems" = "06000",
                                               "Human Diseases" = "07000"
                                           ),
                                           selected = "01100"),
                                br(),
                                div(style = "font-size: 12px; color: #666;",
                                    "Pathways are cached locally to reduce server load")
                            ),
                            column(6,
                                h4("Available Pathways"),
                                DT::dataTableOutput("pathway_table"),
                                br(),
                                actionButton("load_pathway", "Load Selected Pathway", 
                                           class = "btn-success", icon = icon("download"))
                            )
                        )
                    )
                ),
                
                fluidRow(
                    box(
                        title = "Pathway Information", status = "info", solidHeader = TRUE,
                        width = 12, collapsible = TRUE,
                        
                        verbatimTextOutput("pathway_info")
                    )
                )
            ),
            
            # Network Visualization Tab
            tabItem(tabName = "network",
                fluidRow(
                    box(
                        title = "Pathway Network", status = "success", solidHeader = TRUE,
                        width = 8, collapsible = TRUE,
                        
                        visNetworkOutput("pathway_network", height = "700px")
                    ),
                    
                    box(
                        title = "Network Options", status = "primary", solidHeader = TRUE,
                        width = 4, collapsible = TRUE,
                        
                        h5("Node Coloring"),
                        radioButtons("node_coloring", "Color nodes by:",
                                   choices = list(
                                       "KEGG Default" = "kegg_default",
                                       "Phylostratum" = "phylostratum"
                                   ),
                                   selected = "kegg_default"),
                        
                        conditionalPanel(
                            condition = "input.node_coloring == 'phylostratum'",
                            helpText("Colors nodes by evolutionary age (phylostratum)."),
                            br()
                        ),
                        
                        # Phylostratum Legend
                        conditionalPanel(
                            condition = "input.node_coloring == 'phylostratum'",
                            hr(),
                            h5("Phylostratum Legend"),
                            div(style = "max-height: 400px; overflow-y: auto;",
                                DT::dataTableOutput("phylostratum_legend_table")
                            )
                        )
                    )
                ),
                
                fluidRow(
                    box(
                        title = "Edge Legend", status = "warning", solidHeader = TRUE,
                        width = 6, collapsible = TRUE, height = "500px",
                        
                        h5("Relationship Types"),
                        p("Different colors and styles represent different biological relationships:"),
                        div(style = "max-height: 400px; overflow-y: auto;",
                            DT::dataTableOutput("edge_legend_table")
                        )
                    ),
                    
                    box(
                        title = "Selected Node Information", status = "info", solidHeader = TRUE,
                        width = 6, collapsible = TRUE,
                        
                        verbatimTextOutput("node_info")
                    )
                )
            ),
            
            # Gene Annotations Tab
            tabItem(tabName = "annotations",
                fluidRow(
                    box(
                        title = "Gene Annotation Upload", status = "primary", solidHeader = TRUE,
                        width = 6, collapsible = TRUE,
                        
                        h4("Upload Annotation File"),
                        fileInput("annotation_file", "Choose CSV/TSV File:",
                                 accept = c(".csv", ".tsv", ".txt")),
                        
                        helpText("File should contain columns: gene_id, annotation_value"),
                        
                        h4("Annotation Settings"),
                        selectInput("annotation_column", "Annotation Column:",
                                   choices = NULL),
                        
                        selectInput("color_scheme", "Color Scheme:",
                                   choices = list(
                                       "Viridis" = "viridis",
                                       "Plasma" = "plasma",
                                       "Red-Blue" = "RdBu",
                                       "Red-Yellow-Blue" = "RdYlBu",
                                       "Spectral" = "Spectral"
                                   ),
                                   selected = "viridis"),
                        
                        actionButton("apply_annotations", "Apply Annotations", 
                                    class = "btn-success", icon = icon("paint-brush"))
                    ),
                    
                    box(
                        title = "Annotation Preview", status = "info", solidHeader = TRUE,
                        width = 6, collapsible = TRUE,
                        
                        DT::dataTableOutput("annotation_preview")
                    )
                ),
                
                fluidRow(
                    box(
                        title = "Color Legend", status = "success", solidHeader = TRUE,
                        width = 12, collapsible = TRUE,
                        
                        plotlyOutput("color_legend", height = "200px")
                    )
                )
            )
        )
    )
)
