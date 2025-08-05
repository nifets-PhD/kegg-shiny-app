# Define UI
ui <- dashboardPage(
    dashboardHeader(title = "KEGG Pathway Visualizer"),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem("Gene Set", tabName = "geneset", icon = icon("dna")),
            menuItem("KEGG Enrichment", tabName = "enrichment", icon = icon("chart-bar")),
            menuItem("Pathway Explorer", tabName = "explorer", icon = icon("search")),
            menuItem("Network Visualization", tabName = "network", icon = icon("project-diagram"))
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
            tabItem(tabName = "explorer",
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
                                
                                h4("Gene-Based Pathway Search"),
                                tabsetPanel(
                                    tabPanel("Upload File",
                                        br(),
                                        fileInput("gene_file", "Upload Gene List (TXT/CSV):",
                                                 accept = c(".txt", ".csv", ".tsv")),
                                        helpText("REQUIRED: One gene symbol per line only (HGNC format, e.g., TP53, BRCA1)"),
                                        div(style = "color: #d9534f; font-size: 11px;",
                                            "⚠️ Files must contain newline-separated gene symbols only!")
                                    ),
                                    tabPanel("Text Input",
                                        br(),
                                        textAreaInput("gene_text", "Enter Gene Symbols (ONE PER LINE):",
                                                     value = "TP53\nBRCA1\nEGFR\nMYC\nGAPDH\nMTOR\nAKT1\nPIK3CA\nPTEN\nTSC1\nINS\nIGF1\nIRS1\nFOXO1\nGSK3B\nVEGFA\nHIF1A\nMAPK1\nRAF1\nRAS",
                                                     placeholder = "REQUIRED: Enter gene symbols separated by newlines ONLY!\nTP53\nBRCA1\nEGFR\nMYC\nGAPDH",
                                                     height = "120px"),
                                        br(),
                                        actionButton("load_genes_text", "Load Genes", class = "btn-info")
                                    )
                                ),
                                
                                conditionalPanel(
                                    condition = "output.genes_loaded",
                                    hr(),
                                    h5("Find Pathways Containing Your Genes"),
                                    actionButton("find_pathways_with_genes", "Find Pathways", 
                                               class = "btn-success", icon = icon("search")),
                                    br(), br()
                                ),
                                
                                h5("Or select from categories:"),
                                selectInput("pathway_category", "Pathway Category:",
                                           choices = list(
                                               "Metabolism" = "01100",
                                               "Genetic Information Processing" = "03000", 
                                               "Environmental Information Processing" = "04000",
                                               "Cellular Processes" = "05000"
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
                    conditionalPanel(
                        condition = "output.genes_loaded",
                        box(
                            title = "Loaded Gene List", status = "info", solidHeader = TRUE,
                            width = 12, collapsible = TRUE, collapsed = TRUE,
                            
                            fluidRow(
                                column(8,
                                    DT::dataTableOutput("loaded_genes_table")
                                ),
                                column(4,
                                    h5("Gene Statistics"),
                                    verbatimTextOutput("gene_stats"),
                                    br(),
                                    actionButton("clear_genes", "Clear Gene List", 
                                               class = "btn-warning", icon = icon("trash"))
                                )
                            )
                        )
                    )
                ),
                
                fluidRow(
                    conditionalPanel(
                        condition = "output.genes_loaded",
                        box(
                            title = "Your Selected Genes", status = "success", solidHeader = TRUE,
                            width = 12, collapsible = TRUE, collapsed = FALSE,
                            
                            div(style = "max-height: 200px; overflow-y: auto; background: #f8f9fa; padding: 10px; border-radius: 5px;",
                                verbatimTextOutput("selected_genes_display", placeholder = TRUE)
                            ),
                            br(),
                            p(style = "font-size: 12px; color: #666;",
                              "These genes will be highlighted in pathway networks when you load a pathway.")
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
                                       "Phylostratum" = "phylostratum",
                                       "Highlight Uploaded Genes" = "highlight_genes"
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
                ),
                
                # Gene information boxes at the bottom
                fluidRow(
                    conditionalPanel(
                        condition = "output.genes_loaded",
                        box(
                            title = "Your Selected Genes", status = "success", solidHeader = TRUE,
                            width = 6, collapsible = TRUE, collapsed = TRUE,
                            
                            div(style = "max-height: 150px; overflow-y: auto; background: #f8f9fa; padding: 10px; border-radius: 5px;",
                                verbatimTextOutput("selected_genes_display_network", placeholder = TRUE)
                            ),
                            p(style = "font-size: 12px; color: #666; margin-top: 10px;",
                              "⭐ Your genes are highlighted in red in the network visualization above.")
                        ),
                        
                        box(
                            title = "Genes Found in Current Pathway", status = "warning", solidHeader = TRUE,
                            width = 6, collapsible = TRUE, collapsed = FALSE,
                            
                            conditionalPanel(
                                condition = "output.pathway_loaded",
                                div(style = "max-height: 150px; overflow-y: auto; background: #fff3cd; padding: 10px; border-radius: 5px;",
                                    verbatimTextOutput("genes_in_pathway", placeholder = TRUE)
                                ),
                                p(style = "font-size: 12px; color: #856404; margin-top: 10px;",
                                  "🎯 These are your genes that appear in the current pathway network.")
                            ),
                            
                            conditionalPanel(
                                condition = "!output.pathway_loaded",
                                p(style = "color: #666; font-style: italic;",
                                  "Load a pathway to see which of your genes are present in it.")
                            )
                        )
                    )
                )
            )
        )
    )
)
