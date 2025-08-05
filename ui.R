# Define UI
ui <- dashboardPage(
    dashboardHeader(title = "Evo KEGG Pathways"),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem("Gene Set", tabName = "geneset", icon = icon("list")),
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
            # Gene Set Tab
            tabItem(tabName = "geneset",
                fluidRow(
                    box(
                        title = "Gene Set Management", status = "primary", solidHeader = TRUE,
                        width = 12, collapsible = TRUE,
                        
                        h4("Upload or Enter Gene Symbols"),
                        tabsetPanel(
                            tabPanel("Text Input",
                                br(),
                                div(style = "padding: 10px; background-color: #d1ecf1; border: 1px solid #bee5eb; border-radius: 5px; margin-bottom: 15px;",
                                    h5(style = "color: #0c5460; margin-top: 0;", "📋 Important: Use HGNC Gene IDs"),
                                    p(style = "color: #0c5460; margin-bottom: 5px;", 
                                      "Please use official HGNC (Human Gene Nomenclature Committee) gene symbols."),
                                    p(style = "color: #0c5460; margin-bottom: 0; font-size: 12px;",
                                      "Examples: TP53, BRCA1, EGFR, MYC (not Entrez IDs or Ensembl IDs)")
                                ),
                                textAreaInput("gene_text", "Enter Gene Symbols (ONE PER LINE):",
                                             value = "TP53\nBRCA1\nEGFR\nMYC\nGAPDH\nMTOR\nAKT1\nPIK3CA\nPTEN\nTSC1\nINS\nIGF1\nIRS1\nFOXO1\nGSK3B\nVEGFA\nHIF1A\nMAPK1\nRAF1\nRAS",
                                             placeholder = "REQUIRED: Enter gene symbols separated by newlines ONLY!\nTP53\nBRCA1\nEGFR\nMYC\nGAPDH",
                                             height = "200px"),
                                br(),
                                actionButton("load_genes_text", "Load Genes", class = "btn-info", icon = icon("upload"))
                            ),
                            tabPanel("Upload File",
                                br(),
                                div(style = "padding: 10px; background-color: #d1ecf1; border: 1px solid #bee5eb; border-radius: 5px; margin-bottom: 15px;",
                                    h5(style = "color: #0c5460; margin-top: 0;", "📋 Important: Use HGNC Gene IDs"),
                                    p(style = "color: #0c5460; margin-bottom: 5px;", 
                                      "Please use official HGNC (Human Gene Nomenclature Committee) gene symbols."),
                                    p(style = "color: #0c5460; margin-bottom: 0; font-size: 12px;",
                                      "Examples: TP53, BRCA1, EGFR, MYC (not Entrez IDs or Ensembl IDs)")
                                ),
                                fileInput("gene_file", "Upload Gene List (TXT/CSV):",
                                         accept = c(".txt", ".csv", ".tsv")),
                                helpText("REQUIRED: One gene symbol per line only (HGNC format, e.g., TP53, BRCA1)"),
                                div(style = "color: #d9534f; font-size: 11px;",
                                    "⚠️ Files must contain newline-separated gene symbols only!")
                            )
                        ),
                        
                        br(),
                        actionButton("clear_genes", "Clear Gene List", 
                                   class = "btn-warning", icon = icon("trash"))
                    )
                ),
                
                fluidRow(
                    conditionalPanel(
                        condition = "output.genes_loaded",
                        box(
                            title = "Loaded Gene List", status = "success", solidHeader = TRUE,
                            width = 8, collapsible = TRUE,
                            
                            DT::dataTableOutput("loaded_genes_table")
                        ),
                        
                        box(
                            title = "Gene Validation Results", status = "info", solidHeader = TRUE,
                            width = 4, collapsible = TRUE,
                            
                            verbatimTextOutput("gene_validation")
                        )
                    )
                )
            ),
            
            # KEGG Enrichment Tab
            tabItem(tabName = "enrichment",
                fluidRow(
                    box(
                        title = "KEGG Enrichment Analysis", status = "primary", solidHeader = TRUE,
                        width = 4, collapsible = TRUE,
                        
                        conditionalPanel(
                            condition = "output.genes_loaded",
                            h5("Analysis Parameters"),
                            numericInput("pvalue_cutoff", "P-value Cutoff:", 
                                       value = 0.05, min = 0.001, max = 0.1, step = 0.01),
                            numericInput("qvalue_cutoff", "Q-value Cutoff (FDR):", 
                                       value = 0.2, min = 0.05, max = 0.5, step = 0.05),
                            br(),
                            actionButton("run_enrichment", "Run KEGG Enrichment", 
                                       class = "btn-success", icon = icon("play")),
                            br(), br(),
                            
                            conditionalPanel(
                                condition = "output.enrichment_done",
                                h5("Visualization Options"),
                                numericInput("plot_top_n", "Show Top N Pathways:", 
                                           value = 20, min = 5, max = 50, step = 5),
                                actionButton("update_plot", "Update Plot", 
                                           class = "btn-info", icon = icon("refresh"))
                            )
                        ),
                        
                        conditionalPanel(
                            condition = "!output.genes_loaded",
                            div(style = "text-align: center; color: #666; padding: 20px;",
                                icon("dna", style = "font-size: 48px;"),
                                br(), br(),
                                h5("No Gene Set Loaded"),
                                p("Please load a gene set in the 'Gene Set' tab first."),
                                br(),
                                actionButton("goto_geneset", "Go to Gene Set Tab", 
                                           class = "btn-primary", icon = icon("arrow-right"))
                            )
                        )
                    ),
                    
                    box(
                        title = "Enrichment Results Summary", status = "info", solidHeader = TRUE,
                        width = 8, collapsible = TRUE,
                        
                        conditionalPanel(
                            condition = "output.enrichment_done",
                            verbatimTextOutput("enrichment_summary")
                        ),
                        
                        conditionalPanel(
                            condition = "!output.enrichment_done && output.genes_loaded",
                            div(style = "text-align: center; color: #666; padding: 40px;",
                                icon("chart-bar", style = "font-size: 48px;"),
                                br(), br(),
                                h5("Ready for Enrichment Analysis"),
                                p("Click 'Run KEGG Enrichment' to analyze your gene set.")
                            )
                        )
                    )
                ),
                
                fluidRow(
                    conditionalPanel(
                        condition = "output.enrichment_done",
                        box(
                            title = "Enrichment Plot", status = "success", solidHeader = TRUE,
                            width = 12, collapsible = TRUE,
                            
                            plotlyOutput("enrichment_plot", height = "600px")
                        )
                    )
                ),
                
                fluidRow(
                    conditionalPanel(
                        condition = "output.enrichment_done",
                        box(
                            title = "Detailed Enrichment Results", status = "warning", solidHeader = TRUE,
                            width = 12, collapsible = TRUE,
                            
                            DT::dataTableOutput("enrichment_table"),
                            br(),
                            fluidRow(
                                column(6,
                                    conditionalPanel(
                                        condition = "typeof input.enrichment_table_rows_selected !== 'undefined' && input.enrichment_table_rows_selected.length > 0",
                                        div(style = "background: #e7f3ff; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
                                            h5(style = "margin: 0; color: #0066cc;", "Selected Pathway:"),
                                            verbatimTextOutput("selected_enrichment_pathway", placeholder = FALSE)
                                        ),
                                        actionButton("load_enriched_pathway", "Load Selected Pathway", 
                                                   class = "btn-success", icon = icon("arrow-right"),
                                                   style = "margin-right: 10px;")
                                    )
                                ),
                                column(6,
                                    downloadButton("download_enrichment", "Download Results", 
                                                 class = "btn-primary", icon = icon("download"))
                                )
                            )
                        )
                    )
                )
            ),
            
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
                    box(
                        title = "Pathway Information", status = "info", solidHeader = TRUE,
                        width = 12, collapsible = TRUE,
                        
                        htmlOutput("pathway_info")
                    )
                )
            ),
            
            # Network Visualization Tab
            tabItem(tabName = "network",
                # Pathway header with title and KEGG link
                fluidRow(
                    conditionalPanel(
                        condition = "output.pathway_loaded",
                        box(
                            title = NULL, status = "primary", solidHeader = FALSE,
                            width = 12, 
                            style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; margin-bottom: 10px;",
                            
                            div(style = "padding: 10px;",
                                htmlOutput("pathway_header")
                            )
                        )
                    )
                ),
                
                fluidRow(
                    box(
                        title = "Pathway Network", status = "success", solidHeader = TRUE,
                        width = 8, collapsible = TRUE,
                        
                        conditionalPanel(
                            condition = "output.pathway_loaded",
                            visNetworkOutput("pathway_network", height = "700px")
                        ),
                        
                        conditionalPanel(
                            condition = "!output.pathway_loaded",
                            div(style = "text-align: center; color: #666; padding: 100px 20px;",
                                icon("project-diagram", style = "font-size: 64px; color: #ccc;"),
                                br(), br(),
                                h4("No Pathway Loaded"),
                                p("Load a pathway from the Pathway Explorer tab or select one from enrichment results."),
                                br(),
                                actionButton("goto_explorer", "Go to Pathway Explorer", 
                                           class = "btn-primary", icon = icon("search"))
                            )
                        )
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
