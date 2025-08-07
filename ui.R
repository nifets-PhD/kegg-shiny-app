# Load color configuration
source("R/color_config.R")

# Define UI
ui <- dashboardPage(
    skin = "purple",  # Use purple skin
    dashboardHeader(title = "Pathway Evolution Analysis"),
    
    dashboardSidebar(
        sidebarMenu(id = "tabs",
            menuItem("Integrated Visualisation", tabName = "network", icon = icon("project-diagram")),
            menuItem("Pathway Explorer", tabName = "explorer", icon = icon("search")),
            menuItem("Your Gene Set", tabName = "geneset", icon = icon("list")),
            menuItem("KEGG Enrichment", tabName = "enrichment", icon = icon("chart-bar")),
            menuItem("Your Expression Data", tabName = "evolution", icon = icon("dna"))
        ),
        
        # Add footer to sidebar
        div(
            style = "position: fixed; bottom: 10px; left: 0; right: 230px; text-align: center; 
                     color: #999; font-size: 11px; padding: 0 10px; z-index: 1000;
                     background-color: rgba(28, 23, 47, 0.9); margin-left: 0; width: 230px;",
            hr(style = "border-color: #444; margin: 5px 0;"),
            div(
                style = "display: flex; align-items: center; justify-content: center; gap: 5px;",
                span("Powered by"),
                tags$a(
                    href = "https://github.com/drostlab/myTAI",
                    target = "_blank",
                    style = "color: #bb86fc; font-weight: bold; text-decoration: none;",
                    "myTAI"
                )
            )
        )
    ),
    
    dashboardBody(
        tags$head(
            tags$style(HTML(paste0("
                .content-wrapper, .right-side {
                    background-color: #f4f4f4;
                }
                .box {
                    margin-bottom: 20px;
                }
                ",
                generate_custom_css()  # Add custom purple/green theme
            ))),
            
            # JavaScript for EMERALD alignment functionality
            tags$script(HTML("
                Shiny.addCustomMessageHandler('updateEmeraldLink', function(message) {
                    var container = document.getElementById('emerald_link_container');
                    if (container) {
                        container.innerHTML = message.html;
                    }
                });
                
                Shiny.addCustomMessageHandler('updateEmeraldContainer', function(data) {
                    var container = document.getElementById(data.containerId);
                    if (container) {
                        container.innerHTML = data.html;
                    }
                });
            "))
        ),
        
        tabItems(
            # Gene Set Tab
            tabItem(tabName = "geneset",
                fluidRow(
                    box(
                        title = "Gene Set Management", status = "primary", solidHeader = TRUE,
                        width = 12, collapsible = TRUE,
                        
                        h4("Upload or Enter Gene IDs"),
                        
                        # Gene ID Type Selection
                        fluidRow(
                            column(6,
                                selectInput("gene_id_type", "Select Gene ID Type:",
                                           choices = list(
                                               "Entrez Gene ID (Recommended for KEGG)" = "entrez",
                                               "Gene Symbol (HGNC)" = "symbol", 
                                               "Ensembl Gene ID" = "ensembl",
                                               "UniProt ID" = "uniprot"
                                           ),
                                           selected = "entrez")
                            ),
                            column(6,
                                div(style = "margin-top: 25px; font-size: 12px; color: #666;",
                                    "💡 Entrez IDs provide best KEGG pathway mapping")
                            )
                        ),
                        
                        tabsetPanel(
                            tabPanel("Text Input",
                                br(),
                                # Dynamic instruction panel based on selected ID type
                                uiOutput("gene_id_instructions"),
                                # Dynamic text area that updates examples based on gene ID type
                                uiOutput("gene_text_input"),
                                br(),
                                actionButton("load_genes_text", "Load Genes", class = "btn-info", icon = icon("upload"))
                            ),
                            tabPanel("Upload File",
                                br(),
                                # Dynamic instruction panel for file upload
                                uiOutput("file_upload_instructions"),
                                fileInput("gene_file", "Upload Gene List (TXT/CSV):",
                                         accept = c(".txt", ".csv", ".tsv")),
                                helpText("One gene ID per line. File format will be validated based on selected ID type."),
                                div(style = "color: #d9534f; font-size: 11px;",
                                    "⚠️ Files must contain newline-separated gene IDs only!")
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
                            width = 12, collapsible = TRUE,
                            
                            DT::dataTableOutput("loaded_genes_table")
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
                            title = "Detailed Enrichment Results", status = "info", solidHeader = TRUE,
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
            
            # Pathway Visualisation Tab
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
                    # Column 1: Pathway Network and Expression Profiles
                    column(8,
                        # Pathway Network
                        box(
                            title = "Pathway Network", status = "success", solidHeader = TRUE,
                            width = 12, collapsible = TRUE,
                            
                            conditionalPanel(
                                condition = "output.pathway_loaded",
                                visNetworkOutput("pathway_network", height = "450px")
                            ),
                            
                            conditionalPanel(
                                condition = "!output.pathway_loaded",
                                div(style = "text-align: center; color: #666; padding: 60px 20px;",
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
                        
                        # Expression Profiles - always show
                        box(
                            title = "Gene Expression Profiles", status = "info", solidHeader = TRUE,
                            width = 12, collapsible = TRUE,
                            
                            # Show plot when both phyloset and pathway are loaded
                            conditionalPanel(
                                condition = "output.phyloset_created && output.pathway_loaded",
                                plotOutput("pathway_expression_plot", height = "600px")
                            ),
                            
                            # Direct to evolution tab when no phyloset created
                            conditionalPanel(
                                condition = "!output.phyloset_created",
                                div(style = "text-align: center; color: #666; padding: 40px 20px;",
                                    icon("dna", style = "font-size: 48px; color: #ccc;"),
                                    br(), br(),
                                    h4("No Expression Dataset Loaded"),
                                    p("Load an expression dataset to view gene expression profiles colored by evolutionary age."),
                                    br(),
                                    actionButton("goto_evolution_from_network", "Go to Your Expression Data", 
                                               class = "btn-success", icon = icon("dna"))
                                )
                            ),
                            
                            # Show when phyloset created but no pathway loaded
                            conditionalPanel(
                                condition = "output.phyloset_created && !output.pathway_loaded",
                                div(style = "text-align: center; color: #666; padding: 40px 20px;",
                                    icon("project-diagram", style = "font-size: 48px; color: #ccc;"),
                                    br(), br(),
                                    h4("No Pathway Selected"),
                                    p("Load a pathway to view expression profiles of pathway genes."),
                                    br(),
                                    actionButton("goto_explorer_from_expression", "Go to Pathway Explorer", 
                                               class = "btn-success", icon = icon("search"))
                                )
                            )
                        )
                    ),
                    
                    # Column 2: Network Options, Selected Node Info, Edge Legend
                    column(4,
                        # Options
                        box(
                            title = "Options", status = "primary", solidHeader = TRUE,
                            width = 12, collapsible = TRUE,
                            
                            h5("Gene Colouring"),
                            radioButtons("node_coloring", "Colour genes by:",
                                       choices = list(
                                           "KEGG Default" = "kegg_default",
                                           "Phylostratum" = "phylostratum",
                                           "Highlight Uploaded Genes" = "highlight_genes"
                                       ),
                                       selected = "kegg_default"),
                            
                            # Show selected genes when highlighting uploaded genes
                            conditionalPanel(
                                condition = "input.node_coloring == 'highlight_genes' && output.genes_loaded",
                                hr(),
                                h6("Highlighted Genes:"),
                                div(style = "background: #f8f9fa; padding: 4px; border-radius: 4px; max-height: 150px; overflow-y: auto;",
                                    DT::dataTableOutput("network_highlighted_genes", height = "120px")
                                )
                            ),
                            
                            # Direct to gene set tab when highlight_genes selected but no genes loaded
                            conditionalPanel(
                                condition = "input.node_coloring == 'highlight_genes' && !output.genes_loaded",
                                hr(),
                                div(style = "background: #fff3cd; padding: 12px; border-radius: 4px; border-left: 4px solid #856404;",
                                    h6(style = "margin: 0 0 8px 0; color: #856404;", "No Genes to Highlight"),
                                    p(style = "margin: 0 0 10px 0; font-size: 13px; color: #856404;", 
                                      "Upload your gene set first to highlight genes in the pathway."),
                                    actionButton("goto_geneset_from_network", "Go to Your Gene Set Tab", 
                                               class = "btn-warning btn-sm", icon = icon("upload"))
                                )
                            ),
                            
                            conditionalPanel(
                                condition = "input.node_coloring == 'phylostratum'",
                                hr(),
                                h6("Phylostratum Legend"),
                                helpText("Colors genes by evolutionary age (phylostratum)."),
                                div(style = "background: #f8f9fa; padding: 4px; border-radius: 4px; max-height: 180px;",
                                    DT::dataTableOutput("phylostratum_legend_table", height = "160px")
                                )
                            )
                        ),
                        
                        # Selected Node Information
                        box(
                            title = "Selected Node Information", status = "info", solidHeader = TRUE,
                            width = 12, collapsible = TRUE,
                            
                            htmlOutput("node_info")
                        )
                    )
                ),
                
                # Edge Legend at bottom
                fluidRow(
                    box(
                        title = "Edge Legend", status = "info", solidHeader = TRUE,
                        width = 12, collapsible = TRUE, collapsed = TRUE,
                        
                        h5("Relationship Types"),
                        p("Different colours and styles represent different biological relationships:"),
                        div(style = "max-height: 300px; overflow-y: auto;",
                            DT::dataTableOutput("edge_legend_table")
                        )
                    )
                )
            ),
            
            # Evolutionary Transcriptomics Tab
            tabItem(tabName = "evolution",
                fluidRow(
                    box(
                        title = "Expression Data Upload", status = "primary", solidHeader = TRUE,
                        width = 6, collapsible = TRUE, height = "auto", style = "min-height: 650px;",
                        
                        h4("Upload Gene Expression Dataset"),
                        
                        # Gene ID Type Selection for expression data
                        selectInput("expr_gene_id_type", "Gene ID Type in Expression Data:",
                                   choices = list(
                                       "Gene Symbol (HGNC)" = "symbol",
                                       "Entrez Gene ID" = "entrez",
                                       "Ensembl Gene ID" = "ensembl",
                                       "UniProt ID" = "uniprot"
                                   ),
                                   selected = "symbol"),
                        
                        fileInput("expression_file", "Choose Expression File (CSV/TSV):",
                                 accept = c(".csv", ".tsv", ".txt")),
                        
                        helpText("Upload a gene expression dataset with genes as rows and samples as columns."),
                        helpText("First column should contain gene IDs, remaining columns should be expression values."),
                        
                        br(),
                        h5("Example Data"),
                        actionButton("load_example_expr", "Load Example Dataset", 
                                   class = "btn-info", icon = icon("download")),
                        helpText("Loads brain cell type expression data for demonstration."),
                        
                        br(),
                        conditionalPanel(
                            condition = "output.expression_loaded",
                            div(style = "background: #d4edda; padding: 10px; border-radius: 5px; border-left: 4px solid #28a745;",
                                h5(style = "margin: 0; color: #155724;", "✓ Expression Data Loaded"),
                                verbatimTextOutput("expression_summary", placeholder = FALSE)
                            ),
                            
                            br(),
                            actionButton("toggle_expression_preview", 
                                       "Show Expression Matrix Preview (Top 50 genes)", 
                                       class = "btn-info btn-sm",
                                       icon = icon("eye")),
                            br(), br(),
                            conditionalPanel(
                                condition = "input.toggle_expression_preview % 2 == 1",
                                div(style = "margin-bottom: 15px;",
                                    DT::dataTableOutput("expression_preview", height = "350px")
                                )
                            )
                        )
                    ),
                    
                    box(
                        title = "Phylostratum Mapping", status = "success", solidHeader = TRUE,
                        width = 6, collapsible = TRUE, height = "auto", style = "min-height: 650px;",
                        
                        conditionalPanel(
                            condition = "output.expression_loaded",
                            h4("Create Phylostratum Mapping"),
                            p("Map your expression data to evolutionary phylostrata for transcriptomic age analysis."),
                            
                            # Group assignment for samples (optional)
                            h5("Sample Groups (Optional)"),
                            helpText("Assign sample groups for comparative analysis. Leave empty to use sample names."),
                            textAreaInput("sample_groups", "Sample Groups (comma-separated):",
                                        placeholder = "e.g., Control,Control,Treatment,Treatment...",
                                        height = "80px"),
                            
                            br(),
                            actionButton("create_phyloset", "Create PhyloExpressionSet", 
                                       class = "btn-success", icon = icon("play")),
                                       
                            # Show mapping results after phyloset creation
                            conditionalPanel(
                                condition = "output.phyloset_created",
                                br(),
                                div(style = "background: #d1ecf1; padding: 15px; border-radius: 5px; border-left: 4px solid #17a2b8;",
                                    h5(style = "margin-top: 0; color: #0c5460;", "✓ PhyloExpressionSet Created"),
                                    verbatimTextOutput("phyloset_mapping_info", placeholder = FALSE)
                                )
                            )
                        ),
                        
                        conditionalPanel(
                            condition = "!output.expression_loaded",
                            div(style = "text-align: center; color: #666; padding: 20px;",
                                icon("upload", style = "font-size: 48px;"),
                                br(), br(),
                                h5("No Expression Data"),
                                p("Please upload expression data or load example dataset first.")
                            )
                        )
                    )
                ),
                
                # Evolutionary Transcriptomics Visualizations
                fluidRow(
                    conditionalPanel(
                        condition = "output.phyloset_created",
                        box(
                            title = "Evolutionary Transcriptomics Analysis", status = "primary", solidHeader = TRUE,
                            width = 12, collapsible = TRUE,
                            
                            fluidRow(
                                column(6,
                                    h5("TAI Signature", style = "text-align: center; margin-bottom: 15px;"),
                                    plotOutput("tai_signature_plot", height = "350px")
                                ),
                                column(6,
                                    h5("Distribution of Phylostrata", style = "text-align: center; margin-bottom: 15px;"),
                                    plotOutput("phylostrata_distribution_plot", height = "350px")
                                )
                            ),
                            
                            br(),
                            
                            fluidRow(
                                column(6,
                                    h5("Sample Space Analysis", style = "text-align: center; margin-bottom: 15px;"),
                                    plotOutput("sample_space_plot", height = "350px")
                                ),
                                column(6,
                                    h5("Gene Expression Heatmap", style = "text-align: center; margin-bottom: 15px;"),
                                    plotOutput("gene_heatmap_plot", height = "350px")
                                )
                            )
                        )
                    )
                )
            )
        )
    )
)
