# Define UI
ui <- dashboardPage(
    dashboardHeader(title = "Evo KEGG Pathways"),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem("Gene Set", tabName = "geneset", icon = icon("list")),
            menuItem("KEGG Enrichment", tabName = "enrichment", icon = icon("chart-bar")),
            menuItem("Pathway Explorer", tabName = "explorer", icon = icon("search")),
            menuItem("Network Visualization", tabName = "network", icon = icon("project-diagram")),
            menuItem("Evolutionary Transcriptomics", tabName = "evolution", icon = icon("dna"))
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
                        
                        htmlOutput("node_info")
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
            ),
            
            # Evolutionary Transcriptomics Tab
            tabItem(tabName = "evolution",
                fluidRow(
                    box(
                        title = "Expression Data Upload", status = "primary", solidHeader = TRUE,
                        width = 6, collapsible = TRUE,
                        
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
                            )
                        )
                    ),
                    
                    box(
                        title = "Phylostratum Mapping", status = "info", solidHeader = TRUE,
                        width = 6, collapsible = TRUE,
                        
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
                            
                            br(), br(),
                            conditionalPanel(
                                condition = "output.phyloset_created",
                                div(style = "background: #d1ecf1; padding: 10px; border-radius: 5px; border-left: 4px solid #17a2b8;",
                                    h5(style = "margin: 0; color: #0c5460;", "✓ PhyloExpressionSet Created"),
                                    verbatimTextOutput("phyloset_summary", placeholder = FALSE)
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
                
                # Visualization Options and Plots
                fluidRow(
                    conditionalPanel(
                        condition = "output.phyloset_created",
                        box(
                            title = "Evolutionary Transcriptomics Visualizations", 
                            status = "success", solidHeader = TRUE,
                            width = 12, collapsible = TRUE,
                            
                            h4("Select Visualization Type"),
                            p("Explore different aspects of transcriptomic age in your dataset:"),
                            
                            fluidRow(
                                column(4,
                                    h5("Primary Analyses"),
                                    actionButton("plot_tai_signature", "TAI Signature", 
                                               class = "btn-primary btn-block", 
                                               icon = icon("line-chart"),
                                               style = "margin-bottom: 5px;"),
                                    helpText("Transcriptome Age Index over developmental stages/conditions"),
                                    
                                    actionButton("plot_distribution_strata", "Distribution of Phylostrata", 
                                               class = "btn-primary btn-block",
                                               icon = icon("bar-chart"),
                                               style = "margin-bottom: 5px;"),
                                    helpText("Distribution of genes across evolutionary ages")
                                ),
                                column(4,
                                    h5("Expression Patterns"),
                                    actionButton("plot_gene_heatmap", "Gene Expression Heatmap", 
                                               class = "btn-warning btn-block",
                                               icon = icon("th"),
                                               style = "margin-bottom: 5px;"),
                                    helpText("Heatmap of gene expression by phylostratum"),
                                    
                                    actionButton("plot_contribution", "Phylostratum Contribution", 
                                               class = "btn-warning btn-block",
                                               icon = icon("pie-chart"),
                                               style = "margin-bottom: 5px;"),
                                    helpText("Contribution of each phylostratum to overall TAI")
                                ),
                                column(4,
                                    h5("Advanced Analyses"),
                                    actionButton("plot_gene_space", "Gene Space", 
                                               class = "btn-info btn-block",
                                               icon = icon("chart-scatter"),
                                               style = "margin-bottom: 5px;"),
                                    helpText("PCA of genes colored by phylostratum"),
                                    
                                    actionButton("plot_sample_space", "Sample Space", 
                                               class = "btn-info btn-block",
                                               icon = icon("object-group"),
                                               style = "margin-bottom: 5px;"),
                                    helpText("PCA of samples with phylostratum information"),
                                    
                                    actionButton("plot_gene_profiles", "Gene Expression Profiles", 
                                               class = "btn-info btn-block",
                                               icon = icon("line"),
                                               style = "margin-bottom: 5px;"),
                                    helpText("Individual gene expression profiles by phylostratum")
                                )
                            )
                        )
                    )
                ),
                
                # Plot Display Area
                fluidRow(
                    conditionalPanel(
                        condition = "output.evolution_plot_ready",
                        box(
                            title = NULL, status = "success", solidHeader = TRUE,
                            width = 12, collapsible = TRUE,
                            
                            # Dynamic title
                            htmlOutput("evolution_plot_title"),
                            
                            # Plot output
                            plotOutput("evolution_plot", height = "600px"),
                            
                            br(),
                            fluidRow(
                                column(6,
                                    downloadButton("download_evolution_plot", "Download Plot", 
                                                 class = "btn-primary", icon = icon("download"))
                                ),
                                column(6,
                                    div(style = "text-align: right;",
                                        actionButton("clear_evolution_plot", "Clear Plot", 
                                                   class = "btn-secondary", icon = icon("eraser"))
                                    )
                                )
                            )
                        )
                    )
                ),
                
                # Data Summary and Legend
                fluidRow(
                    conditionalPanel(
                        condition = "output.phyloset_created",
                        box(
                            title = "Phylostratum Legend", status = "warning", solidHeader = TRUE,
                            width = 6, collapsible = TRUE, collapsed = TRUE,
                            
                            h5("Evolutionary Timeline"),
                            p("Each phylostratum represents a major evolutionary transition:"),
                            div(style = "max-height: 400px; overflow-y: auto;",
                                DT::dataTableOutput("evolution_strata_legend")
                            )
                        ),
                        
                        box(
                            title = "Dataset Information", status = "info", solidHeader = TRUE,
                            width = 6, collapsible = TRUE, collapsed = TRUE,
                            
                            h5("Current Dataset"),
                            verbatimTextOutput("evolution_dataset_info"),
                            
                            br(),
                            h5("Analysis Notes"),
                            p("• TAI (Transcriptome Age Index) measures the relative evolutionary age of gene expression"),
                            p("• Lower TAI = predominantly ancient genes expressed"),  
                            p("• Higher TAI = more evolutionarily recent genes expressed"),
                            p("• Phylostrata range from 1 (most ancient) to 28 (Homo sapiens specific)")
                        )
                    )
                )
            )
        )
    )
)
