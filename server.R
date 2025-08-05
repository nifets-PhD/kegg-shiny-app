# Define Server
server <- function(input, output, session) {
    
    # Reactive values to store data
    values <- reactiveValues(
        pathways_list = NULL,
        selected_pathway = NULL,
        pathway_graph = NULL,
        nodes = NULL,
        edges = NULL,
        sequence_data = NULL
    )
    
    # HSA gene data for current pathway (reactive)
    pathway_hsa_data <- reactiveVal(data.frame())
    
    # Gene upload functionality
    uploaded_genes <- reactiveVal(character(0))
    
    # Load phylomap data for gene annotation
    phylomap_data <- NULL
    tryCatch({
        phylomap_file <- "data/phylomap_hgnc.tsv"
        if (file.exists(phylomap_file)) {
            phylomap_data <- read.delim(phylomap_file, stringsAsFactors = FALSE)
            cat("Loaded phylomap data with", nrow(phylomap_data), "genes\n")
        } else {
            cat("Phylomap file not found:", phylomap_file, "\n")
        }
    }, error = function(e) {
        cat("Error loading phylomap data:", e$message, "\n")
        phylomap_data <- NULL
    })
    
    # Load KEGG pathways and sequence data on app start
    observe({
        showNotification("Loading KEGG pathways...", type = "message")
        values$pathways_list <- get_kegg_pathways()
        showNotification("KEGG pathways loaded successfully!", type = "message")
    })
    
    # Handle file upload
    observeEvent(input$gene_file, {
        req(input$gene_file)
        
        tryCatch({
            genes <- load_gene_file(input$gene_file$datapath)
            uploaded_genes(genes)
            
            # Automatically validate genes
            validation <- validate_gene_symbols(genes)
            gene_validation_results(validation)
            
            showNotification(
                paste("Loaded", length(genes), "genes from file"),
                type = "message"
            )
        }, error = function(e) {
            showNotification(
                paste("Error loading file:", e$message),
                type = "error"
            )
        })
    })
    
    # Handle text input - don't auto-process, wait for button click
    observeEvent(input$gene_text, {
        # Just clear genes if text area is emptied
        if (is.null(input$gene_text) || input$gene_text == "") {
            uploaded_genes(character(0))
        }
    })
    
    # Handle Load Genes button
    observeEvent(input$load_genes_text, {
        if (!is.null(input$gene_text) && input$gene_text != "") {
            tryCatch({
                genes <- parse_gene_text(input$gene_text)
                uploaded_genes(genes)
                
                # Automatically validate genes
                validation <- validate_gene_symbols(genes)
                gene_validation_results(validation)
                
                showNotification(
                    paste("Loaded", length(genes), "genes from text"),
                    type = "message"
                )
            }, error = function(e) {
                showNotification(
                    paste("Error parsing genes:", e$message),
                    type = "error"
                )
            })
        } else {
            showNotification(
                "Please enter some gene symbols first",
                type = "warning",
                duration = 3
            )
        }
    })
    
    # Clear genes button
    observeEvent(input$clear_genes, {
        uploaded_genes(character(0))
        updateTextAreaInput(session, "gene_text", value = "")
        
        showNotification(
            "Cleared uploaded genes",
            type = "message",
            duration = 2
        )
    })
    
    # Display uploaded genes table
    output$loaded_genes_table <- renderDT({
        genes <- uploaded_genes()
        
        if (length(genes) == 0) {
            return(datatable(
                data.frame(Message = "No genes uploaded yet"),
                options = list(pageLength = 5, searching = FALSE, info = FALSE, paging = FALSE),
                rownames = FALSE
            ))
        }
        
        # Create a summary table of uploaded genes
        gene_df <- data.frame(
            Gene_Symbol = genes,
            Status = "Uploaded",
            stringsAsFactors = FALSE
        )
        
        # Check if genes are in phylomap (optional)
        if (!is.null(phylomap_data)) {
            phylomap_genes <- toupper(phylomap_data$GeneID)
            gene_df$In_Phylomap <- ifelse(toupper(gene_df$Gene_Symbol) %in% phylomap_genes, "Yes", "No")
        }
        
        datatable(
            gene_df,
            options = list(
                pageLength = 10,
                searching = TRUE,
                info = TRUE,
                lengthChange = FALSE,
                columnDefs = list(list(className = 'dt-center', targets = 1:ncol(gene_df)-1))
            ),
            rownames = FALSE
        ) %>%
        formatStyle(
            'In_Phylomap',
            backgroundColor = styleEqual('Yes', '#d4edda'),
            color = styleEqual('Yes', '#155724')
        )
    })
    
    # Gene statistics
    output$gene_stats <- renderText({
        genes <- uploaded_genes()
        if (length(genes) == 0) {
            return("No genes uploaded")
        }
        
        stats <- paste0(
            "Total genes: ", length(genes),
            " | Unique genes: ", length(unique(genes))
        )
        
        # Add phylomap overlap if available
        if (!is.null(phylomap_data)) {
            phylomap_genes <- toupper(phylomap_data$GeneID)
            overlap <- sum(toupper(genes) %in% phylomap_genes)
            stats <- paste0(stats, " | In phylomap: ", overlap, " (", round(overlap/length(genes)*100, 1), "%)")
        }
        
        return(stats)
    })
    
    # Selected genes display for Pathway Explorer tab
    output$selected_genes_display <- renderText({
        genes <- uploaded_genes()
        if (length(genes) == 0) {
            return("No genes selected yet. Upload a file or enter genes in the text box above.")
        }
        
        # Format genes in columns for better display
        genes_per_row <- 5
        gene_rows <- split(genes, ceiling(seq_along(genes) / genes_per_row))
        formatted_rows <- sapply(gene_rows, function(row) paste(row, collapse = "  |  "))
        
        paste(c(
            paste("Selected Genes (", length(genes), " total):"),
            "",
            formatted_rows
        ), collapse = "\n")
    })
    
    # Selected genes display for Network tab (same content, different output)
    output$selected_genes_display_network <- renderText({
        genes <- uploaded_genes()
        if (length(genes) == 0) {
            return("No genes selected. Go to Pathway Explorer tab to upload genes.")
        }
        
        # Format genes in columns for better display
        genes_per_row <- 6
        gene_rows <- split(genes, ceiling(seq_along(genes) / genes_per_row))
        formatted_rows <- sapply(gene_rows, function(row) paste(row, collapse = "  |  "))
        
        paste(c(
            paste("Your", length(genes), "selected genes:"),
            "",
            formatted_rows
        ), collapse = "\n")
    })
    
    # Show which uploaded genes are present in the current pathway
    output$genes_in_pathway <- renderText({
        genes <- uploaded_genes()
        nodes <- values$nodes
        
        if (length(genes) == 0) {
            return("No genes selected")
        }
        
        if (is.null(nodes) || nrow(nodes) == 0) {
            return("No pathway loaded")
        }
        
        # Find which uploaded genes are present in the pathway
        # Check against both HGNC symbols and labels
        genes_upper <- toupper(genes)
        found_genes <- character(0)
        
        # Check HGNC symbols
        if (!is.null(nodes$hgnc_symbol)) {
            hgnc_matches <- genes_upper[genes_upper %in% toupper(nodes$hgnc_symbol)]
            found_genes <- c(found_genes, hgnc_matches)
        }
        
        # Check labels (node labels)  
        if (!is.null(nodes$label)) {
            label_matches <- genes_upper[genes_upper %in% toupper(nodes$label)]
            found_genes <- c(found_genes, label_matches)
        }
        
        # Check gene names (original KEGG names)
        if (!is.null(nodes$gene_name)) {
            # Extract individual gene symbols from KEGG gene names like "hsa:2475 hsa:57521"
            all_kegg_genes <- unlist(strsplit(nodes$gene_name, "\\s+"))
            all_kegg_genes <- gsub("hsa:", "", all_kegg_genes)
            # Try to match against these KEGG IDs as well
            kegg_id_matches <- genes_upper[genes_upper %in% toupper(all_kegg_genes)]
            found_genes <- c(found_genes, kegg_id_matches)
        }
        
        # Remove duplicates and get original case
        found_genes <- unique(found_genes)
        
        if (length(found_genes) == 0) {
            return(paste0("None of your ", length(genes), " genes were found in this pathway.\n\n",
                         "Your genes: ", paste(genes, collapse = ", "), "\n\n",
                         "Try loading a different pathway that might contain your genes."))
        }
        
        # Format the found genes
        genes_per_row <- 4
        gene_rows <- split(found_genes, ceiling(seq_along(found_genes) / genes_per_row))
        formatted_rows <- sapply(gene_rows, function(row) paste(row, collapse = "  |  "))
        
        missing_genes <- genes_upper[!genes_upper %in% found_genes]
        
        result <- paste(c(
            paste("🎯 FOUND:", length(found_genes), "out of", length(genes), "genes in this pathway:"),
            "",
            formatted_rows
        ), collapse = "\n")
        
        if (length(missing_genes) > 0) {
            result <- paste0(result, "\n\n❌ Not found: ", paste(missing_genes, collapse = ", "))
        }
        
        return(result)
    })
    
    # Pathway loaded flag
    output$pathway_loaded <- reactive({
        !is.null(values$nodes) && nrow(values$nodes) > 0
    })
    outputOptions(output, "pathway_loaded", suspendWhenHidden = FALSE)
    
    # Conditional display flag for genes loaded
    output$genes_loaded <- reactive({
        length(uploaded_genes()) > 0
    })
    outputOptions(output, "genes_loaded", suspendWhenHidden = FALSE)
    
    # HSA genes table
    output$hsa_genes_table <- renderDT({
        hsa_data <- pathway_hsa_data()
        
        if (nrow(hsa_data) == 0) {
            return(datatable(
                data.frame(Message = "No HSA gene data available. Load a pathway to see gene information."),
                options = list(pageLength = 5, searching = FALSE, info = FALSE, paging = FALSE),
                rownames = FALSE
            ))
        }
        
        # Create summary table using the HSA utility function
        tryCatch({
            source("R/hsa_utils.R")
            summary_table <- create_hsa_summary_table(hsa_data)
            
            datatable(
                summary_table,
                options = list(
                    pageLength = 15,
                    searching = TRUE,
                    info = TRUE,
                    lengthChange = TRUE,
                    scrollX = TRUE,
                    columnDefs = list(
                        list(className = 'dt-center', targets = c(4, 6)),  # Chromosome and Has_Sequence
                        list(width = '100px', targets = 0),  # HSA_ID
                        list(width = '100px', targets = 1),  # Gene_Symbol
                        list(width = '200px', targets = 2),  # Gene_Name
                        list(width = '300px', targets = 3)   # Description
                    )
                ),
                rownames = FALSE,
                colnames = c("HSA ID", "Gene Symbol", "Gene Name", "Description", "Chromosome", "Gene Type", "Seq Length")
            ) %>%
            formatStyle(
                'Has_Sequence',
                backgroundColor = styleEqual(TRUE, '#d4edda'),
                color = styleEqual(TRUE, '#155724')
            ) %>%
            formatStyle(
                'Gene_Type',
                backgroundColor = styleEqual('protein_coding', '#e2f3ff'),
                color = styleEqual('protein_coding', '#0066cc')
            )
            
        }, error = function(e) {
            cat("Error creating HSA summary table:", e$message, "\n")
            return(datatable(
                data.frame(Error = paste("Error creating table:", e$message)),
                options = list(pageLength = 5, searching = FALSE, info = FALSE, paging = FALSE),
                rownames = FALSE
            ))
        })
    })
    
    # Search pathways
    observeEvent(input$search_pathways, {
        req(input$pathway_search)
        filtered_pathways <- search_pathways(values$pathways_list, input$pathway_search)
        output$pathway_table <- DT::renderDataTable({
            filtered_pathways
        }, selection = 'single', options = list(pageLength = 10))
    })
    
    # Display pathways by category
    observe({
        req(values$pathways_list, input$pathway_category)
        category_pathways <- filter_pathways_by_category(values$pathways_list, input$pathway_category)
        output$pathway_table <- DT::renderDataTable({
            category_pathways
        }, selection = 'single', options = list(pageLength = 10))
    })
    
    # Load selected pathway
    observeEvent(input$load_pathway, {
        req(input$pathway_table_rows_selected)
        
        showNotification("Loading pathway data...", type = "message")
        
        # Get selected pathway
        if (!is.null(input$pathway_search) && input$pathway_search != "") {
            filtered_pathways <- search_pathways(values$pathways_list, input$pathway_search)
            selected_row <- filtered_pathways[input$pathway_table_rows_selected, ]
        } else {
            category_pathways <- filter_pathways_by_category(values$pathways_list, input$pathway_category)
            selected_row <- category_pathways[input$pathway_table_rows_selected, ]
        }
        
        values$selected_pathway <- selected_row
        
        # Load pathway graph with HSA gene data
        tryCatch({
            pathway_data <- parse_kegg_pathway_with_hsa(selected_row$pathway_id, use_cached = TRUE, fetch_hsa_data = FALSE)
            values$pathway_graph <- pathway_data$kegg_data
            values$nodes <- pathway_data$nodes
            values$edges <- pathway_data$edges
            
            # Store HSA gene data
            pathway_hsa_data(pathway_data$hsa_data)
            
            showNotification(
                paste("Pathway loaded successfully with", nrow(pathway_data$nodes), "network nodes!"), 
                type = "message"
            )
            
            # Switch to network tab
            updateTabItems(session, "tabs", "network")
            
        }, error = function(e) {
            showNotification(paste("Error loading pathway:", e$message), type = "warning")
        })
    })
    
    # Display pathway information
    output$pathway_info <- renderUI({
        if (!is.null(values$selected_pathway)) {
            coord_info <- ""
            if (!is.null(values$nodes) && !is.null(values$nodes$x)) {
                coord_count <- sum(!is.na(values$nodes$x))
                coord_info <- paste0("KEGG Layout Coordinates: ", coord_count, " nodes positioned<br>")
            }
            
            # Create KEGG pathway URL
            pathway_id <- values$selected_pathway$pathway_id
            kegg_url <- paste0("https://www.kegg.jp/pathway/", pathway_id)
            
            HTML(paste0(
                "<div style='font-family: monospace; font-size: 14px; line-height: 1.6;'>",
                "<strong>=== PATHWAY INFORMATION ===</strong><br>",
                "<strong>Pathway ID:</strong> ", pathway_id, "<br>",
                "<strong>Name:</strong> ", values$selected_pathway$pathway_name, "<br>",
                "<strong>Description:</strong> ", values$selected_pathway$description, "<br>",
                "<strong>KEGG URL:</strong> <a href='", kegg_url, "' target='_blank' style='color: #007bff;'>", kegg_url, "</a><br><br>",
                
                "<strong>=== NETWORK STATISTICS ===</strong><br>",
                "<strong>Total Nodes:</strong> ", ifelse(!is.null(values$nodes), nrow(values$nodes), "Not loaded"), "<br>",
                "<strong>Total Edges:</strong> ", ifelse(!is.null(values$edges), nrow(values$edges), "Not loaded"), "<br>",
                coord_info, "<br>",
                
                "<strong>=== GENE ID INFORMATION ===</strong><br>",
                "• Nodes show KEGG Gene IDs prominently<br>",
                "• HGNC symbols shown in parentheses when available<br>",
                "• Click on nodes to see detailed ID mappings<br>",
                "• Use 'KEGG Original Layout' for pathway-specific positioning<br><br>",
                
                "<em style='color: #666;'>Tip: KEGG IDs are more specific for pathway analysis than gene symbols!</em>",
                "</div>"
            ))
        } else {
            HTML("<div style='font-style: italic; color: #666; padding: 20px; text-align: center;'>
                  No pathway selected - choose a pathway from the Pathway Explorer tab or from enrichment results
                  </div>")
        }
    })
    
    # Render network visualization
    output$pathway_network <- renderVisNetwork({
        req(values$nodes, values$edges)
        
        # Use user's selected coloring mode
        coloring_mode <- input$node_coloring
        highlight_genes <- uploaded_genes()
        
        # Only pass highlight_genes if the user selected the highlight_genes mode
        highlight_genes_to_use <- if(coloring_mode == "highlight_genes") highlight_genes else NULL
        
        create_kegg_network_visualization(
            values$nodes, 
            values$edges,
            show_labels = TRUE,
            show_edges = TRUE,
            coloring_mode = coloring_mode,
            highlight_genes = highlight_genes_to_use
        )
    })
    
    # Render edge legend
    output$edge_legend_table <- DT::renderDataTable({
        req(values$edges)
        
        legend_data <- create_edge_legend(values$edges)
        if (!is.null(legend_data)) {
            # Create a formatted table with colored dots
            legend_display <- data.frame(
                Type = legend_data$relationship,
                Description = legend_data$description,
                stringsAsFactors = FALSE
            )
            legend_display
        } else {
            data.frame(
                Type = "No relationships",
                Description = "No edge relationships found in this pathway",
                stringsAsFactors = FALSE
            )
        }
    }, options = list(
        pageLength = 20,  # Show more entries per page
        paging = FALSE,   # Disable paging to show all entries
        searching = FALSE,
        ordering = FALSE,
        info = FALSE,
        dom = 't',        # Only show the table
        scrollY = '350px',
        scrollCollapse = TRUE,
        columnDefs = list(
            list(width = '120px', targets = 0),  # Type column
            list(width = '250px', targets = 1)   # Description column
        )
    ), escape = FALSE)
    
    # Handle node selection
    observe({
        input$pathway_network_selected
        
        if (!is.null(input$pathway_network_selected) && !is.null(values$nodes)) {
            selected_node <- values$nodes[values$nodes$id == input$pathway_network_selected, ]
            
            output$node_info <- renderText({
                if (nrow(selected_node) > 0) {
                    # Find edges connected to this node
                    node_id <- selected_node$id
                    incoming_edges <- values$edges[values$edges$to == node_id, ]
                    outgoing_edges <- values$edges[values$edges$from == node_id, ]
                    
                    # Build relationship info
                    relationship_info <- ""
                    if (nrow(incoming_edges) > 0 || nrow(outgoing_edges) > 0) {
                        relationship_info <- paste0(
                            "\n=== RELATIONSHIPS ===", "\n",
                            "Incoming connections: ", nrow(incoming_edges), "\n",
                            "Outgoing connections: ", nrow(outgoing_edges), "\n"
                        )
                        
                        # Show some specific relationships
                        if (nrow(incoming_edges) > 0) {
                            sample_incoming <- head(incoming_edges, 3)
                            for (i in seq_len(nrow(sample_incoming))) {
                                edge <- sample_incoming[i, ]
                                rel_type <- if (!is.null(edge$relation_type) && !is.na(edge$relation_type)) edge$relation_type else "Unknown"
                                subtype <- if (!is.null(edge$subtype) && !is.na(edge$subtype) && edge$subtype != "") edge$subtype else ""
                                
                                # Find the source node name
                                source_node <- values$nodes[values$nodes$id == edge$from, ]
                                source_name <- if (nrow(source_node) > 0) source_node$hgnc_symbol else edge$from
                                
                                relationship_info <- paste0(relationship_info, 
                                    "← ", source_name, " (", rel_type, 
                                    if (subtype != "") paste0(": ", subtype) else "", ")", "\n")
                            }
                        }
                        
                        if (nrow(outgoing_edges) > 0) {
                            sample_outgoing <- head(outgoing_edges, 3)
                            for (i in seq_len(nrow(sample_outgoing))) {
                                edge <- sample_outgoing[i, ]
                                rel_type <- if (!is.null(edge$relation_type) && !is.na(edge$relation_type)) edge$relation_type else "Unknown"
                                subtype <- if (!is.null(edge$subtype) && !is.na(edge$subtype) && edge$subtype != "") edge$subtype else ""
                                
                                # Find the target node name
                                target_node <- values$nodes[values$nodes$id == edge$to, ]
                                target_name <- if (nrow(target_node) > 0) target_node$hgnc_symbol else edge$to
                                
                                relationship_info <- paste0(relationship_info, 
                                    "→ ", target_name, " (", rel_type, 
                                    if (subtype != "") paste0(": ", subtype) else "", ")", "\n")
                            }
                        }
                    }
                    
                    # Get HSA gene information
                    hsa_info_text <- ""
                    current_hsa_data <- pathway_hsa_data()
                    if (!is.null(current_hsa_data) && nrow(current_hsa_data) > 0) {
                        # Try to match by HGNC symbol or HSA ID
                        matching_gene <- NULL
                        
                        if (!is.null(selected_node$hgnc_symbol) && selected_node$hgnc_symbol != "") {
                            matching_gene <- current_hsa_data[current_hsa_data$hgnc_symbol == selected_node$hgnc_symbol, ]
                        }
                        
                        # If no match by HGNC, try to extract HSA ID from node
                        if (is.null(matching_gene) || nrow(matching_gene) == 0) {
                            # Extract HSA code from gene_name if it contains hsa:
                            if (!is.null(selected_node$gene_name) && grepl("hsa:", selected_node$gene_name)) {
                                hsa_codes <- regmatches(selected_node$gene_name, gregexpr("hsa:\\d+", selected_node$gene_name))[[1]]
                                if (length(hsa_codes) > 0) {
                                    hsa_numeric <- gsub("hsa:", "", hsa_codes[1])
                                    matching_gene <- current_hsa_data[current_hsa_data$hsa_code == hsa_numeric, ]
                                }
                            }
                        }
                        
                        if (!is.null(matching_gene) && nrow(matching_gene) > 0) {
                            gene_info <- matching_gene[1, ]
                            
                            hsa_info_text <- paste0("\n\n=== HSA GENE DATA ===\n",
                                "HSA ID: ", gene_info$hsa_code, "\n",
                                "Gene Symbol: ", gene_info$hgnc_symbol, "\n",
                                "Gene Name: ", gene_info$external_gene_name, "\n",
                                "Gene Type: ", gene_info$gene_biotype, "\n",
                                "Chromosome: ", gene_info$chromosome_name, "\n",
                                "Protein ID: ", ifelse(is.na(gene_info$protein_id), "N/A", gene_info$protein_id), "\n"
                            )
                            
                            if (!is.na(gene_info$peptide) && gene_info$peptide != "") {
                                hsa_info_text <- paste0(hsa_info_text,
                                    "Amino Acid Sequence Length: ", gene_info$peptide_length, " residues\n",
                                    "Sequence Preview: ", substr(gene_info$peptide, 1, 50), 
                                    ifelse(nchar(gene_info$peptide) > 50, "...", ""), "\n"
                                )
                            } else {
                                hsa_info_text <- paste0(hsa_info_text, "Amino Acid Sequence: Not available\n")
                            }
                            
                            if (!is.na(gene_info$description)) {
                                hsa_info_text <- paste0(hsa_info_text,
                                    "Description: ", substr(gene_info$description, 1, 200), 
                                    ifelse(nchar(gene_info$description) > 200, "...", ""), "\n"
                                )
                            }
                        }
                    }
                    
                    paste(
                        "=== GENE INFORMATION ===", "\n",
                        "KEGG Gene ID:", selected_node$kegg_id, "\n",
                        "HGNC Symbol:", selected_node$hgnc_symbol, "\n",
                        "Original KEGG ID:", selected_node$gene_name, "\n",
                        "Node Type:", selected_node$type, "\n",
                        "Display Label:", selected_node$label, "\n",
                        relationship_info, "\n",
                        "=== PATHWAY CONTEXT ===", "\n",
                        "In KEGG, this gene is referenced as:", selected_node$kegg_id, "\n",
                        "Common name/symbol:", selected_node$hgnc_symbol, "\n\n",
                        "Tip: Search for '", selected_node$kegg_id, "' in KEGG database", "\n",
                        "or '", selected_node$hgnc_symbol, "' in gene databases",
                        hsa_info_text
                    )
                } else {
                    "No node selected - click on a node to see detailed gene information"
                }
            })
        }
    })
    
    # Phylostratum legend table
    output$phylostratum_legend_table <- DT::renderDataTable({
        legend_df <- generate_phylostratum_legend()
        
        if (is.null(legend_df)) {
            return(data.frame(Message = "Phylostratum legend data not available"))
        }
        
        # Create a formatted data table with colored cells
        DT::datatable(
            legend_df[, c("Rank", "Name")],  # Show rank and name columns
            options = list(
                pageLength = 28,  # Show all strata
                paging = FALSE,   # Disable paging
                searching = FALSE, # Disable search
                info = FALSE,     # Hide table info
                dom = 't',        # Only show table
                scrollY = '350px',
                scrollCollapse = TRUE,
                columnDefs = list(
                    list(width = '40px', targets = 0, className = 'text-center'),  # Rank column
                    list(width = '160px', targets = 1)  # Name column
                )
            ),
            rownames = FALSE,
            colnames = c("PS", "Evolutionary Stage")
        ) %>%
        DT::formatStyle(
            "Rank",
            backgroundColor = DT::styleEqual(
                legend_df$Rank,
                legend_df$Color
            ),
            color = DT::styleEqual(
                legend_df$Rank,
                sapply(legend_df$Color, function(color) {
                    if (is.na(color)) return("#000000")
                    # Calculate brightness to determine text color
                    rgb_vals <- col2rgb(color)
                    brightness <- (rgb_vals[1] * 0.299 + rgb_vals[2] * 0.587 + rgb_vals[3] * 0.114) / 255
                    if (brightness < 0.5) "#FFFFFF" else "#000000"
                })
            ),
            fontWeight = "bold"
        )
    })
    
    # =========================================================================
    # ENRICHMENT ANALYSIS FUNCTIONALITY
    # =========================================================================
    
    # Reactive values for enrichment analysis
    enrichment_results <- reactiveVal(NULL)
    gene_validation_results <- reactiveVal(NULL)
    
    # Gene validation
    observeEvent(input$validate_genes, {
        genes <- uploaded_genes()
        if (length(genes) > 0) {
            showNotification("Validating genes...", type = "message")
            
            validation <- validate_gene_symbols(genes)
            gene_validation_results(validation)
            
            showNotification(
                paste("Validation complete:", length(validation$valid_genes), "valid,", 
                      length(validation$invalid_genes), "invalid"),
                type = "message"
            )
        }
    })
    
    # Gene validation output
    output$gene_validation <- renderText({
        validation <- gene_validation_results()
        if (is.null(validation)) {
            return("Click 'Validate Genes' to check gene symbols")
        }
        
        total_genes <- length(validation$valid_genes) + length(validation$invalid_genes)
        conversion_rate <- round(validation$conversion_rate * 100, 1)
        
        result <- paste0(
            "GENE VALIDATION RESULTS\n",
            "========================\n",
            "Total genes: ", total_genes, "\n",
            "Valid genes: ", length(validation$valid_genes), "\n",
            "Invalid genes: ", length(validation$invalid_genes), "\n",
            "Conversion rate: ", conversion_rate, "%\n"
        )
        
        if (length(validation$invalid_genes) > 0) {
            result <- paste0(result, "\nInvalid genes:\n", 
                           paste(validation$invalid_genes, collapse = ", "))
        }
        
        return(result)
    })
    
    # Navigation to Gene Set tab
    observeEvent(input$goto_geneset, {
        updateTabItems(session, "tabs", selected = "geneset")
    })
    
    # Navigation to Pathway Explorer tab
    observeEvent(input$goto_explorer, {
        updateTabItems(session, "tabs", selected = "explorer")
    })
    
    # Run KEGG enrichment analysis
    observeEvent(input$run_enrichment, {
        genes <- uploaded_genes()
        if (length(genes) == 0) {
            showNotification("Please load genes first", type = "warning")
            return()
        }
        
        showNotification("Running KEGG enrichment analysis...", type = "message")
        
        # Run enrichment analysis
        result <- perform_kegg_enrichment(
            genes, 
            organism = "hsa",
            pvalue_cutoff = input$pvalue_cutoff,
            qvalue_cutoff = input$qvalue_cutoff
        )
        
        enrichment_results(result)
        
        if (is.null(result)) {
            showNotification("No significant pathways found. Try adjusting p-value cutoffs.", 
                           type = "warning")
        } else {
            showNotification(paste("Found", nrow(result@result), "significant pathways!"), 
                           type = "message")
        }
    })
    
    # Enrichment completion flag
    output$enrichment_done <- reactive({
        !is.null(enrichment_results())
    })
    outputOptions(output, "enrichment_done", suspendWhenHidden = FALSE)
    
    # Enrichment summary
    output$enrichment_summary <- renderText({
        result <- enrichment_results()
        genes <- uploaded_genes()
        
        if (is.null(result) || length(genes) == 0) {
            return("No enrichment results available")
        }
        
        n_pathways <- nrow(result@result)
        n_genes <- length(genes)
        
        if (n_pathways == 0) {
            return(paste0(
                "ENRICHMENT ANALYSIS RESULTS\n",
                "============================\n",
                "Input genes: ", n_genes, "\n",
                "Significant pathways: 0\n\n",
                "No significantly enriched pathways found.\n",
                "Try increasing the p-value or q-value cutoffs."
            ))
        }
        
        top_pathway <- result@result[1, ]
        min_pval <- min(result@result$p.adjust)
        max_genes <- max(result@result$Count)
        
        summary_text <- paste0(
            "ENRICHMENT ANALYSIS RESULTS\n",
            "============================\n",
            "Input genes: ", n_genes, "\n",
            "Significant pathways: ", n_pathways, "\n",
            "Best p-value: ", format(min_pval, scientific = TRUE, digits = 3), "\n",
            "Max genes per pathway: ", max_genes, "\n\n",
            "Top pathway:\n",
            top_pathway$Description, "\n",
            "P-adjust: ", format(top_pathway$p.adjust, scientific = TRUE, digits = 3), "\n",
            "Genes: ", top_pathway$Count, "/", 
            strsplit(top_pathway$BgRatio, "/")[[1]][2]
        )
        
        return(summary_text)
    })
    
    # Enrichment results table
    output$enrichment_table <- DT::renderDataTable({
        result <- enrichment_results()
        genes <- uploaded_genes()
        
        if (is.null(result)) {
            return(data.frame(Message = "No enrichment results available"))
        }
        
        formatted_results <- format_enrichment_results(result, genes)
        
        if (nrow(formatted_results) == 0) {
            return(data.frame(Message = "No significant pathways found"))
        }
        
        DT::datatable(
            formatted_results,
            options = list(
                pageLength = 25,
                searching = TRUE,
                scrollX = TRUE,
                columnDefs = list(
                    list(width = '80px', targets = c(0)),  # Pathway ID
                    list(width = '300px', targets = c(1)), # Description
                    list(width = '80px', targets = c(2, 3, 7)), # Ratios and Count
                    list(width = '100px', targets = c(4, 5, 6, 8)), # P-values and Fold
                    list(width = '200px', targets = c(9)) # Genes
                )
            ),
            selection = 'single',  # Enable single row selection
            rownames = FALSE
        ) %>%
        DT::formatSignif(columns = c('pvalue', 'p.adjust', 'qvalue'), digits = 3) %>%
        DT::formatRound(columns = 'FoldEnrichment', digits = 2) %>%
        DT::formatStyle(
            'p.adjust',
            backgroundColor = DT::styleInterval(
                cuts = c(0.001, 0.01, 0.05),
                values = c('#d4edda', '#fff3cd', '#f8d7da', '#ffffff')
            )
        ) %>%
        DT::formatStyle(
            'FoldEnrichment',
            backgroundColor = DT::styleInterval(
                cuts = c(1.5, 2, 5),
                values = c('#ffffff', '#e7f3ff', '#cce7ff', '#99d6ff')
            )
        )
    })
    
    # Enrichment plot
    output$enrichment_plot <- renderPlotly({
        result <- enrichment_results()
        
        if (is.null(result)) {
            return(NULL)
        }
        
        plot_data <- create_enrichment_plot_data(result, input$plot_top_n)
        
        if (is.null(plot_data) || nrow(plot_data) == 0) {
            return(NULL)
        }
        
        # Create the plot
        p <- ggplot(plot_data, aes(x = FoldEnrichment, y = PathwayShort)) +
            geom_point(aes(size = Count, color = p.adjust), alpha = 0.7) +
            scale_color_gradient(low = "red", high = "blue", name = "Adj. P-value") +
            scale_size_continuous(name = "Gene Count", range = c(3, 10)) +
            labs(
                title = paste("Top", nrow(plot_data), "Enriched KEGG Pathways"),
                x = "Fold Enrichment",
                y = "Pathway"
            ) +
            theme_minimal() +
            theme(
                axis.text.y = element_text(size = 10),
                plot.title = element_text(size = 14, hjust = 0.5),
                legend.position = "right"
            )
        
        ggplotly(p, tooltip = c("x", "y", "size", "colour"), height = 600) %>%
            layout(margin = list(l = 300))
    })
    
    # Update plot when parameters change
    observeEvent(input$update_plot, {
        # The plot will automatically update due to reactive dependencies
        showNotification("Plot updated", type = "message", duration = 2)
    })
    
    # Download enrichment results
    output$download_enrichment <- downloadHandler(
        filename = function() {
            paste0("kegg_enrichment_results_", Sys.Date(), ".csv")
        },
        content = function(file) {
            result <- enrichment_results()
            genes <- uploaded_genes()
            
            if (!is.null(result)) {
                formatted_results <- format_enrichment_results(result, genes)
                write.csv(formatted_results, file, row.names = FALSE)
            }
        }
    )
    
    # Show selected pathway info from enrichment table
    output$selected_enrichment_pathway <- renderText({
        if (is.null(input$enrichment_table_rows_selected)) {
            return("")
        }
        
        result <- enrichment_results()
        if (!is.null(result)) {
            formatted_results <- format_enrichment_results(result, uploaded_genes())
            selected_row <- formatted_results[input$enrichment_table_rows_selected, ]
            
            paste0(
                "Pathway: ", selected_row$Pathway, "\n",
                "Description: ", selected_row$Description, "\n",
                "P-adjust: ", format(selected_row$p.adjust, scientific = TRUE, digits = 3), "\n",
                "Genes: ", selected_row$Count, " genes"
            )
        }
    })
    
    # Load pathway selected from enrichment results with button
    observeEvent(input$load_enriched_pathway, {
        req(input$enrichment_table_rows_selected)
        
        result <- enrichment_results()
        if (!is.null(result)) {
            formatted_results <- format_enrichment_results(result, uploaded_genes())
            selected_row <- formatted_results[input$enrichment_table_rows_selected, ]
            
            # Extract pathway ID (should be like "hsa04110")
            pathway_id <- selected_row$Pathway
            
            # Show progress notification
            progress <- Progress$new()
            progress$set(message = paste("Loading pathway", pathway_id), value = 0.1)
            on.exit(progress$close())
            
            showNotification(paste("Loading enriched pathway:", pathway_id), type = "message")
            
            # Create a mock pathway selection row similar to pathway_table format
            selected_pathway <- data.frame(
                pathway_id = pathway_id,
                pathway_name = selected_row$Description,
                description = selected_row$Description,
                stringsAsFactors = FALSE
            )
            
            values$selected_pathway <- selected_pathway
            progress$set(value = 0.3)
            
            # Load pathway graph with HSA gene data
            tryCatch({
                progress$set(message = "Parsing pathway data...", value = 0.5)
                pathway_data <- parse_kegg_pathway_with_hsa(pathway_id, use_cached = TRUE, fetch_hsa_data = FALSE)
                
                progress$set(message = "Building network...", value = 0.7)
                values$pathway_graph <- pathway_data$kegg_data
                values$nodes <- pathway_data$nodes
                values$edges <- pathway_data$edges
                
                # Store HSA gene data
                pathway_hsa_data(pathway_data$hsa_data)
                
                progress$set(message = "Finalizing...", value = 0.9)
                
                showNotification(
                    paste("Enriched pathway loaded successfully with", nrow(pathway_data$nodes), "network nodes!"), 
                    type = "message"
                )
                
                progress$set(value = 1.0)
                
                # Switch to network tab
                updateTabItems(session, "tabs", "network")
                
            }, error = function(e) {
                showNotification(paste("Error loading enriched pathway:", e$message), type = "warning")
            })
        }
    })
    
    # Remove the old row selection observer since we now use a button
    # observeEvent(input$enrichment_table_rows_selected, { ... }) - REMOVED
    
    # Pathway header for network visualization
    output$pathway_header <- renderUI({
        if (!is.null(values$selected_pathway)) {
            pathway_id <- values$selected_pathway$pathway_id
            kegg_url <- paste0("https://www.kegg.jp/pathway/", pathway_id)
            
            HTML(paste0(
                "<div style='text-align: center;'>",
                "<h3 style='margin: 5px 0; color: white; text-shadow: 1px 1px 2px rgba(0,0,0,0.5);'>",
                values$selected_pathway$pathway_name,
                "</h3>",
                "<p style='margin: 5px 0; font-size: 16px; color: #e8f4f8;'>",
                "<strong>Pathway ID:</strong> ", pathway_id, " | ",
                "<a href='", kegg_url, "' target='_blank' style='color: #fff; text-decoration: underline;'>",
                "View on KEGG Database ↗",
                "</a>",
                "</p>",
                "</div>"
            ))
        }
    })
}
