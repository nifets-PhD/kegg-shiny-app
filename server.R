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
    
    # Gene upload functionality
    # Reactive values for uploaded genes
    uploaded_genes <- reactiveVal(character(0))
    converted_entrez_genes <- reactiveVal(character(0))  # Store Entrez IDs for KEGG operations
    gene_id_type <- reactiveVal("entrez")  # Track the current gene ID type
    
    # Reactive values for gene validation
    gene_validation_results <- reactiveVal(NULL)
    
    # Dynamic instruction panels based on gene ID type
    output$gene_id_instructions <- renderUI({
        id_type <- input$gene_id_type
        
        instructions <- switch(id_type,
            "entrez" = list(
                title = "📋 Entrez Gene IDs (Recommended)",
                desc = "Use numeric Entrez Gene IDs for best KEGG pathway mapping. No conversion needed.",
                examples = "Examples: 1956 (EGFR), 7157 (TP53), 672 (BRCA1), 4233 (MET)"
            ),
            "symbol" = list(
                title = "📋 Gene Symbols (HGNC)",
                desc = "Official HGNC gene symbols. Will be automatically converted to Entrez IDs for KEGG analysis.",
                examples = "Examples: TP53, BRCA1, EGFR, MYC"
            ),
            "ensembl" = list(
                title = "📋 Ensembl Gene IDs",
                desc = "Ensembl gene identifiers. Will be automatically converted to Entrez IDs for KEGG analysis.",
                examples = "Examples: ENSG00000141510, ENSG00000012048, ENSG00000146648"
            ),
            "uniprot" = list(
                title = "📋 UniProt IDs",
                desc = "UniProt protein identifiers. Will be automatically converted to Entrez IDs for KEGG analysis.",
                examples = "Examples: P04637, P38398, P00533, P01106"
            )
        )
        
        div(style = "padding: 10px; background-color: #d1ecf1; border: 1px solid #bee5eb; border-radius: 5px; margin-bottom: 15px;",
            h5(style = "color: #0c5460; margin-top: 0;", instructions$title),
            p(style = "color: #0c5460; margin-bottom: 5px;", instructions$desc),
            p(style = "color: #0c5460; margin-bottom: 0; font-size: 12px;", instructions$examples)
        )
    })
    
    output$file_upload_instructions <- renderUI({
        id_type <- input$gene_id_type
        
        instructions <- switch(id_type,
            "entrez" = list(
                title = "📋 Entrez Gene IDs (Recommended)",
                desc = "Upload files containing numeric Entrez Gene IDs for best KEGG pathway mapping.",
                examples = "Examples: 1956, 7157, 672, 4233"
            ),
            "symbol" = list(
                title = "📋 Gene Symbols (HGNC)",
                desc = "Upload files containing official HGNC gene symbols.",
                examples = "Examples: TP53, BRCA1, EGFR, MYC"
            ),
            "ensembl" = list(
                title = "📋 Ensembl Gene IDs",
                desc = "Upload files containing Ensembl gene identifiers.",
                examples = "Examples: ENSG00000141510, ENSG00000012048"
            ),
            "uniprot" = list(
                title = "📋 UniProt IDs",
                desc = "Upload files containing UniProt protein identifiers.",
                examples = "Examples: P04637, P38398, P00533"
            )
        )
        
        div(style = "padding: 10px; background-color: #d1ecf1; border: 1px solid #bee5eb; border-radius: 5px; margin-bottom: 15px;",
            h5(style = "color: #0c5460; margin-top: 0;", instructions$title),
            p(style = "color: #0c5460; margin-bottom: 5px;", instructions$desc),
            p(style = "color: #0c5460; margin-bottom: 0; font-size: 12px;", instructions$examples)
        )
    })
    
    # Dynamic text area with examples that change based on gene ID type
    output$gene_text_input <- renderUI({
        id_type <- input$gene_id_type
        
        # Define examples for each ID type
        examples <- switch(id_type,
            "entrez" = list(
                value = "1956\n672\n1950\n4233\n2597\n2475\n207\n5290\n5728\n7249\n3630\n3479\n3667\n2308\n2932\n7422\n3091\n5594\n5894",
                placeholder = "Enter Entrez Gene IDs separated by newlines!\n1956\n672\n1950\n4233\n2597"
            ),
            "symbol" = list(
                value = "TP53\nBRCA1\nEGFR\nMYC\nGAPDH\nMTOR\nAKT1\nPIK3CA\nPTEN\nTSC1\nINS\nIGF1\nIRS1\nFOXO1\nGSK3B\nVEGFA\nHIF1A\nMAPK1\nRAF1\nRAS",
                placeholder = "Enter gene symbols separated by newlines!\nTP53\nBRCA1\nEGFR\nMYC\nGAPDH"
            ),
            "ensembl" = list(
                value = "ENSG00000141510\nENSG00000012048\nENSG00000146648\nENSG00000136997\nENSG00000111640\nENSG00000198793\nENSG00000142208\nENSG00000171862\nENSG00000171791\nENSG00000165699\nENSG00000254647\nENSG00000140443\nENSG00000169032\nENSG00000149311\nENSG00000168036\nENSG00000112715\nENSG00000100644\nENSG00000100030\nENSG00000132155\nENSG00000174775",
                placeholder = "Enter Ensembl Gene IDs separated by newlines!\nENSG00000141510\nENSG00000012048\nENSG00000146648"
            ),
            "uniprot" = list(
                value = "P04637\nP38398\nP00533\nP01106\nP04406\nP42345\nP31749\nP42336\nP60484\nQ92574\nP01308\nP05019\nP35568\nO14757\nP49841\nP15692\nQ16665\nP28482\nP04049\nP01112",
                placeholder = "Enter UniProt IDs separated by newlines!\nP04637\nP38398\nP00533\nP01106"
            )
        )
        
        textAreaInput("gene_text", "Enter Gene IDs (ONE PER LINE):",
                     value = examples$value,
                     placeholder = examples$placeholder,
                     height = "200px")
    })
    
    # Load comprehensive phylomap data for gene annotation
    phylomap_data <- NULL
    tryCatch({
        phylomap_file <- "data/phylomap.tsv"  # Use comprehensive phylomap
        if (file.exists(phylomap_file)) {
            phylomap_data <- read.delim(phylomap_file, stringsAsFactors = FALSE)
            cat("Loaded comprehensive phylomap data with", nrow(phylomap_data), "protein entries\n")
        } else {
            cat("Comprehensive phylomap file not found:", phylomap_file, "\n")
        }
    }, error = function(e) {
        cat("Error loading comprehensive phylomap data:", e$message, "\n")
        phylomap_data <- NULL
    })
    
    # Load comprehensive gene ID mapping
    comprehensive_mapping <- NULL
    tryCatch({
        comprehensive_mapping <- download_hgnc_uniprot_mapping()
        if (nrow(comprehensive_mapping) > 0) {
            cat("Loaded comprehensive gene mapping with", nrow(comprehensive_mapping), "entries\n")
        }
    }, error = function(e) {
        cat("Error loading comprehensive gene mapping:", e$message, "\n")
        comprehensive_mapping <- NULL
    })
    
    # Load KEGG pathways and sequence data on app start
    observe({
        showNotification("Loading KEGG pathways...", type = "message")
        values$pathways_list <- get_kegg_pathways()
        showNotification("KEGG pathways loaded successfully!", type = "message")
    })
    
    # Handle file upload
    observeEvent(input$gene_file, {
        req(input$gene_file, input$gene_id_type)
        
        tryCatch({
            genes <- load_gene_file(input$gene_file$datapath)
            uploaded_genes(genes)
            gene_id_type(input$gene_id_type)
            
            # Automatically validate genes with the selected ID type
            validation <- validate_gene_ids(genes, input$gene_id_type)
            gene_validation_results(validation)
            
            # Store converted Entrez IDs for KEGG pathway operations
            if (!is.null(validation$entrez_mapping) && nrow(validation$entrez_mapping) > 0) {
                entrez_ids <- validation$entrez_mapping$ENTREZID[!is.na(validation$entrez_mapping$ENTREZID)]
                converted_entrez_genes(entrez_ids)
                cat("Stored", length(entrez_ids), "Entrez IDs from uploaded file for KEGG pathway operations\n")
            } else {
                converted_entrez_genes(character(0))
                cat("No valid Entrez IDs found from uploaded file\n")
            }
            
            showNotification(
                paste("Loaded", length(genes), input$gene_id_type, "IDs from file"),
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
        req(input$gene_id_type)
        
        if (!is.null(input$gene_text) && input$gene_text != "") {
            tryCatch({
                genes <- parse_gene_text(input$gene_text)
                uploaded_genes(genes)
                gene_id_type(input$gene_id_type)
                
                # Automatically validate genes with the selected ID type
                validation <- validate_gene_ids(genes, input$gene_id_type)
                gene_validation_results(validation)
                
                # Store converted Entrez IDs for KEGG pathway operations
                if (!is.null(validation$entrez_mapping) && nrow(validation$entrez_mapping) > 0) {
                    entrez_ids <- validation$entrez_mapping$ENTREZID[!is.na(validation$entrez_mapping$ENTREZID)]
                    converted_entrez_genes(entrez_ids)
                    cat("Stored", length(entrez_ids), "Entrez IDs for KEGG pathway operations\n")
                } else {
                    converted_entrez_genes(character(0))
                    cat("No valid Entrez IDs found for KEGG pathway operations\n")
                }
                
                showNotification(
                    paste("Loaded", length(genes), input$gene_id_type, "IDs from text"),
                    type = "message"
                )
            }, error = function(e) {
                showNotification(
                    paste("Error parsing gene IDs:", e$message),
                    type = "error"
                )
            })
        } else {
            showNotification(
                "Please enter some gene IDs first",
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
                data.frame(Message = "No gene IDs uploaded yet"),
                options = list(pageLength = 5, searching = FALSE, info = FALSE, paging = FALSE),
                rownames = FALSE
            ))
        }
        
        # Create a summary table of uploaded genes
        id_type <- gene_id_type()
        id_type_label <- switch(id_type,
                               "entrez" = "Entrez_ID",
                               "symbol" = "Gene_Symbol",
                               "ensembl" = "Ensembl_ID",
                               "uniprot" = "UniProt_ID",
                               "Gene_ID")
        
        gene_df <- data.frame(
            ID = genes,
            Type = id_type_label,
            Status = "Uploaded",
            stringsAsFactors = FALSE
        )
        names(gene_df)[1] <- id_type_label
        
        # Enhanced phylomap matching using the fixed map_genes_to_phylostrata function
        if (!is.null(phylomap_data)) {
            # Use the enhanced map_genes_to_phylostrata function that handles duplicates properly
            phylostrata_result <- map_genes_to_phylostrata(genes, id_type)
            
            # Check if genes were successfully mapped (have non-NA phylostrata)
            gene_df$In_Phylomap <- ifelse(
                !is.na(phylostrata_result), 
                "Yes", 
                "No"
            )
        } else {
            # If phylomap data is not available
            gene_df$In_Phylomap <- "Unknown"
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
        entrez_ids <- converted_entrez_genes()
        
        if (length(genes) == 0) {
            return("No genes uploaded")
        }
        
        id_type_label <- switch(gene_id_type(),
                               "entrez" = "Entrez IDs",
                               "symbol" = "Gene Symbols",
                               "ensembl" = "Ensembl IDs", 
                               "uniprot" = "UniProt IDs",
                               "Unknown")
        
        stats <- paste0(
            "Total genes: ", length(genes), " (", id_type_label, ")",
            " | Unique genes: ", length(unique(genes))
        )
        
        # Add Entrez conversion info if not already Entrez IDs
        if (gene_id_type() != "entrez" && length(entrez_ids) > 0) {
            conversion_rate <- round(length(entrez_ids) / length(genes) * 100, 1)
            stats <- paste0(stats, " | KEGG-ready: ", length(entrez_ids), " (", conversion_rate, "%)")
        }
        
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
        # Use both original genes and Entrez IDs for matching
        original_genes <- uploaded_genes()
        entrez_ids <- converted_entrez_genes()
        nodes <- values$nodes
        
        if (length(original_genes) == 0) {
            return("No genes selected")
        }
        
        if (is.null(nodes) || nrow(nodes) == 0) {
            return("No pathway loaded")
        }
        
        found_genes <- character(0)
        found_by_entrez <- character(0)
        
        # First, try to match using Entrez IDs (most reliable for KEGG pathways)
        if (length(entrez_ids) > 0 && !is.null(nodes$kegg_id)) {
            # Match Entrez IDs against KEGG IDs in pathway nodes
            for (i in seq_len(nrow(nodes))) {
                node_kegg_ids <- character(0)
                
                # Extract all Entrez IDs from this node
                if (!is.na(nodes$kegg_id[i]) && nodes$kegg_id[i] != "") {
                    node_kegg_ids <- c(node_kegg_ids, nodes$kegg_id[i])
                }
                
                # Also check gene_name field for Entrez IDs
                if (!is.null(nodes$gene_name) && !is.na(nodes$gene_name[i])) {
                    # Extract numeric IDs from entry names like "hsa:2475 hsa:57521"
                    numeric_ids <- regmatches(nodes$gene_name[i], gregexpr("\\d+", nodes$gene_name[i]))[[1]]
                    node_kegg_ids <- c(node_kegg_ids, numeric_ids)
                }
                
                # Check if any of our Entrez IDs match this node
                matches <- intersect(entrez_ids, node_kegg_ids)
                if (length(matches) > 0) {
                    # Map back to original gene symbols for display
                    validation <- gene_validation_results()
                    if (!is.null(validation$entrez_mapping)) {
                        for (entrez_match in matches) {
                            original_gene <- validation$entrez_mapping[
                                !is.na(validation$entrez_mapping$ENTREZID) & 
                                validation$entrez_mapping$ENTREZID == entrez_match, 
                                names(validation$entrez_mapping)[1]  # Get original ID column
                            ]
                            if (length(original_gene) > 0 && !is.na(original_gene[1])) {
                                found_by_entrez <- c(found_by_entrez, original_gene[1])
                            }
                        }
                    }
                }
            }
        }
        
        # Additionally, try traditional matching for any that weren't found by Entrez ID
        # This provides backward compatibility and catches edge cases
        genes_upper <- toupper(original_genes)
        traditional_found <- character(0)
        
        # Check HGNC symbols
        if (!is.null(nodes$hgnc_symbol)) {
            hgnc_matches <- original_genes[genes_upper %in% toupper(nodes$hgnc_symbol)]
            traditional_found <- c(traditional_found, hgnc_matches)
        }
        
        # Check labels (node labels)  
        if (!is.null(nodes$label)) {
            label_matches <- original_genes[genes_upper %in% toupper(nodes$label)]
            traditional_found <- c(traditional_found, label_matches)
        }
        
        # Combine all found genes and remove duplicates
        all_found <- unique(c(found_by_entrez, traditional_found))
        
        if (length(all_found) == 0) {
            id_type_label <- switch(gene_id_type(),
                                   "entrez" = "Entrez IDs",
                                   "symbol" = "Gene Symbols", 
                                   "ensembl" = "Ensembl IDs",
                                   "uniprot" = "UniProt IDs",
                                   "Gene IDs")
            
            return(paste0("None of your ", length(original_genes), " genes were found in this pathway.\n\n",
                         "Your genes (", id_type_label, "): ", paste(head(original_genes, 10), collapse = ", "), 
                         if (length(original_genes) > 10) "..." else "", "\n\n",
                         "💡 Tip: Gene matching works best when you upload Entrez Gene IDs.\n",
                         "Check the Gene Set tab to see conversion rates for your IDs."))
        }
        
        # Format the found genes
        genes_per_row <- 4
        gene_rows <- split(all_found, ceiling(seq_along(all_found) / genes_per_row))
        formatted_rows <- sapply(gene_rows, function(row) paste(row, collapse = "  |  "))
        
        missing_genes <- setdiff(original_genes, all_found)
        
        result <- paste(c(
            paste("🎯 FOUND:", length(all_found), "out of", length(original_genes), "genes in this pathway:"),
            "",
            formatted_rows
        ), collapse = "\n")
        
        if (length(missing_genes) > 0 && length(missing_genes) <= 10) {
            result <- paste0(result, "\n\n❌ Not found: ", paste(missing_genes, collapse = ", "))
        } else if (length(missing_genes) > 10) {
            result <- paste0(result, "\n\n❌ Not found: ", paste(missing_genes[1:10], collapse = ", "), "... (and ", length(missing_genes) - 10, " more)")
        }
        
        if (length(found_by_entrez) > 0) {
            result <- paste0(result, "\n\n✨ ", length(found_by_entrez), " genes matched using Entrez ID conversion")
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
        
        # Get selected pathway
        if (!is.null(input$pathway_search) && input$pathway_search != "") {
            filtered_pathways <- search_pathways(values$pathways_list, input$pathway_search)
            selected_row <- filtered_pathways[input$pathway_table_rows_selected, ]
        } else {
            category_pathways <- filter_pathways_by_category(values$pathways_list, input$pathway_category)
            selected_row <- category_pathways[input$pathway_table_rows_selected, ]
        }
        
        pathway_id <- selected_row$pathway_id
        
        # Show progress notification
        progress <- Progress$new()
        progress$set(message = paste("Loading pathway", pathway_id), value = 0.1)
        on.exit(progress$close())
        
        showNotification(paste("Loading pathway data:", pathway_id), type = "message")
        
        values$selected_pathway <- selected_row
        progress$set(value = 0.3)
        
        # Load pathway graph with HSA gene data
        tryCatch({
            progress$set(message = "Parsing pathway data...", value = 0.5)
            pathway_data <- parse_kegg_pathway_with_hsa(pathway_id, use_cached = TRUE)
            
            progress$set(message = "Building network...", value = 0.7)
            values$pathway_graph <- pathway_data$kegg_data
            values$nodes <- pathway_data$nodes
            values$edges <- pathway_data$edges
            
            progress$set(message = "Finalizing...", value = 0.9)
            
            showNotification(
                paste("Pathway loaded successfully with", nrow(pathway_data$nodes), "network nodes!"), 
                type = "message"
            )
            
            progress$set(value = 1.0)
            
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
        
        # Use Entrez IDs for gene highlighting instead of original IDs
        # This ensures genes are found in KEGG pathways regardless of input ID type
        highlight_genes_to_use <- NULL
        if (coloring_mode == "highlight_genes") {
            entrez_ids <- converted_entrez_genes()
            if (length(entrez_ids) > 0) {
                highlight_genes_to_use <- entrez_ids
                cat("Using", length(entrez_ids), "Entrez IDs for gene highlighting in pathway\n")
            } else {
                # Fallback to original genes if no Entrez conversion available
                highlight_genes_to_use <- uploaded_genes()
                cat("No Entrez IDs available, using original", length(highlight_genes_to_use), "IDs for highlighting\n")
            }
        }
        
        create_kegg_network_visualization(
            values$nodes, 
            values$edges,
            show_labels = TRUE,
            show_edges = TRUE,
            coloring_mode = coloring_mode,
            highlight_genes = highlight_genes_to_use,
            id_type = gene_id_type()  # Pass the current gene ID type
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
            
            output$node_info <- renderUI({
                if (nrow(selected_node) > 0) {
                    # Find edges connected to this node
                    node_id <- selected_node$id
                    incoming_edges <- values$edges[values$edges$to == node_id, ]
                    outgoing_edges <- values$edges[values$edges$from == node_id, ]
                    
                    # Get comprehensive gene ID information
                    hgnc_symbol <- selected_node$hgnc_symbol
                    kegg_id <- selected_node$kegg_id
                    gene_name <- selected_node$gene_name
                    display_label <- selected_node$label
                    node_type <- selected_node$type
                    
                    # Initialize comprehensive gene ID info
                    gene_id_info <- ""
                    uniprot_info <- NULL
                    uniprot_link <- ""
                    
                    # Get comprehensive mapping data for this gene
                    if (!is.null(hgnc_symbol) && length(hgnc_symbol) > 0 && !is.na(hgnc_symbol) && hgnc_symbol != "") {
                        tryCatch({
                            # Try to get comprehensive mapping from the loaded data
                            if (!is.null(comprehensive_mapping)) {
                                # Look up by HGNC symbol
                                gene_info <- comprehensive_mapping[
                                    !is.na(comprehensive_mapping$hgnc_symbol) & 
                                    comprehensive_mapping$hgnc_symbol == hgnc_symbol, 
                                ]
                                
                                if (nrow(gene_info) > 0) {
                                    gene_record <- gene_info[1, ]  # Take first match if multiple
                                    
                                    # Build comprehensive ID information
                                    gene_id_info <- paste0(
                                        "<strong>=== GENE IDENTIFIERS ===</strong><br>",
                                        if (!is.null(gene_record$hgnc_symbol) && !is.na(gene_record$hgnc_symbol) && gene_record$hgnc_symbol != "") 
                                            paste0("<strong>Gene Symbol (HGNC):</strong> <code>", gene_record$hgnc_symbol, "</code><br>") else "",
                                        if (!is.null(gene_record$entrezgene_id) && !is.na(gene_record$entrezgene_id) && gene_record$entrezgene_id != "") 
                                            paste0("<strong>Entrez Gene ID:</strong> <code>", gene_record$entrezgene_id, "</code><br>") else "",
                                        if (!is.null(gene_record$ensembl_gene_id) && !is.na(gene_record$ensembl_gene_id) && gene_record$ensembl_gene_id != "") 
                                            paste0("<strong>Ensembl Gene ID:</strong> <code>", gene_record$ensembl_gene_id, "</code><br>") else "",
                                        if (!is.null(gene_record$uniprotswissprot) && !is.na(gene_record$uniprotswissprot) && gene_record$uniprotswissprot != "") 
                                            paste0("<strong>UniProt ID:</strong> <code>", gene_record$uniprotswissprot, "</code><br>") else "",
                                        "<strong>KEGG Gene ID:</strong> <code>", kegg_id, "</code><br>",
                                        "<strong>Original KEGG Entry:</strong> <code>", gene_name, "</code><br>",
                                        "<strong>Node Type:</strong> ", node_type, "<br>",
                                        "<strong>Display Label:</strong> ", display_label, "<br>"
                                    )
                                    
                                    # Set uniprot_info for additional processing
                                    if (!is.null(gene_record$uniprotswissprot) && !is.na(gene_record$uniprotswissprot) && gene_record$uniprotswissprot != "") {
                                        uniprot_info <- gene_record
                                        # Create structure link if possible
                                        uniprot_link <- create_uniprot_structure_link(gene_record$uniprotswissprot)
                                    }
                                }
                            }
                        }, error = function(e) {
                            cat("Error getting comprehensive mapping:", e$message, "\n")
                        })
                    }
                    
                    # Fallback to basic information if comprehensive mapping failed
                    if (gene_id_info == "") {
                        gene_id_info <- paste0(
                            "<strong>=== GENE INFORMATION ===</strong><br>",
                            "<strong>Gene Symbol:</strong> <code>", if(!is.null(hgnc_symbol) && !is.na(hgnc_symbol)) hgnc_symbol else "Not available", "</code><br>",
                            "<strong>KEGG Gene ID:</strong> <code>", kegg_id, "</code><br>",
                            "<strong>Original KEGG Entry:</strong> <code>", gene_name, "</code><br>",
                            "<strong>Node Type:</strong> ", node_type, "<br>",
                            "<strong>Display Label:</strong> ", display_label, "<br>"
                        )
                        
                        # Try simpler UniProt lookup as fallback
                        if (!is.null(hgnc_symbol) && !is.na(hgnc_symbol) && hgnc_symbol != "") {
                            uniprot_id <- get_uniprot_id(hgnc_symbol, comprehensive_mapping)
                            if (!is.null(uniprot_id)) {
                                uniprot_link <- create_uniprot_structure_link(uniprot_id)
                            }
                        }
                    }
                    
                    # Get phylostratum information if available
                    phylo_info <- ""
                    if (!is.null(selected_node$phylostratum) && !is.na(selected_node$phylostratum)) {
                        phylo_stratum <- selected_node$phylostratum
                        # Load phylostratum legend to get the name
                        legend_data <- tryCatch({
                            source("R/kegg_utils.R", local = TRUE)
                            load_phylostratum_legend()
                        }, error = function(e) NULL)
                        
                        if (!is.null(legend_data)) {
                            phylo_name <- legend_data$Name[legend_data$Rank == phylo_stratum]
                            if (length(phylo_name) > 0) {
                                phylo_info <- paste0(
                                    "<strong>Phylostratum:</strong> ", phylo_stratum, " - ", phylo_name[1], "<br>"
                                )
                            } else {
                                phylo_info <- paste0("<strong>Phylostratum:</strong> ", phylo_stratum, "<br>")
                            }
                        } else {
                            phylo_info <- paste0("<strong>Phylostratum:</strong> ", phylo_stratum, "<br>")
                        }
                    }
                    
                    # Build enhanced UniProt information section
                    uniprot_section <- ""
                    if (!is.null(uniprot_info)) {
                        # Build UniProt section content (excluding ID which is already shown above)
                        uniprot_content_parts <- c()
                        
                        # Add protein name if available
                        if (!is.null(uniprot_info$name) && !is.na(uniprot_info$name) && uniprot_info$name != "") {
                            uniprot_content_parts <- c(uniprot_content_parts, 
                                paste0("<strong>Protein Name:</strong> ", uniprot_info$name, "<br>"))
                        }
                        
                        # Add aliases if available
                        if (!is.null(uniprot_info$alias_symbol) && !is.na(uniprot_info$alias_symbol) && uniprot_info$alias_symbol != "") {
                            uniprot_content_parts <- c(uniprot_content_parts, 
                                paste0("<strong>Aliases:</strong> ", uniprot_info$alias_symbol, "<br>"))
                        }
                        
                        # Add structure link if available
                        if (uniprot_link != "") {
                            uniprot_content_parts <- c(uniprot_content_parts, 
                                paste0("<strong>3D Structure:</strong> ", uniprot_link, "<br>"))
                        }
                        
                        # Only create the UniProt section if we have additional information beyond the ID
                        if (length(uniprot_content_parts) > 0) {
                            uniprot_section <- paste0(
                                "<br><strong>=== PROTEIN INFORMATION ===</strong><br>",
                                paste(uniprot_content_parts, collapse = "")
                            )
                        }
                    } else if (uniprot_link != "") {
                        # If we only have a structure link, show minimal section
                        uniprot_section <- paste0(
                            "<br><strong>=== PROTEIN INFORMATION ===</strong><br>",
                            "<strong>3D Structure:</strong> ", uniprot_link, "<br>"
                        )
                    }
                    
                    # Build relationship info
                    relationship_info <- ""
                    if (nrow(incoming_edges) > 0 || nrow(outgoing_edges) > 0) {
                        relationship_info <- paste0(
                            "<br><strong>=== RELATIONSHIPS ===</strong><br>",
                            "Incoming connections: ", nrow(incoming_edges), "<br>",
                            "Outgoing connections: ", nrow(outgoing_edges), "<br>"
                        )
                        
                        # Show some specific relationships
                        if (nrow(incoming_edges) > 0) {
                            sample_incoming <- head(incoming_edges, 3)
                            for (i in seq_len(nrow(sample_incoming))) {
                                edge <- sample_incoming[i, ]
                                rel_type <- if (!is.null(edge$relation_type) && !is.na(edge$relation_type)) edge$relation_type else "Unknown"
                                subtype <- if (!is.null(edge$subtype) && !is.na(edge$subtype) && edge$subtype != "") edge$subtype else ""
                                
                                # Find the source node name - prefer symbol over ID
                                source_node <- values$nodes[values$nodes$id == edge$from, ]
                                source_name <- edge$from  # default fallback
                                if (nrow(source_node) > 0) {
                                    if (!is.null(source_node$hgnc_symbol) && !is.na(source_node$hgnc_symbol) && source_node$hgnc_symbol != "") {
                                        source_name <- source_node$hgnc_symbol
                                    } else if (!is.null(source_node$label) && !is.na(source_node$label) && source_node$label != "") {
                                        source_name <- source_node$label
                                    }
                                }
                                
                                relationship_info <- paste0(relationship_info, 
                                    "← ", source_name, " (", rel_type, 
                                    if (subtype != "") paste0(": ", subtype) else "", ")", "<br>")
                            }
                        }
                        
                        if (nrow(outgoing_edges) > 0) {
                            sample_outgoing <- head(outgoing_edges, 3)
                            for (i in seq_len(nrow(sample_outgoing))) {
                                edge <- sample_outgoing[i, ]
                                rel_type <- if (!is.null(edge$relation_type) && !is.na(edge$relation_type)) edge$relation_type else "Unknown"
                                subtype <- if (!is.null(edge$subtype) && !is.na(edge$subtype) && edge$subtype != "") edge$subtype else ""
                                
                                # Find the target node name - prefer symbol over ID
                                target_node <- values$nodes[values$nodes$id == edge$to, ]
                                target_name <- edge$to  # default fallback
                                if (nrow(target_node) > 0) {
                                    if (!is.null(target_node$hgnc_symbol) && !is.na(target_node$hgnc_symbol) && target_node$hgnc_symbol != "") {
                                        target_name <- target_node$hgnc_symbol
                                    } else if (!is.null(target_node$label) && !is.na(target_node$label) && target_node$label != "") {
                                        target_name <- target_node$label
                                    }
                                }
                                
                                relationship_info <- paste0(relationship_info, 
                                    "→ ", target_name, " (", rel_type, 
                                    if (subtype != "") paste0(": ", subtype) else "", ")", "<br>")
                            }
                        }
                    }
                    
                    # Create comprehensive HTML content
                    html_content <- paste0(
                        gene_id_info,
                        if (phylo_info != "") paste0("<br>", phylo_info) else "",
                        uniprot_section,
                        relationship_info, 
                        "<br><strong>=== PATHWAY CONTEXT ===</strong><br>",
                        "In KEGG, this gene is referenced as: <code>", kegg_id, "</code><br>",
                        if (!is.null(hgnc_symbol) && !is.na(hgnc_symbol) && hgnc_symbol != "") 
                            paste0("Common name/symbol: <code>", hgnc_symbol, "</code><br>") else "",
                        "<br>",
                        "<em>Tip: Search for '", kegg_id, "' in KEGG database<br>",
                        if (!is.null(hgnc_symbol) && !is.na(hgnc_symbol) && hgnc_symbol != "") 
                            paste0("or '", hgnc_symbol, "' in gene databases</em>") else "</em>"
                    )
                    
                    # Return as HTML
                    HTML(html_content)
                } else {
                    HTML("<em>No node selected - click on a node to see detailed gene information</em>")
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
            return("Click 'Load Genes' to validate gene IDs")
        }
        
        total_genes <- length(validation$valid_genes) + length(validation$invalid_genes)
        conversion_rate <- round(validation$conversion_rate * 100, 1)
        id_type_label <- switch(validation$id_type,
                               "entrez" = "Entrez IDs",
                               "symbol" = "Gene Symbols",
                               "ensembl" = "Ensembl IDs", 
                               "uniprot" = "UniProt IDs",
                               "Gene IDs")
        
        # Count converted Entrez IDs
        entrez_count <- 0
        if (!is.null(validation$entrez_mapping)) {
            entrez_count <- sum(!is.na(validation$entrez_mapping$ENTREZID))
        }
        
        result <- paste0(
            "GENE ID VALIDATION RESULTS\n",
            "==========================\n",
            "Input Type: ", id_type_label, "\n",
            "Total IDs: ", total_genes, "\n",
            "Valid IDs: ", length(validation$valid_genes), "\n",
            "Invalid IDs: ", length(validation$invalid_genes), "\n",
            "Validation rate: ", conversion_rate, "%\n"
        )
        
        if (validation$id_type != "entrez") {
            result <- paste0(result, "\n🎯 KEGG PATHWAY COMPATIBILITY:\n")
            result <- paste0(result, "Converted to Entrez IDs: ", entrez_count, "\n")
            result <- paste0(result, "KEGG conversion rate: ", conversion_rate, "%\n")
            result <- paste0(result, "\n💡 ", entrez_count, " genes available for KEGG pathway analysis\n")
        } else {
            result <- paste0(result, "\n✅ Using Entrez IDs (optimal for KEGG pathways)\n")
        }
        
        if (length(validation$invalid_genes) > 0 && length(validation$invalid_genes) <= 10) {
            result <- paste0(result, "\n❌ Could not convert:\n", 
                           paste(validation$invalid_genes, collapse = ", "), "\n")
            result <- paste0(result, "\n💡 These genes won't be highlighted in pathway visualizations")
        } else if (length(validation$invalid_genes) > 10) {
            result <- paste0(result, "\n❌ Could not convert (showing first 10):\n", 
                           paste(validation$invalid_genes[1:10], collapse = ", "), "...\n")
            result <- paste0(result, "\n💡 ", length(validation$invalid_genes), " genes won't be highlighted in pathways")
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
        original_genes <- uploaded_genes()
        entrez_ids <- converted_entrez_genes()
        
        if (length(original_genes) == 0) {
            showNotification("Please load genes first", type = "warning")
            return()
        }
        
        if (length(entrez_ids) == 0) {
            showNotification("No valid Entrez IDs found for enrichment analysis. Check your gene IDs.", type = "error")
            return()
        }
        
        # Show progress notification
        progress <- Progress$new()
        progress$set(message = "Starting KEGG enrichment analysis...", value = 0.1)
        on.exit(progress$close())
        
        showNotification(paste("Running KEGG enrichment with", length(entrez_ids), "Entrez IDs..."), type = "message")
        
        progress$set(message = paste("Using", length(entrez_ids), "converted Entrez IDs..."), value = 0.3)
        
        # Run enrichment analysis directly with Entrez IDs for optimal performance
        result <- perform_kegg_enrichment(
            entrez_ids,  # Use converted Entrez IDs directly
            id_type = "entrez",  # Always use Entrez for KEGG analysis
            organism = "hsa",
            pvalue_cutoff = input$pvalue_cutoff,
            qvalue_cutoff = input$qvalue_cutoff,
            progress = progress
        )
        
        progress$set(message = "Processing results...", value = 0.9)
        
        enrichment_results(result)
        
        progress$set(value = 1.0)
        
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
                pathway_data <- parse_kegg_pathway_with_hsa(pathway_id, use_cached = TRUE)
                
                progress$set(message = "Building network...", value = 0.7)
                values$pathway_graph <- pathway_data$kegg_data
                values$nodes <- pathway_data$nodes
                values$edges <- pathway_data$edges
                
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
