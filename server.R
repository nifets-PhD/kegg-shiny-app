# Define Server
server <- function(input, output, session) {
    
    # Reactive values to store data
    values <- reactiveValues(
        pathways_list = NULL,
        selected_pathway = NULL,
        pathway_graph = NULL,
        nodes = NULL,
        edges = NULL
    )
    
    # Centralized gene management - all genes stored internally as Entrez IDs
    # Load centralized gene ID management system
    source("R/gene_id_manager.R")
    
    # Reactive values for gene management
    internal_genes <- reactiveVal(character(0))  # Always Entrez IDs (internal standard)
    gene_mapping_info <- reactiveVal(data.frame())  # Mapping info for display
    gene_conversion_stats <- reactiveVal(list())  # Conversion statistics
    
    # Helper functions for backward compatibility with existing code
    # These return Entrez IDs (the internal standard) for all operations
    uploaded_genes <- reactive({
        # Returns Entrez IDs (internal standard)
        internal_genes()
    })
    
    converted_entrez_genes <- reactive({
        # Returns Entrez IDs (internal standard) - same as uploaded_genes
        internal_genes()
    })
    
    gene_id_type <- reactive({
        # Always return "entrez" since internally everything is Entrez
        "entrez"
    })
    
    # Function to get display symbols for user-friendly output
    get_display_symbols <- reactive({
        entrez_ids <- internal_genes()
        if (length(entrez_ids) == 0) return(character(0))
        return(entrez_to_symbols(entrez_ids, comprehensive_mapping))
    })
    
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
        comprehensive_mapping <- load_comprehensive_gene_mapping()
    }, error = function(e) {
        cat("Error loading comprehensive gene mapping:", e$message, "\n")
        comprehensive_mapping <- NULL
    })
    
    # Load KEGG pathways on app start
    observe({
        showNotification("Loading KEGG pathways...", type = "message")
        values$pathways_list <- get_kegg_pathways()
        showNotification("KEGG pathways loaded successfully!", type = "message")
    })
    
    # Handle file upload
    observeEvent(input$gene_file, {
        req(input$gene_file, input$gene_id_type)
        
        tryCatch({
            # Load raw genes from file
            raw_genes <- load_gene_file(input$gene_file$datapath)
            
            # Convert to internal standard (Entrez IDs) immediately
            conversion_result <- convert_to_internal_standard(
                gene_ids = raw_genes,
                input_type = input$gene_id_type,
                comprehensive_mapping = comprehensive_mapping
            )
            
            # Store results in reactive values
            internal_genes(conversion_result$entrez_ids)
            gene_mapping_info(conversion_result$mapping_info)
            gene_conversion_stats(conversion_result$conversion_stats)
            
            # Show notification with conversion stats
            stats <- conversion_result$conversion_stats
            showNotification(
                paste("Loaded", stats$input_count, input$gene_id_type, "IDs from file.",
                      "Converted", stats$converted_count, "to Entrez IDs",
                      paste0("(", round(stats$conversion_rate * 100, 1), "% success rate)")),
                type = "message"
            )
            
            cat("Gene upload summary:\n")
            cat("  Input:", stats$input_count, input$gene_id_type, "IDs\n")
            cat("  Converted:", stats$converted_count, "to Entrez IDs\n")
            cat("  Success rate:", round(stats$conversion_rate * 100, 1), "%\n")
            
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
            internal_genes(character(0))
            gene_mapping_info(data.frame())
            gene_conversion_stats(list())
        }
    })
    
    # Handle Load Genes button
    observeEvent(input$load_genes_text, {
        req(input$gene_id_type)
        
        if (!is.null(input$gene_text) && input$gene_text != "") {
            tryCatch({
                # Parse genes from text
                raw_genes <- parse_gene_text(input$gene_text)
                
                # Convert to internal standard (Entrez IDs) immediately
                conversion_result <- convert_to_internal_standard(
                    gene_ids = raw_genes,
                    input_type = input$gene_id_type,
                    comprehensive_mapping = comprehensive_mapping
                )
                
                # Store results in reactive values
                internal_genes(conversion_result$entrez_ids)
                gene_mapping_info(conversion_result$mapping_info)
                gene_conversion_stats(conversion_result$conversion_stats)
                
                # Show notification with conversion stats
                stats <- conversion_result$conversion_stats
                showNotification(
                    paste("Loaded", stats$input_count, input$gene_id_type, "IDs from text.",
                          "Converted", stats$converted_count, "to Entrez IDs",
                          paste0("(", round(stats$conversion_rate * 100, 1), "% success rate)")),
                    type = "message"
                )
                
                cat("Text gene input summary:\n")
                cat("  Input:", stats$input_count, input$gene_id_type, "IDs\n")
                cat("  Converted:", stats$converted_count, "to Entrez IDs\n")
                cat("  Success rate:", round(stats$conversion_rate * 100, 1), "%\n")
            }, error = function(e) {
                showNotification(
                    paste("Error processing genes:", e$message),
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
        internal_genes(character(0))
        gene_mapping_info(data.frame())
        gene_conversion_stats(list())
        updateTextAreaInput(session, "gene_text", value = "")
        
        showNotification(
            "Cleared uploaded genes",
            type = "message",
            duration = 2
        )
    })
    
    # Display uploaded genes table
    output$loaded_genes_table <- renderDT({
        entrez_genes <- internal_genes()
        mapping_info <- gene_mapping_info()
        
        if (length(entrez_genes) == 0) {
            return(datatable(
                data.frame(Message = "No gene IDs uploaded yet"),
                options = list(pageLength = 5, searching = FALSE, info = FALSE, paging = FALSE),
                rownames = FALSE
            ))
        }
        
        # Create display table with original IDs and their Entrez mappings
        if (nrow(mapping_info) > 0) {
            gene_df <- data.frame(
                Original_ID = mapping_info$input_id,
                Input_Type = mapping_info$input_type,
                Entrez_ID = mapping_info$entrez_id,
                Gene_Symbol = mapping_info$symbol,
                Status = "Converted",
                stringsAsFactors = FALSE
            )
        } else {
            # Fallback if no mapping info (shouldn't happen with centralized system)
            gene_df <- data.frame(
                Original_ID = entrez_genes,
                Input_Type = "entrez",
                Entrez_ID = entrez_genes,
                Gene_Symbol = entrez_genes,
                Status = "Direct",
                stringsAsFactors = FALSE
            )
        }
        
        # Enhanced phylomap matching using Entrez IDs (internal standard)
        if (!is.null(phylomap_data)) {
            # Convert Entrez IDs to gene symbols for phylomap matching
            gene_symbols <- entrez_to_symbols(entrez_genes, comprehensive_mapping)
            phylostrata_result <- map_genes_to_phylostrata(gene_symbols, "symbol")
            
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
                pageLength = 50,
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
        entrez_genes <- internal_genes()
        mapping_info <- gene_mapping_info()
        stats <- gene_conversion_stats()
        
        if (length(entrez_genes) == 0) {
            return("No genes uploaded")
        }
        
        # Create comprehensive statistics using centralized data
        if (length(stats) > 0) {
            input_type_label <- if (nrow(mapping_info) > 0) {
                unique_types <- unique(mapping_info$input_type)
                if (length(unique_types) == 1) {
                    switch(unique_types[1],
                        "entrez" = "Entrez IDs",
                        "symbol" = "Gene Symbols", 
                        "ensembl" = "Ensembl IDs",
                        "uniprot" = "UniProt IDs",
                        "Mixed types"
                    )
                } else {
                    "Mixed types"
                }
            } else {
                "Unknown"
            }
            
            stat_text <- paste0(
                "Input: ", stats$input_count, " (", input_type_label, ")",
                " | Converted to Entrez: ", stats$converted_count, " (", round(stats$conversion_rate * 100, 1), "%)",
                " | KEGG-ready: ", length(entrez_genes)
            )
            
            # Add phylomap overlap if available
            if (!is.null(phylomap_data) && length(entrez_genes) > 0) {
                gene_symbols <- entrez_to_symbols(entrez_genes, comprehensive_mapping)
                phylomap_genes <- toupper(phylomap_data$GeneID)
                overlap <- sum(toupper(gene_symbols) %in% phylomap_genes)
                stat_text <- paste0(stat_text, " | In phylomap: ", overlap, " (", round(overlap/length(entrez_genes)*100, 1), "%)")
            }
            
            return(stat_text)
        } else {
            return(paste0("Internal genes: ", length(entrez_genes), " (Entrez IDs)"))
        }
    })
    
    # Helper function to get readable gene labels from network nodes
    get_readable_gene_labels <- function(input_genes, nodes = NULL) {
        if (length(input_genes) == 0) return(character(0))
        
        # If no pathway is loaded, return original genes
        if (is.null(nodes) || nrow(nodes) == 0) {
            return(input_genes)
        }
        
        readable_labels <- character(length(input_genes))
        entrez_ids <- converted_entrez_genes()
        
        for (i in seq_along(input_genes)) {
            gene <- input_genes[i]
            readable_labels[i] <- gene  # Default fallback
            
            # Try to find this gene in the network nodes
            # Method 1: Match by Entrez ID if available
            if (length(entrez_ids) >= i && !is.na(entrez_ids[i])) {
                for (j in seq_len(nrow(nodes))) {
                    node_kegg_ids <- character(0)
                    
                    # Extract Entrez IDs from this node
                    if (!is.na(nodes$kegg_id[j]) && nodes$kegg_id[j] != "") {
                        node_kegg_ids <- c(node_kegg_ids, nodes$kegg_id[j])
                    }
                    
                    if (!is.null(nodes$gene_name) && !is.na(nodes$gene_name[j])) {
                        numeric_ids <- regmatches(nodes$gene_name[j], gregexpr("\\d+", nodes$gene_name[j]))[[1]]
                        node_kegg_ids <- c(node_kegg_ids, numeric_ids)
                    }
                    
                    if (entrez_ids[i] %in% node_kegg_ids) {
                        if (!is.null(nodes$label) && !is.na(nodes$label[j]) && nodes$label[j] != "") {
                            readable_labels[i] <- nodes$label[j]
                            break
                        }
                    }
                }
            }
            
            # Method 2: Try direct matching by gene symbol/label
            if (readable_labels[i] == gene) {  # Still using fallback, try direct match
                gene_upper <- toupper(gene)
                
                # Check HGNC symbols
                if (!is.null(nodes$hgnc_symbol)) {
                    match_idx <- which(toupper(nodes$hgnc_symbol) == gene_upper)
                    if (length(match_idx) > 0 && !is.null(nodes$label) && 
                        !is.na(nodes$label[match_idx[1]]) && nodes$label[match_idx[1]] != "") {
                        readable_labels[i] <- nodes$label[match_idx[1]]
                    }
                }
                
                # Check labels directly
                if (readable_labels[i] == gene && !is.null(nodes$label)) {
                    match_idx <- which(toupper(nodes$label) == gene_upper)
                    if (length(match_idx) > 0) {
                        readable_labels[i] <- nodes$label[match_idx[1]]
                    }
                }
            }
        }
        
        return(readable_labels)
    }
    
    # Selected genes display for Pathway Explorer tab
    output$selected_genes_display <- renderText({
        genes <- uploaded_genes()
        if (length(genes) == 0) {
            return("No genes selected yet. Upload a file or enter genes in the text box above.")
        }
        
        # Get readable labels (will return original if no pathway loaded)
        readable_genes <- get_readable_gene_labels(genes, values$nodes)
        
        # Format genes in columns for better display
        genes_per_row <- 5
        gene_rows <- split(readable_genes, ceiling(seq_along(readable_genes) / genes_per_row))
        formatted_rows <- sapply(gene_rows, function(row) paste(row, collapse = "  |  "))
        
        paste(c(
            paste("Selected Genes (", length(genes), " total):"),
            "",
            formatted_rows
        ), collapse = "\n")
    })
    
    # Selected genes display for Network tab
    output$selected_genes_display_network <- renderText({
        genes <- uploaded_genes()
        if (length(genes) == 0) {
            return("No genes selected. Go to Pathway Explorer tab to upload genes.")
        }
        
        # Get readable labels from network nodes
        readable_genes <- get_readable_gene_labels(genes, values$nodes)
        
        # Format genes in columns for better display
        genes_per_row <- 6
        gene_rows <- split(readable_genes, ceiling(seq_along(readable_genes) / genes_per_row))
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
        
        found_genes_with_labels <- list()  # Store both original gene and readable label
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
                    # Map back to original gene symbols for display, but use node label for readable display
                    validation <- gene_validation_results()
                    node_label <- if (!is.null(nodes$label) && !is.na(nodes$label[i]) && nodes$label[i] != "") {
                        nodes$label[i]
                    } else if (!is.null(nodes$hgnc_symbol) && !is.na(nodes$hgnc_symbol[i]) && nodes$hgnc_symbol[i] != "") {
                        nodes$hgnc_symbol[i]
                    } else {
                        matches[1]  # Fallback to Entrez ID
                    }
                    
                    if (!is.null(validation$entrez_mapping)) {
                        for (entrez_match in matches) {
                            original_gene <- validation$entrez_mapping[
                                !is.na(validation$entrez_mapping$ENTREZID) & 
                                validation$entrez_mapping$ENTREZID == entrez_match, 
                                names(validation$entrez_mapping)[1]  # Get original ID column
                            ]
                            if (length(original_gene) > 0 && !is.na(original_gene[1])) {
                                found_genes_with_labels[[original_gene[1]]] <- node_label
                                found_by_entrez <- c(found_by_entrez, original_gene[1])
                            }
                        }
                    }
                }
            }
        }
        
        # Additionally, try traditional matching for any that weren't found by Entrez ID
        genes_upper <- toupper(original_genes)
        traditional_found <- character(0)
        
        # Check HGNC symbols
        if (!is.null(nodes$hgnc_symbol)) {
            for (j in seq_along(original_genes)) {
                gene <- original_genes[j]
                if (!(gene %in% names(found_genes_with_labels))) {  # Not already found
                    match_idx <- which(toupper(nodes$hgnc_symbol) == toupper(gene))
                    if (length(match_idx) > 0) {
                        node_label <- if (!is.null(nodes$label) && !is.na(nodes$label[match_idx[1]]) && nodes$label[match_idx[1]] != "") {
                            nodes$label[match_idx[1]]
                        } else {
                            nodes$hgnc_symbol[match_idx[1]]
                        }
                        found_genes_with_labels[[gene]] <- node_label
                        traditional_found <- c(traditional_found, gene)
                    }
                }
            }
        }
        
        # Check labels (node labels)  
        if (!is.null(nodes$label)) {
            for (j in seq_along(original_genes)) {
                gene <- original_genes[j]
                if (!(gene %in% names(found_genes_with_labels))) {  # Not already found
                    match_idx <- which(toupper(nodes$label) == toupper(gene))
                    if (length(match_idx) > 0) {
                        found_genes_with_labels[[gene]] <- nodes$label[match_idx[1]]
                        traditional_found <- c(traditional_found, gene)
                    }
                }
            }
        }
        
        # Get all found genes
        all_found <- names(found_genes_with_labels)
        
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
        
        # Format the found genes using readable labels
        readable_labels <- unlist(found_genes_with_labels)
        genes_per_row <- 4
        gene_rows <- split(readable_labels, ceiling(seq_along(readable_labels) / genes_per_row))
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
    
    # Network highlighted genes display (for highlight_genes option)
    output$network_highlighted_genes <- renderDT({
        entrez_ids <- internal_genes()
        if (length(entrez_ids) == 0) {
            return(datatable(
                data.frame(Message = "No genes loaded"),
                options = list(pageLength = 5, searching = FALSE, info = FALSE, 
                              paging = FALSE, dom = 't'),
                rownames = FALSE
            ))
        }
        
        # Convert internal Entrez IDs to symbols for display
        symbols <- entrez_to_symbols(entrez_ids, comprehensive_mapping)
        
        # Calculate columns: min(num_genes, 5)
        num_cols <- min(length(symbols), 5)
        num_rows <- ceiling(length(symbols) / num_cols)
        
        # Create matrix to arrange genes
        gene_matrix <- matrix(NA, nrow = num_rows, ncol = num_cols)
        for (i in seq_along(symbols)) {
            row <- ((i - 1) %/% num_cols) + 1
            col <- ((i - 1) %% num_cols) + 1
            gene_matrix[row, col] <- symbols[i]
        }
        
        # Convert to data frame
        gene_df <- as.data.frame(gene_matrix, stringsAsFactors = FALSE)
        colnames(gene_df) <- paste0("Col", 1:num_cols)
        
        # Replace NA with empty strings
        gene_df[is.na(gene_df)] <- ""
        
        datatable(
            gene_df,
            options = list(
                pageLength = 25, 
                searching = FALSE, 
                info = FALSE, 
                paging = FALSE,
                dom = 't',
                columnDefs = list(list(className = 'dt-center', targets = '_all'))
            ),
            rownames = FALSE,
            colnames = rep("", ncol(gene_df))
        )
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
    
    # Load pathway from node selection (map nodes)
    observeEvent(input$load_pathway_from_node, {
        req(input$load_pathway_from_node)
        
        pathway_id <- input$load_pathway_from_node
        
        # Show progress notification
        progress <- Progress$new()
        progress$set(message = paste("Loading pathway", pathway_id, "from node"), value = 0.1)
        on.exit(progress$close())
        
        showNotification(paste("Loading pathway from node:", pathway_id), type = "message")
        
        # Create a fake pathway entry for the selected pathway
        # We need to have the pathway in the pathways list to load it properly
        # For now, we'll create a minimal entry
        fake_pathway_entry <- data.frame(
            pathway_id = pathway_id,
            pathway_name = paste("Pathway", pathway_id),
            description = paste("Loaded from map node -", pathway_id),
            stringsAsFactors = FALSE
        )
        
        values$selected_pathway <- fake_pathway_entry
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
                paste("Pathway", pathway_id, "loaded successfully with", nrow(pathway_data$nodes), "network nodes!"), 
                type = "message"
            )
            
            progress$set(value = 1.0)
            
            # Switch to network tab to show the loaded pathway
            updateTabItems(session, "tabs", "network")
            
        }, error = function(e) {
            showNotification(paste("Error loading pathway from node:", e$message), type = "warning")
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
            
            # Create the data table with color styling
            DT::datatable(
                legend_display,
                options = list(
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
                ), 
                rownames = FALSE,
                escape = FALSE
            ) %>%
            DT::formatStyle(
                'Type',
                backgroundColor = DT::styleEqual(
                    legend_data$relationship,
                    legend_data$color
                ),
                color = DT::styleEqual(
                    legend_data$relationship,
                    sapply(legend_data$color, function(color) {
                        if (is.na(color)) return("#000000")
                        # Calculate brightness to determine text color
                        rgb_vals <- col2rgb(color)
                        brightness <- (rgb_vals[1] * 0.299 + rgb_vals[2] * 0.587 + rgb_vals[3] * 0.114) / 255
                        if (brightness < 0.5) "#FFFFFF" else "#000000"
                    })
                ),
                fontWeight = "bold"
            )
        } else {
            DT::datatable(
                data.frame(
                    Type = "No relationships",
                    Description = "No edge relationships found in this pathway",
                    stringsAsFactors = FALSE
                ),
                options = list(
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
                ), 
                rownames = FALSE,
                escape = FALSE
            )
        }
    })
    
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
                    
                    # DEBUG: Check what types are available in values$nodes
                    cat("DEBUG: All node types in values$nodes:\n")
                    print(table(values$nodes$type, useNA = "always"))
                    cat("DEBUG: Selected node ID:", input$pathway_network_selected, "\n")
                    cat("DEBUG: Selected node data:\n")
                    print(selected_node)
                    
                    # Get basic node information
                    hgnc_symbol <- selected_node$hgnc_symbol
                    kegg_id <- selected_node$kegg_id
                    gene_name <- selected_node$gene_name
                    display_label <- selected_node$label
                    node_type <- selected_node$type
                    cat("DEBUG: Extracted node_type:", node_type, "\n")
                    print(node_type)
                    
                    # Handle different node types
                    if (node_type == "gene") {
                        # ======= GENE NODE HANDLING =======
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
                        
                        # Add EMERALD alignment section if UniProt ID is available
                        if (!is.null(uniprot_info$uniprotswissprot) && !is.na(uniprot_info$uniprotswissprot) && uniprot_info$uniprotswissprot != "") {
                            # Store the selected gene's UniProt ID
                            selected_gene_uniprot(uniprot_info$uniprotswissprot)
                            
                            # Get genes in pathway with UniProt IDs for the dropdown
                            pathway_genes_with_uniprot <- c()
                            if (!is.null(values$nodes) && nrow(values$nodes) > 0 && !is.null(comprehensive_mapping)) {
                                for (i in seq_len(nrow(values$nodes))) {
                                    node <- values$nodes[i, ]
                                    if (!is.null(node$hgnc_symbol) && !is.na(node$hgnc_symbol) && node$hgnc_symbol != "") {
                                        # Look up UniProt for this node
                                        node_mapping <- comprehensive_mapping[
                                            !is.na(comprehensive_mapping$hgnc_symbol) & 
                                            comprehensive_mapping$hgnc_symbol == node$hgnc_symbol & 
                                            !is.na(comprehensive_mapping$uniprotswissprot) & 
                                            comprehensive_mapping$uniprotswissprot != "", 
                                        ]
                                        if (nrow(node_mapping) > 0) {
                                            pathway_genes_with_uniprot <- c(pathway_genes_with_uniprot, node$hgnc_symbol)
                                        }
                                    }
                                }
                            }
                            
                            # Remove the currently selected gene from the list
                            pathway_genes_with_uniprot <- setdiff(pathway_genes_with_uniprot, hgnc_symbol)
                            
                            if (length(pathway_genes_with_uniprot) > 0) {
                                uniprot_content_parts <- c(uniprot_content_parts, 
                                    paste0("<strong>EMERALD Protein Alignment:</strong><br>",
                                           "Select a second protein from this pathway: <br>",
                                           "<select id='second_gene_select' onchange='selectSecondGene(this.value)' style='width: 100%; margin: 5px 0;'>",
                                           "<option value=''>Choose a protein...</option>",
                                           paste0("<option value='", pathway_genes_with_uniprot, "'>", pathway_genes_with_uniprot, "</option>", collapse = ""),
                                           "</select><br>",
                                           "<div id='emerald_link_container' style='margin-top: 5px;'></div>"))
                            }
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
                    
                    # Create comprehensive HTML content with JavaScript for EMERALD alignment
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
                            paste0("or '", hgnc_symbol, "' in gene databases</em>") else "</em>",
                        
                        # Add JavaScript for EMERALD alignment functionality
                        "<script>
                        function selectSecondGene(secondGeneSymbol) {
                            if (!secondGeneSymbol) {
                                document.getElementById('emerald_link_container').innerHTML = '';
                                return;
                            }
                            
                            // Send the selected gene to Shiny
                            Shiny.setInputValue('second_gene_selected', secondGeneSymbol, {priority: 'event'});
                        }
                        </script>"
                    )
                    
                    # Return as HTML
                    HTML(html_content)
                    
                } else if (node_type == "compound") {
                    # ======= COMPOUND NODE HANDLING =======
                    
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
                                
                                source_node <- values$nodes[values$nodes$id == edge$from, ]
                                source_name <- edge$from
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
                                
                                target_node <- values$nodes[values$nodes$id == edge$to, ]
                                target_name <- edge$to
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
                    
                    # Create compound HTML content
                    html_content <- paste0(
                        "<strong>=== COMPOUND INFORMATION ===</strong><br>",
                        "<strong>Node Type:</strong> ", node_type, "<br>",
                        "<strong>KEGG Compound ID:</strong> <code>", kegg_id, "</code><br>",
                        "<strong>Original KEGG Entry:</strong> <code>", gene_name, "</code><br>",
                        "<strong>Display Label:</strong> ", display_label, "<br>",
                        relationship_info,
                        "<br><strong>=== PATHWAY CONTEXT ===</strong><br>",
                        "In KEGG, this compound is referenced as: <code>", kegg_id, "</code><br>",
                        "<br>",
                        "<em>Tip: Search for '", kegg_id, "' in KEGG COMPOUND database</em>"
                    )
                    
                    # Return as HTML
                    HTML(html_content)
                    
                } else if (node_type == "map") {
                    # ======= MAP NODE HANDLING =======
                    
                    # Extract pathway ID from the kegg_name field
                    pathway_id <- ""
                    if (!is.null(gene_name) && gene_name != "") {
                        # Extract pathway ID (e.g., "path:hsa05200" -> "hsa05200")
                        pathway_match <- regmatches(gene_name, regexpr("hsa\\d+", gene_name))
                        if (length(pathway_match) > 0) {
                            pathway_id <- pathway_match[1]
                        }
                    }
                    
                    # Build relationship info
                    relationship_info <- ""
                    if (nrow(incoming_edges) > 0 || nrow(outgoing_edges) > 0) {
                        relationship_info <- paste0(
                            "<br><strong>=== RELATIONSHIPS ===</strong><br>",
                            "Incoming connections: ", nrow(incoming_edges), "<br>",
                            "Outgoing connections: ", nrow(outgoing_edges), "<br>"
                        )
                    }
                    
                    # Create pathway link
                    pathway_link <- ""
                    if (pathway_id != "") {
                        pathway_link <- paste0(
                            "<br><strong>=== PATHWAY NAVIGATION ===</strong><br>",
                            "<strong>KEGG Pathway Link:</strong> ",
                            "<a href='https://www.kegg.jp/kegg-bin/show_pathway?", pathway_id, "' target='_blank' style='color: #007bff; text-decoration: underline;'>",
                            "View in KEGG ↗</a><br>",
                            "<strong>Load in Visualizer:</strong> ",
                            "<button id='load_pathway_btn_", node_id, "' onclick='loadPathwayFromNode(\"", pathway_id, "\")' ",
                            "style='background-color: #007bff; color: white; border: none; padding: 5px 10px; border-radius: 3px; cursor: pointer; margin-top: 5px;'>",
                            "Load ", display_label, "</button><br>"
                        )
                    }
                    
                    # Create map HTML content
                    html_content <- paste0(
                        "<strong>=== PATHWAY MAP INFORMATION ===</strong><br>",
                        if (pathway_id != "") paste0("<strong>Pathway ID:</strong> <code>", pathway_id, "</code><br>") else "",
                        "<strong>Original KEGG Entry:</strong> <code>", gene_name, "</code><br>",
                        "<strong>Display Label:</strong> ", display_label, "<br>",
                        relationship_info,
                        pathway_link,
                        "<br><strong>=== PATHWAY CONTEXT ===</strong><br>",
                        "This node represents a reference to another KEGG pathway.<br>",
                        if (pathway_id != "") paste0("You can load the ", pathway_id, " pathway to explore it in detail.<br>") else "",
                        "<br>",
                        "<em>Tip: Click 'Load Pathway' button to navigate to the referenced pathway</em>",
                        
                        # Add JavaScript for pathway loading functionality
                        "<script>
                        function loadPathwayFromNode(pathwayId) {
                            // Send the pathway ID to Shiny for loading
                            Shiny.setInputValue('load_pathway_from_node', pathwayId, {priority: 'event'});
                            
                            // Show loading feedback
                            document.getElementById('load_pathway_btn_", node_id, "').innerHTML = 'Loading...';
                            document.getElementById('load_pathway_btn_", node_id, "').disabled = true;
                        }
                        </script>"
                    )
                    
                    # Return as HTML
                    HTML(html_content)
                    
                } else {
                    # ======= OTHER NODE TYPES =======
                    
                    # Build relationship info
                    relationship_info <- ""
                    if (nrow(incoming_edges) > 0 || nrow(outgoing_edges) > 0) {
                        relationship_info <- paste0(
                            "<br><strong>=== RELATIONSHIPS ===</strong><br>",
                            "Incoming connections: ", nrow(incoming_edges), "<br>",
                            "Outgoing connections: ", nrow(outgoing_edges), "<br>"
                        )
                    }
                    
                    # Create generic HTML content
                    html_content <- paste0(
                        "<strong>=== NODE INFORMATION ===</strong><br>",
                        "<strong>KEGG ID:</strong> <code>", kegg_id, "</code><br>",
                        "<strong>Original KEGG Entry:</strong> <code>", gene_name, "</code><br>",
                        "<strong>Display Label:</strong> ", display_label, "<br>",
                        relationship_info,
                        "<br><strong>=== PATHWAY CONTEXT ===</strong><br>",
                        "In KEGG, this node is referenced as: <code>", kegg_id, "</code><br>",
                        "<br>",
                        "<em>Tip: Search for '", kegg_id, "' in KEGG database</em>"
                    )
                    
                    # Return as HTML
                    HTML(html_content)
                }
                } else {
                    HTML("<em>No node selected - click on a node to see detailed information</em>")
                }
            })
        }
    })
    
    # Handle second gene selection for EMERALD alignment
    observeEvent(input$second_gene_selected, {
        req(input$second_gene_selected, selected_gene_uniprot())
        
        second_gene_symbol <- input$second_gene_selected
        first_gene_uniprot <- selected_gene_uniprot()
        
        # Get UniProt ID for the second gene
        if (!is.null(comprehensive_mapping)) {
            second_gene_mapping <- comprehensive_mapping[
                !is.na(comprehensive_mapping$hgnc_symbol) & 
                comprehensive_mapping$hgnc_symbol == second_gene_symbol & 
                !is.na(comprehensive_mapping$uniprotswissprot) & 
                comprehensive_mapping$uniprotswissprot != "", 
            ]
            
            if (nrow(second_gene_mapping) > 0) {
                second_gene_uniprot <- second_gene_mapping$uniprotswissprot[1]
                
                # Create EMERALD alignment URL
                emerald_url <- paste0(
                    "https://algbio.github.io/emerald-ui/?seqA=", first_gene_uniprot, 
                    "&seqB=", second_gene_uniprot, "&alpha=0.75&delta=8"
                )
                
                # Update the link container via JavaScript
                session$sendCustomMessage("updateEmeraldLink", list(
                    html = paste0(
                        "<strong>🔗 EMERALD Alignment:</strong> ",
                        "<a href='", emerald_url, "' target='_blank' style='color: #007bff; text-decoration: underline;'>",
                        "Align proteins ↗", 
                        "</a><br>",
                        "<small style='color: #666;'>Compare ", first_gene_uniprot, " vs ", second_gene_uniprot, "</small>"
                    )
                ))
                
                second_gene_for_emerald(second_gene_symbol)
                showNotification(paste("EMERALD alignment ready for", first_gene_uniprot, "vs", second_gene_uniprot), type = "message")
            }
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
                scrollY = '120px',  # Match container height minus padding
                scrollCollapse = TRUE,
                autoWidth = FALSE,  # Prevent auto width calculation
                columnDefs = list(
                    list(width = '30px', targets = 0, className = 'text-center'),  # Rank column narrower
                    list(width = '120px', targets = 1)  # Name column narrower
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
    
    # Reactive values for EMERALD alignment
    selected_gene_uniprot <- reactiveVal(NULL)
    second_gene_for_emerald <- reactiveVal(NULL)
    
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
    
    # Navigation to Gene Set tab from network options
    observeEvent(input$goto_geneset_from_network, {
        updateTabItems(session, "tabs", selected = "geneset")
    })
    
    # Navigation to Evolutionary Transcriptomics tab from network
    observeEvent(input$goto_evolution_from_network, {
        updateTabItems(session, "tabs", selected = "evolution")
    })
    
    # Navigation to Pathway Explorer tab from expression card
    observeEvent(input$goto_explorer_from_expression, {
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
            geom_point(aes(size = Count, colour = p.adjust), alpha = 0.8) +
            scale_colour_gradient(low = "red", high = "blue", name = "P-adjust") +
            scale_size_continuous(name = "Gene Count", range = c(3, 10)) +
            labs(
                title = paste("Top", min(nrow(plot_data), input$plot_top_n), "Enriched KEGG Pathways"),
                x = "Fold Enrichment",
                y = "KEGG Pathway",
                subtitle = paste("Based on", length(uploaded_genes()), "input genes")
            ) +
            theme_minimal() +
            theme(
                axis.text.y = element_text(size = 10),
                plot.title = element_text(size = 14, face = "bold"),
                legend.position = "right"
            )
        
        ggplotly(p, tooltip = c("x", "y", "size", "colour"), height = 600) %>%
            layout(margin = list(l = 300))
    })
    
    # Update plot when parameters change
    observeEvent(input$update_plot, {
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
        req(input$enrichment_table_rows_selected)
        
        result <- enrichment_results()
        if (is.null(result) || length(input$enrichment_table_rows_selected) == 0) {
            return("No pathway selected")
        }
        
        selected_pathway <- result@result[input$enrichment_table_rows_selected, ]
        
        paste0(
            "SELECTED ENRICHED PATHWAY\n",
            "========================\n",
            "Pathway: ", selected_pathway$Description, "\n",
            "ID: ", selected_pathway$ID, "\n",
            "P-value: ", format(selected_pathway$pvalue, scientific = TRUE, digits = 3), "\n",
            "Adjusted p-value: ", format(selected_pathway$p.adjust, scientific = TRUE, digits = 3), "\n",
            "Fold enrichment: ", round(selected_pathway$FoldEnrichment, 2), "\n",
            "Genes in pathway: ", selected_pathway$Count, "/", strsplit(selected_pathway$BgRatio, "/")[[1]][2], "\n",
            "Your genes found: ", selected_pathway$GeneRatio
        )
    })
    
    # Load pathway selected from enrichment results with button
    observeEvent(input$load_enriched_pathway, {
        req(input$enrichment_table_rows_selected)
        
        result <- enrichment_results()
        if (is.null(result)) {
            showNotification("No enrichment results available", type = "warning")
            return()
        }
        
        selected_pathway <- result@result[input$enrichment_table_rows_selected, ]
        pathway_id <- selected_pathway$ID
        
        # Create pathway info for values$selected_pathway
        pathway_info <- list(
            pathway_id = pathway_id,
            pathway_name = selected_pathway$Description,
            description = paste("Enriched pathway (p =", format(selected_pathway$p.adjust, scientific = TRUE, digits = 3), ")")
        )
        
        showNotification(paste("Loading enriched pathway:", pathway_id), type = "message")
        
        # Load pathway data
        tryCatch({
            pathway_data <- parse_kegg_pathway_with_hsa(pathway_id, use_cached = TRUE)
            
            values$selected_pathway <- pathway_info
            values$pathway_graph <- pathway_data$kegg_data
            values$nodes <- pathway_data$nodes
            values$edges <- pathway_data$edges
            
            showNotification(
                paste("Loaded enriched pathway with", nrow(pathway_data$nodes), "nodes!"), 
                type = "message"
            )
            
            # Switch to network tab to view
            updateTabItems(session, "tabs", "network")
            
        }, error = function(e) {
            showNotification(paste("Error loading pathway:", e$message), type = "error")
        })
    })
    
    # Pathway header for network visualization
    output$pathway_header <- renderUI({
        if (!is.null(values$selected_pathway)) {
            pathway_id <- values$selected_pathway$pathway_id
            kegg_url <- paste0("https://www.kegg.jp/pathway/", pathway_id)
            
            HTML(paste0(
                "<div style='text-align: center;'>",
                "<h4 style='margin: 5px 0; color: white; text-shadow: 1px 1px 2px rgba(0,0,0,0.5);'>",
                values$selected_pathway$pathway_name,
                "</h4>",
                "<p style='margin: 5px 0; font-size: 14px; color: #e8f4f8;'>",
                "<strong>Pathway ID:</strong> ", pathway_id, " | ",
                "<a href='", kegg_url, "' target='_blank' style='color: #fff; text-decoration: underline;'>",
                "View on KEGG Database ↗",
                "</a>",
                "</p>",
                "</div>"
            ))
        }
    })
    
    # ================================================================
    # EVOLUTIONARY TRANSCRIPTOMICS TAB
    # ================================================================
    
    # Reactive values for evolutionary transcriptomics
    expression_data <- reactiveVal(NULL)
    phyloexpression_set <- reactiveVal(NULL)
    current_evolution_plot <- reactiveVal(NULL)
    evolution_plot_type <- reactiveVal("")
    
    # Load example expression data
    observeEvent(input$load_example_expr, {
        tryCatch({
            example_data <- load_example_expression_data()
            expression_data(example_data)
            showNotification("Example expression data loaded successfully!", type = "message")
        }, error = function(e) {
            showNotification(paste("Error loading example data:", e$message), type = "error")
        })
    })
    
    # Handle expression file upload
    observeEvent(input$expression_file, {
        req(input$expression_file)
        
        tryCatch({
            file_ext <- tools::file_ext(input$expression_file$datapath)
            
            if (file_ext %in% c("csv")) {
                data <- read.csv(input$expression_file$datapath, stringsAsFactors = FALSE)
            } else if (file_ext %in% c("tsv", "txt")) {
                data <- read.table(input$expression_file$datapath, header = TRUE, sep = "\t", 
                                 stringsAsFactors = FALSE, quote = '"')
            } else {
                stop("Unsupported file format. Please use CSV or TSV files.")
            }
            
            # Validate the data structure
            if (ncol(data) < 2) {
                stop("Expression data must have at least 2 columns (gene IDs + 1 sample)")
            }
            
            expression_data(data)
            showNotification("Expression data uploaded successfully!", type = "message")
            
        }, error = function(e) {
            showNotification(paste("Error uploading expression data:", e$message), type = "error")
        })
    })
    
    # Create PhyloExpressionSet
    observeEvent(input$create_phyloset, {
        req(expression_data())
        
        tryCatch({
            # Parse sample groups if provided
            groups <- NULL
            if (!is.null(input$sample_groups) && nchar(trimws(input$sample_groups)) > 0) {
                groups <- trimws(strsplit(input$sample_groups, ",")[[1]])
                
                # Validate groups length
                n_samples <- ncol(expression_data()) - 1
                if (length(groups) != n_samples) {
                    stop(paste("Number of groups (", length(groups), ") must match number of samples (", n_samples, ")"))
                }
            }
            
            # Create the phyloset using our function
            result <- create_bulk_phyloexpression_set(
                expression_data = expression_data(),
                gene_id_type = input$expr_gene_id_type,
                groups = groups,
                name = "User Expression Dataset"
            )
            
            phyloexpression_set(result)
            showNotification("PhyloExpressionSet created successfully!", type = "message")
            
        }, error = function(e) {
            showNotification(paste("Error creating PhyloExpressionSet:", e$message), type = "error")
        })
    })
    
    # Reactive outputs for evolutionary transcriptomics tab
    
    # Expression data loaded indicator
    output$expression_loaded <- reactive({
        !is.null(expression_data())
    })
    outputOptions(output, "expression_loaded", suspendWhenHidden = FALSE)
    
    # PhyloExpressionSet created indicator
    output$phyloset_created <- reactive({
        !is.null(phyloexpression_set())
    })
    outputOptions(output, "phyloset_created", suspendWhenHidden = FALSE)
    
    # Pathway expression plot for integrated visualization
    output$pathway_expression_plot <- renderPlot({
        req(phyloexpression_set(), values$nodes)
        
        tryCatch({
            phyloset <- phyloexpression_set()$phyloset
            coloring_mode <- input$node_coloring
            selected_node_id <- input$pathway_network_selected
            
            # Get pathway genes (use hgnc_symbol for matching with expression data)
            pathway_genes <- NULL
            if (!is.null(values$nodes) && nrow(values$nodes) > 0) {
                if ("hgnc_symbol" %in% colnames(values$nodes)) {
                    pathway_genes <- values$nodes$hgnc_symbol
                    pathway_genes <- pathway_genes[!is.na(pathway_genes) & pathway_genes != ""]
                }
            }
            
            if (is.null(pathway_genes) || length(pathway_genes) == 0) {
                return(ggplot() + 
                    geom_text(aes(x = 0.5, y = 0.5), 
                             label = "No genes found in current pathway", 
                             size = 5, color = "gray50") +
                    xlim(0, 1) + ylim(0, 1) + theme_void())
            }
            
            # Filter pathway genes to those present in phyloset
            phyloset_genes <- phyloset@gene_ids
            valid_pathway_genes <- intersect(pathway_genes, phyloset_genes)
            
            if (length(valid_pathway_genes) == 0) {
                return(ggplot() + 
                    geom_text(aes(x = 0.5, y = 0.5), 
                             label = "No pathway genes found in expression dataset", 
                             size = 5, color = "gray50") +
                    xlim(0, 1) + ylim(0, 1) + theme_void())
            }
            
            # Determine what genes and colors to use based on selection and coloring mode
            if (!is.null(selected_node_id) && length(selected_node_id) > 0 && selected_node_id != "") {
                # NODE SELECTED: Check node type first
                selected_node <- values$nodes[values$nodes$id == selected_node_id, ]
                
                if (nrow(selected_node) > 0) {
                    # Check if this is a gene node (only gene nodes show interactions)
                    node_type <- selected_node$type[1]
                    
                    if (!is.na(node_type) && node_type == "gene") {
                        # GENE NODE: Show selected gene + connections
                        selected_gene <- selected_node$hgnc_symbol[1]
                        
                        if (!is.na(selected_gene) && selected_gene != "") {
                            # Use the get_interacting_genes function to get proper colors
                            interaction_data <- get_interacting_genes(
                                selected_gene = selected_gene,
                                pathway_edges = values$edges,
                                pathway_nodes = values$nodes,
                                gene_id_type = "symbol"
                            )
                            
                            # Get all genes to plot (selected + interacting)
                            all_genes_to_plot <- c(selected_gene, interaction_data$interacting_genes)
                            
                            # Filter to genes present in phyloset
                            all_genes_to_plot <- intersect(all_genes_to_plot, phyloset_genes)
                            
                            if (length(all_genes_to_plot) > 0) {
                                # Use the interaction colors from get_interacting_genes
                                gene_colors <- interaction_data$interaction_colors
                                
                                # Filter colors to only genes present in phyloset
                                gene_colors <- gene_colors[names(gene_colors) %in% all_genes_to_plot]
                                
                                # Additional validation
                                if (any(is.na(names(gene_colors))) || any(names(gene_colors) == "")) {
                                    valid_color_indices <- !is.na(names(gene_colors)) & names(gene_colors) != ""
                                    gene_colors <- gene_colors[valid_color_indices]
                                }
                                
                                # Create the plot
                                p <- create_gene_profiles_plot(phyloset,
                                                             selected_genes = all_genes_to_plot,
                                                             interaction_colors = gene_colors,
                                                             title = paste0("Expression Profiles: ", selected_gene, " + ", 
                                                                           length(interaction_data$interacting_genes), " connections"))
                                return(p)
                            }
                        }
                    }
                    # For compound or map nodes, fall through to show all pathway genes (same as no selection)
                }
            }
            
            # NO NODE SELECTED OR NON-GENE NODE SELECTED: Show profiles of all pathway genes
                
            if (coloring_mode == "kegg_default") {
                # KEGG Default: Show gene heatmap instead of profiles
                p <- create_gene_heatmap_plot(phyloset,
                                            selected_genes = valid_pathway_genes,
                                            title = paste0("Gene Expression Heatmap (", length(valid_pathway_genes), " pathway genes)"))
                return(p)
                
            } else if (coloring_mode == "highlight_genes") {
                # Highlight uploaded genes: red for uploaded, green for others
                uploaded_genes <- internal_genes()  # Get uploaded Entrez IDs
                uploaded_symbols <- entrez_to_symbols(uploaded_genes, comprehensive_mapping)
                
                gene_colors <- rep("green3", length(valid_pathway_genes))
                names(gene_colors) <- valid_pathway_genes
                
                # Color uploaded genes red
                uploaded_in_pathway <- intersect(uploaded_symbols, valid_pathway_genes)
                if (length(uploaded_in_pathway) > 0) {
                    gene_colors[uploaded_in_pathway] <- "red"
                }
                
                p <- create_gene_profiles_plot(phyloset,
                                                selected_genes = valid_pathway_genes,
                                                interaction_colors = gene_colors,
                                                title = paste0("Expression Profiles (", length(uploaded_in_pathway), " highlighted, ", 
                                                            length(valid_pathway_genes) - length(uploaded_in_pathway), " others)"))
                return(p)
                
            } else {
                # Phylostratum: Show gene heatmap colored by phylostratum
                p <- create_gene_heatmap_plot(phyloset,
                                            selected_genes = valid_pathway_genes,
                                            title = paste0("Gene Expression Heatmap (", length(valid_pathway_genes), " pathway genes)"))
                return(p)
            }
            
            # Fallback
            return(ggplot() + 
                geom_text(aes(x = 0.5, y = 0.5), 
                         label = "Unable to create expression profiles", 
                         size = 5, color = "gray50") +
                xlim(0, 1) + ylim(0, 1) + theme_void())
                
        }, error = function(e) {
            return(ggplot() + 
                geom_text(aes(x = 0.5, y = 0.5), 
                         label = paste0("Error: ", e$message), 
                         size = 4, color = "red") +
                xlim(0, 1) + ylim(0, 1) + theme_void())
        })
    })
    
    # Evolution plot ready indicator
    output$evolution_plot_ready <- reactive({
        !is.null(current_evolution_plot())
    })
    outputOptions(output, "evolution_plot_ready", suspendWhenHidden = FALSE)
    
    # Expression data summary
    output$expression_summary <- renderText({
        req(expression_data())
        data <- expression_data()
        
        paste(
            paste("Genes:", nrow(data)),
            paste("Samples:", ncol(data) - 1),
            paste("Gene ID type:", input$expr_gene_id_type),
            sep = "\n"
        )
    })
    
    # Expression matrix preview (first 50 rows sorted by mean expression)
    output$expression_preview <- DT::renderDataTable({
        req(expression_data())
        data <- expression_data()
        
        # Calculate mean expression for each gene (excluding the gene ID column)
        expression_cols <- names(data)[-1]  # All columns except the first (gene ID)
        data$mean_expression <- rowMeans(data[expression_cols], na.rm = TRUE)
        
        # Sort by mean expression (descending) and show top 50
        data_sorted <- data[order(data$mean_expression, decreasing = TRUE), ]
        preview_data <- head(data_sorted, 50)
        
        # Remove the temporary mean_expression column for display
        preview_data$mean_expression <- NULL
        
        DT::datatable(
            preview_data,
            options = list(
                scrollX = TRUE,
                pageLength = 50,
                searching = FALSE,
                lengthChange = TRUE,
                info = TRUE,
                paging = TRUE,
                ordering = TRUE,
                dom = 'ltp',
                scrollY = "300px",
                autoWidth = TRUE,
                lengthMenu = list(c(10, 25, 50, -1), c("10", "25", "50", "All")),
                columnDefs = list(
                    list(className = 'dt-center', targets = "_all")
                ),
                initComplete = DT::JS(
                    "function(settings, json) {",
                    "$(this.api().table().header()).css({'font-size': '11px', 'padding': '4px'});",
                    "$(this.api().table().body()).css({'font-size': '10px', 'padding': '2px'});",
                    "$(this.api().table().container()).css({'margin': '0px', 'padding': '0px'});",
                    "}"
                )
            ),
            rownames = FALSE,
            class = 'compact stripe hover',
            caption = "Top 50 genes sorted by mean expression (highest first)"
        ) %>%
        DT::formatStyle(columns = 1:ncol(preview_data), fontSize = '10px') %>%
        DT::formatRound(columns = 2:ncol(preview_data), digits = 3)
    }, server = FALSE)
    
    # Toggle button text for expression preview
    observeEvent(input$toggle_expression_preview, {
        if (input$toggle_expression_preview %% 2 == 1) {
            updateActionButton(session, "toggle_expression_preview", 
                             label = "Hide Expression Matrix Preview",
                             icon = icon("eye-slash"))
        } else {
            updateActionButton(session, "toggle_expression_preview", 
                             label = "Show Expression Matrix Preview (Top 50 genes)",
                             icon = icon("eye"))
        }
    })
    
    # Automatic plot generation when phyloset is created
    
    # TAI Signature plot (automatic)
    output$tai_signature_plot <- renderPlot({
        req(phyloexpression_set())
        tryCatch({
            create_tai_signature_plot(phyloexpression_set()$phyloset, 
                                    title = NULL)
        }, error = function(e) {
            ggplot() + 
                geom_text(aes(x = 0.5, y = 0.5), label = paste("Error:", e$message), size = 5) +
                xlim(0, 1) + ylim(0, 1) + theme_void()
        })
    })
    
    # Distribution of phylostrata plot (automatic)
    output$phylostrata_distribution_plot <- renderPlot({
        req(phyloexpression_set())
        tryCatch({
            create_distribution_strata_plot(phyloexpression_set()$phyloset,
                                          title = NULL)
        }, error = function(e) {
            ggplot() + 
                geom_text(aes(x = 0.5, y = 0.5), label = paste("Error:", e$message), size = 5) +
                xlim(0, 1) + ylim(0, 1) + theme_void()
        })
    })
    
    # Sample space plot (automatic)
    output$sample_space_plot <- renderPlot({
        req(phyloexpression_set())
        tryCatch({
            create_sample_space_plot(phyloexpression_set()$phyloset,
                                   title = NULL)
        }, error = function(e) {
            ggplot() + 
                geom_text(aes(x = 0.5, y = 0.5), label = paste("Error:", e$message), size = 5) +
                xlim(0, 1) + ylim(0, 1) + theme_void()
        })
    })
    
    # Gene heatmap plot (automatic)
    output$gene_heatmap_plot <- renderPlot({
        req(phyloexpression_set())
        tryCatch({
            create_gene_heatmap_plot(phyloexpression_set()$phyloset,
                                   title = NULL)
        }, error = function(e) {
            ggplot() + 
                geom_text(aes(x = 0.5, y = 0.5), label = paste("Error:", e$message), size = 5) +
                xlim(0, 1) + ylim(0, 1) + theme_void()
        })
    })

    # PhyloExpressionSet mapping information
    output$phyloset_mapping_info <- renderText({
        req(phyloexpression_set())
        
        data <- expression_data()
        stats <- phyloexpression_set()$mapping_stats
        
        lines <- c(
            paste("✓ PhyloExpressionSet Successfully Created"),
            paste(""),
            paste("Dataset Information:"),
            paste("  • Total genes:", stats$n_total),
            paste("  • Genes mapped to phylostrata:", stats$n_mapped, paste0("(", stats$mapping_rate, "%)")),
            paste("  • Total samples:", ncol(data) - 1),
            paste("  • Sample names:", paste(names(data)[-1], collapse = ", "))
        )
        
        paste(lines, collapse = "\n")
    })
    
    # Evolution plot title
    output$evolution_plot_title <- renderText({
        req(current_evolution_plot())
        paste0("<h3 style='margin-bottom: 20px; color: #2c3e50;'>", evolution_plot_type(), "</h3>")
    })
    
    # Evolution plot
    output$evolution_plot <- renderPlot({
        req(current_evolution_plot())
        current_evolution_plot()
    })
    
    # Evolution strata legend table
    output$evolution_strata_legend <- DT::renderDataTable({
        strata_legend <- load_strata_legend()
        if (!is.null(strata_legend)) {
            # Generate colors for phylostrata like in network tab
            max_stratum <- max(strata_legend$Rank)
            all_colors <- PS_colours(max_stratum)
            
            # Create legend data frame with colors
            legend_df <- data.frame(
                Rank = strata_legend$Rank,
                Name = strata_legend$Name,
                Color = all_colors[strata_legend$Rank],
                stringsAsFactors = FALSE
            )
            
            # Create colored table like in network tab
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
        } else {
            DT::datatable(data.frame(Message = "Phylostratum legend data not available"))
        }
    })
    
    # Gene selection status - DISABLED for streamlined UI
    # output$evolution_selection_status <- renderText({
    #     # This was used for the old manual gene selection interface
    # })
    
    # Dataset information - DISABLED for streamlined UI  
    # output$evolution_dataset_info <- renderText({
    #     # This was used for the old dataset info card
    # })
    
    # Download evolution plot
    output$download_evolution_plot <- downloadHandler(
        filename = function() {
            plot_type <- gsub("[^A-Za-z0-9_-]", "_", evolution_plot_type())
            paste0("evolution_", plot_type, "_", Sys.Date(), ".png")
        },
        content = function(file) {
            req(current_evolution_plot())
            ggsave(file, plot = current_evolution_plot(), 
                   width = 12, height = 8, dpi = 300, bg = "white")
        }
    )
}
