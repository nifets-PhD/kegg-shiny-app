# Define Server
server <- function(input, output, session) {
    
    # Reactive values to store data
    values <- reactiveValues(
        pathways_list = NULL,
        selected_pathway = NULL,
        pathway_graph = NULL,
        nodes = NULL,
        edges = NULL,
        annotations = NULL
    )
    
    # Load KEGG pathways on app start
    observe({
        showNotification("Loading KEGG pathways...", type = "message")
        values$pathways_list <- get_kegg_pathways()
        showNotification("KEGG pathways loaded successfully!", type = "message")
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
        
        # Debug print
        cat("Selected row:\n")
        print(selected_row)
        cat("Pathway ID:", selected_row$pathway_id, "\n")
        
        values$selected_pathway <- selected_row
        
        # Load pathway graph
        tryCatch({
            pathway_data <- load_kegg_pathway(selected_row$pathway_id)
            values$pathway_graph <- pathway_data$graph
            values$nodes <- pathway_data$nodes
            values$edges <- pathway_data$edges
            
            showNotification("Pathway loaded successfully!", type = "message")
            
            # Switch to network tab
            updateTabItems(session, "tabs", "network")
            
        }, error = function(e) {
            showNotification(paste("Error loading pathway:", e$message), type = "error")
        })
    })
    
    # Display pathway information
    output$pathway_info <- renderText({
        if (!is.null(values$selected_pathway)) {
            coord_info <- ""
            if (!is.null(values$nodes) && !is.null(values$nodes$x)) {
                coord_count <- sum(!is.na(values$nodes$x))
                coord_info <- paste0("KEGG Layout Coordinates: ", coord_count, " nodes positioned\n")
            }
            
            paste(
                "=== PATHWAY INFORMATION ===", "\n",
                "Pathway ID:", values$selected_pathway$pathway_id, "\n",
                "Name:", values$selected_pathway$pathway_name, "\n",
                "Description:", values$selected_pathway$description, "\n\n",
                "=== NETWORK STATISTICS ===", "\n",
                "Total Nodes:", ifelse(!is.null(values$nodes), nrow(values$nodes), "Not loaded"), "\n",
                "Total Edges:", ifelse(!is.null(values$edges), nrow(values$edges), "Not loaded"), "\n",
                coord_info, "\n",
                "=== GENE ID INFORMATION ===", "\n",
                "• Nodes show KEGG Gene IDs prominently", "\n",
                "• HGNC symbols shown in parentheses when available", "\n",
                "• Click on nodes to see detailed ID mappings", "\n",
                "• Use 'KEGG Original Layout' for pathway-specific positioning", "\n\n",
                "Tip: KEGG IDs are more specific for pathway analysis than gene symbols!"
            )
        } else {
            "No pathway selected - choose a pathway from the Pathway Explorer tab"
        }
    })
    
    # Render network visualization
    output$pathway_network <- renderVisNetwork({
        req(values$nodes, values$edges)
        
        create_kegg_network_visualization(
            values$nodes, 
            values$edges,
            show_labels = TRUE,
            show_edges = TRUE,
            coloring_mode = ifelse(is.null(input$node_coloring), "kegg_default", input$node_coloring)
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
                        "or '", selected_node$hgnc_symbol, "' in gene databases"
                    )
                } else {
                    "No node selected - click on a node to see detailed gene information"
                }
            })
        }
    })
    
    # Handle annotation file upload
    observeEvent(input$annotation_file, {
        req(input$annotation_file)
        
        tryCatch({
            ext <- tools::file_ext(input$annotation_file$datapath)
            if (ext == "csv") {
                values$annotations <- read.csv(input$annotation_file$datapath, stringsAsFactors = FALSE)
            } else if (ext %in% c("tsv", "txt")) {
                values$annotations <- read.delim(input$annotation_file$datapath, stringsAsFactors = FALSE)
            }
            
            # Update column choices
            updateSelectInput(session, "annotation_column",
                             choices = names(values$annotations))
            
            showNotification("Annotation file loaded successfully!", type = "message")
            
        }, error = function(e) {
            showNotification(paste("Error loading file:", e$message), type = "error")
        })
    })
    
    # Preview annotations
    output$annotation_preview <- DT::renderDataTable({
        req(values$annotations)
        values$annotations
    }, options = list(pageLength = 10, scrollX = TRUE))
    
    # Apply annotations to network
    observeEvent(input$apply_annotations, {
        req(values$nodes, values$annotations, input$annotation_column)
        
        tryCatch({
            # Apply annotations to nodes
            annotated_nodes <- apply_gene_annotations(
                values$nodes, 
                values$annotations, 
                input$annotation_column,
                input$color_scheme
            )
            
            # Update network
            visNetworkProxy("pathway_network") %>%
                visUpdateNodes(nodes = annotated_nodes$nodes)
            
            # Create color legend
            output$color_legend <- renderPlotly({
                create_color_legend(annotated_nodes$color_scale, input$annotation_column)
            })
            
            showNotification("Annotations applied successfully!", type = "message")
            
        }, error = function(e) {
            showNotification(paste("Error applying annotations:", e$message), type = "error")
        })
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
}
