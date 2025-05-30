# Main server logic -----------------------------------------------------------

# Load dependencies
source("global.R")
source("heatmap.R")
source("modules/utilities.R")
source("modules/data_processing.R")
source("modules/heatmap_module.R")
source("modules/rendering_functions.R")

# Initialize reactive values
data_cache <- reactiveVal(list(name = NULL, df = NULL))
cluster_map <- reactiveVal(NULL) 
selected_file <- reactiveVal(NULL)

server <- function(input, output, session) {
  data <- reactiveVal()
  heatmap_obj <- reactiveVal()
  display_labels <- reactiveVal()
  output$sub_heatmaply <- renderPlotly(NULL)
  
  # Reset inputs when new file is selected
  observeEvent(input$file, {
    updateSelectInput(session, "highlight_category", choices = NULL, selected = "Top5")
    updateRadioButtons(session, "matrix_type", selected = "Full Matrix")
    updateRadioButtons(session, "filter_mode", selected = "None")
    updateNumericInput(session, "matrix_size", value = 50)
    updateNumericInput(session, "top_n_genes", value = 100)
    updateSliderInput(session, "inside_range", value = c(-0.5, 0.5))
    updateSliderInput(session, "outside_range", value = c(-0.5, 0.5))
  })
  
  # Process data when submit button is clicked
  observeEvent(input$submit, {
    selected_file(input$file)
    req(selected_file())
    showNotification("Loading data.", duration = NULL, id = "load_msg")
    
    df <- if (!is.null(data_cache()$name) && identical(data_cache()$name, input$file)) {
      data_cache()$df
    } else {
      df <- readRDS(url(file_paths[[input$file]]))
      data_cache(list(name = input$file, df = df))
      df
    }
    
    df <- process_matrix(df, input)
    if (is.null(df)) {
      showNotification("No genes left after filtering.", type = "error")
      removeNotification("load_msg")
      return()
    }
    
    result <- prepare_heatmap_data(input, df, session)
    
    data(list(
      mat = result$mat,
      meta = result$meta,
      hc_rows = result$hc_rows,
      row_categories = result$row_categories,
      row_categories_full = result$row_categories_full,
      top_categories = result$top_categories
    ))
    
    cluster_map(result$clusters)
    display_labels(result$display_labels)
    
    heatmap_obj(result$heatmap)
    
    removeNotification("load_msg")
  })
  
  # Handle heatmap interactions
  observeEvent(input$submit, {
    req(heatmap_obj())
    
    ht <- draw(heatmap_obj(), test = FALSE)
    makeInteractiveComplexHeatmap(
      input, output, session, ht,
      heatmap_id = "heatmap_output",
      
      hover_action = function(df, input, output, session) {
        d <- data()
        mat <- d$mat
        row_categories <- d$row_categories
        labels <- display_labels()
        
        get_display <- function(gene_id) {
          if (!is.null(labels) && gene_id %in% names(labels)) {
            return(labels[[gene_id]])
          } else {
            return(gene_id)
          }
        }
        
        if (!is.null(df)) {
          value <- round(mat[df$row_index, df$column_index], 3)
          row_id <- df$row_label
          col_id <- df$column_label
          row_display <- get_display(row_id)
          col_display <- get_display(col_id)
          category <- if (!is.null(row_categories)) row_categories[[row_id]] else "—"
          
          msg <- paste0(
            "Row: ", row_display, "\n",
            "Column: ", col_display, "\n",
            "Value: ", value, "\n",
            "Category: ", category
          )
          
          shinyjs::runjs(sprintf(
            "document.getElementById('heatmap_tooltip').innerText = `%s`; document.getElementById('heatmap_tooltip').style.display = 'block';",
            gsub("`", "\\`", msg)
          ))
        } else {
          shinyjs::runjs("document.getElementById('heatmap_tooltip').style.display = 'none';")
        }
      },
      
      brush_action = function(df, input, output, session) {
        if (isolate(input$subheatmap_mode) != "Zoom") return()
        
        d <- data()
        mat <- d$mat
        row_categories <- d$row_categories
        top_categories <- d$top_categories
        labels <- display_labels()
        
        selected_ids <- rownames(mat)[unique(unlist(df$row_index))]
        selected_cols <- colnames(mat)[unique(unlist(df$column_index))]
        
        if (length(selected_ids) > 350) {
          showNotification("Too many genes selected (>350). Sub-heatmap, table, and dotplot will not be displayed.", type = "warning", duration = 5000, id = "sub1")
          output$sub_heatmaply <- renderPlotly(NULL)
          output$selected_table <- DT::renderDT(NULL)
          output$dotplot <- renderPlot({ ggplot() + theme_void() })
          return()
        }
        removeNotification("sub1")
        removeNotification("sub2")
        output$sub_heatmaply <- renderPlotly({
          render_sub_heatmap(
            mat = mat,
            row_ids = selected_ids,
            col_ids = selected_cols,
            row_categories = row_categories,
            top_categories = top_categories,
            labels = labels,
            selected_category = input$highlight_category
          )
        })
        
        render_selected_table(output, mat, labels, selected_ids, selected_cols)
        render_dotplot(output, selected_file(), selected_ids, labels)
      }
    )
  })
  
  # Handle cluster selection
  observeEvent(input$selected_cluster, {
    req(input$subheatmap_mode == "Cluster", cluster_map(), data())
    
    if (is.null(input$selected_cluster) || input$selected_cluster == "") return()
    
    clusters <- cluster_map()
    selected_cluster <- as.numeric(input$selected_cluster)
    
    if (!selected_cluster %in% clusters) return()
    
    gene_ids <- names(clusters[clusters == selected_cluster])
    
    if (length(gene_ids) > 350) {
      showNotification("Selected cluster has too many genes (>350). Sub-heatmap, table, and dotplot will not be displayed.", type = "warning", duration = 5000, id= "sub2")
      output$sub_heatmaply <- renderPlotly(NULL)
      output$selected_table <- DT::renderDT(NULL)
      output$dotplot <- renderPlot({ ggplot() + theme_void() })  # Empty plot
      return()
    }
    removeNotification("sub2")
    removeNotification("sub1")
    d <- data()
    mat <- d$mat
    labels <- display_labels()
    
    output$sub_heatmaply <- renderPlotly({
      render_sub_heatmap(
        mat = mat,
        row_ids = gene_ids,
        col_ids = gene_ids,
        row_categories = d$row_categories,
        top_categories = d$top_categories,
        labels = labels,
        selected_category = input$highlight_category
      )
    })
    
    render_selected_table(output, mat, labels, gene_ids)
    render_dotplot(output, input$file, gene_ids, labels)
  })
}