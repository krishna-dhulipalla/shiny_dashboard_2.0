# Heatmap preparation and rendering logic --------------------------------------

#' Prepare heatmap data and metadata
prepare_heatmap_data <- function(input, df, session) {
  meta_cols <- intersect(c("GeneID", "category", "GeneName"), colnames(df))
  meta <- df[, meta_cols, drop = FALSE]
  mat <- as.matrix(df[, !(colnames(df) %in% meta_cols), drop = FALSE])
  
  rownames(mat) <- meta$GeneID
  colnames(mat) <- meta$GeneID
  display_vec <- get_display_labels(meta, input$label_type == "GeneName")
  
  d <- dist(mat)
  hc_rows <- fastcluster::hclust(d)
  clusters <- cutree(hc_rows, k = input$num_clusters)
  updateSelectInput(session, "selected_cluster", 
                    choices = c("None", sort(unique(clusters))))
  
  row_categories_full <- NULL
  row_categories_grouped <- NULL
  top_categories <- character(0)
  
  if ("category" %in% colnames(df)) {
    row_categories_full <- normalize_categories(setNames(df$category, df$GeneID))
    cat_counts <- sort(table(row_categories_full), decreasing = TRUE)
    top_categories <- names(cat_counts)[1:5]
    
    row_categories_grouped <- row_categories_full
    row_categories_grouped[!row_categories_grouped %in% top_categories] <- "Other"
    if (input$highlight_category != "Top5") {
      row_categories_grouped <- row_categories_full
    }
    
    category_choices <- c("Top 5 Categories" = "Top5",
                          setNames(names(cat_counts), paste0(names(cat_counts), " (", cat_counts, ")")))
    if ("Other" %in% row_categories_grouped) {
      category_choices <- c(category_choices, "Other (aggregated)" = "Other")
    }
    updateSelectInput(session, "highlight_category", choices = category_choices, selected = input$highlight_category)
  } else {
    updateSelectInput(session, "highlight_category", choices = NULL, selected = "Top5")
  }
  
  highlight_final <- if (input$highlight_category == "Top5") NULL else input$highlight_category
  
  ht <- create_heatmap(
    mat = mat,
    row_categories = row_categories_grouped,
    highlight_category = highlight_final,
    num_clusters = input$num_clusters,
    hc_rows = hc_rows
  )
  
  list(
    heatmap = ht,
    mat = mat,
    meta = meta,
    hc_rows = hc_rows,
    clusters = clusters,
    display_labels = display_vec,
    row_categories = row_categories_grouped,
    row_categories_full = row_categories_full,
    top_categories = top_categories
  )
}