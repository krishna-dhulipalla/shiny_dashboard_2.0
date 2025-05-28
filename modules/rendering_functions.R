# Output rendering functions ---------------------------------------------------

#' Render sub-heatmap plot
render_sub_heatmap <- function(mat, gene_ids, row_categories, top_categories, labels, selected_category) {
  if (length(gene_ids) < 2) return(NULL)
  
  sub_mat <- mat[gene_ids, gene_ids, drop = FALSE]
  gene_labels <- labels[gene_ids]

  # Prepare category annotations
  if (!is.null(row_categories)) {
    selected_row_categories <- row_categories[gene_ids]
    
    if (!is.null(selected_category) && selected_category != "Top5") {
      selected_row_categories[is.na(selected_row_categories) | selected_row_categories != selected_category] <- "Other"
    } else {
      selected_row_categories[selected_row_categories == ""] <- NA
      selected_row_categories[is.na(selected_row_categories)] <- "Other"
    }
  } else {
    selected_row_categories <- rep("Uncategorized", length(gene_ids))
  }

  names(selected_row_categories) <- gene_labels
  row_ann <- data.frame(Category = factor(selected_row_categories), stringsAsFactors = FALSE)
  rownames(row_ann) <- gene_labels

  # Define colors
  if (!is.null(selected_category) && selected_category != "Top5") {
    category_colors_final <- setNames(c("#E41A1C", "#999999"), c(selected_category, "Other"))
    row_ann$Category <- as.character(row_ann$Category)
    row_ann$Category[is.na(row_ann$Category) | row_ann$Category != selected_category] <- "Other"
    row_ann$Category <- factor(row_ann$Category, levels = names(category_colors_final))
  } else {
    all_levels <- setdiff(top_categories, "Other")
    all_levels <- all_levels[!is.na(all_levels)]
    all_levels <- c(all_levels, "Other")
    row_ann$Category <- factor(row_ann$Category, levels = all_levels)
    category_colors_final <- setNames(my_colors[seq_along(all_levels)], all_levels)
  }

  # Build heatmap
  add_side_annotation <- !(is.null(row_categories) || all(is.na(row_ann$Category)))
  
  p <- heatmaply(
    sub_mat,
    colors = colorRampPalette(c("blue", "white", "red"))(256),
    limits = c(-1, 1),
    dendrogram = "none",
    Rowv = NULL,
    Colv = NULL,
    showticklabels = c(TRUE, TRUE),
    labRow = gene_labels,
    labCol = gene_labels,
    row_side_colors = if (add_side_annotation) row_ann else NULL,
    row_side_palette = if (add_side_annotation) category_colors_final else NULL,
    hide_colorbar = TRUE,
    fontsize_row = 8,
    fontsize_col = 8,
    subplot_widths = if (add_side_annotation) c(0.98, 0.02) else 1,
    margins = c(5, 5, 5, 5),
    plot_method = "plotly"
  )

  for (i in seq_along(p$x$data)) {
    p$x$data[[i]]$showlegend <- FALSE
    p$x$data[[i]]$showscale <- FALSE
  }

  return(p %>% layout(margin = list(t = 5, b = 5, l = 5, r = 5), xaxis = list(tickangle = 270)))
}


# Render Selcted table
render_selected_table <- function(output, mat, labels, selected_ids, selected_cols = NULL) {
  selected_cols <- selected_cols %||% selected_ids
  selected_data <- mat[selected_ids, selected_cols, drop = FALSE]
  
  output$selected_table <- DT::renderDT({
    rownames(selected_data) <- labels[selected_ids]
    colnames(selected_data) <- labels[selected_cols]
    datatable(round(selected_data, 3), options = list(pageLength = 5, scrollX = TRUE))
  })
}

# Render Dot plot
render_dotplot <- function(output, file_name, selected_ids, labels) {
  seurat_file <- seurat_paths[[file_name]]
  seurat_obj <- tryCatch({ readRDS(url(seurat_file)) }, error = function(e) NULL)
  valid_genes <- selected_ids
  if (is.null(seurat_obj)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Seurat object unavailable", size = 6) + theme_void())
  }
  
  output$dotplot <- renderPlot({
    if (!is.null(seurat_obj)) {
      seurat_features <- rownames(seurat_obj)
      valid_genes <- intersect(valid_genes, seurat_features)
      
      if (length(valid_genes) > 0) {
        dp <- DotPlot(seurat_obj, features = valid_genes, dot.scale = 8)
        dp$data <- subset(dp$data, pct.exp > 0)
        
        dp +
          scale_color_gradient(low = "blue", high = "red") +
          scale_size(range = c(2, 8)) +
          scale_x_discrete(labels = labels[valid_genes]) +
          theme_minimal(base_size = 14) +
          theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
          )
      } else {
        ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No valid genes for Seurat DotPlot", size = 6) + theme_void()
      }
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No Seurat object found for this dataset", size = 6) + theme_void()
    }
  })
}