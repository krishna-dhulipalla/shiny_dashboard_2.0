color_mapping <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
my_colors <- c(
  "#E41A1C",  # Red
  "#377EB8",  # Blue
  "#4DAF4A",  # Green
  "#984EA3",  # Purple
  "gray"   # Orange
)

create_heatmap <- function(mat, row_categories = NULL, highlight_category = NULL, num_clusters = NULL, hc_rows = NULL) {
  if (is.null(hc_rows)) {
    d <- dist(mat)
    hc_rows <- fastcluster::hclust(d)
  }
  
  annotations_list <- list()
  color_list <- list()
  split_factor <- NULL  # Initialize split factor
  
  # 1. Category annotation
  if (!is.null(row_categories)) {
    if (!is.null(highlight_category) && highlight_category != "") {
      if (highlight_category == "Other") {
        cat_highlighted <- ifelse(row_categories == "Other", "Other", "Non-Other")
        cat_factor <- factor(cat_highlighted, levels = c("Other", "Non-Other"))
        annotations_list$Category <- cat_factor
        color_list$Category <- setNames(c("gray", "black"), c("Other", "Non-Other"))
      } else {
        cat_highlighted <- ifelse(row_categories == highlight_category, highlight_category, "Other")
        cat_factor <- factor(cat_highlighted, levels = c(highlight_category, "Other"))
        annotations_list$Category <- cat_factor
        color_list$Category <- setNames(c("red", "gray"), c(highlight_category, "Other"))
      }
    } else {
      cat_counts <- sort(table(row_categories[row_categories != "Other"]), decreasing = TRUE)
      top5_levels <- names(cat_counts)[1:min(5, length(cat_counts))]
      category_levels <- c(top5_levels, "Other")
      row_categories <- ifelse(row_categories %in% top5_levels, row_categories, "Other")
      row_categories <- factor(row_categories, levels = category_levels)
      annotations_list$Category <- row_categories
      color_list$Category <- setNames(
        c(my_colors[seq_along(top5_levels)], "gray"), 
        c(top5_levels, "Other")
      )
    }
  }
  
  # 2. Cluster handling and heatmap splitting
  if (!is.null(num_clusters)) {
    clusters <- cutree(hc_rows, k = num_clusters)
    row_order <- hc_rows$order  # Dendrogram order
    clusters_ordered <- clusters[row_order]
    
    cluster_levels <- unique(clusters_ordered)
    cluster_factor <- factor(clusters, levels = cluster_levels)
    
    annotations_list$Cluster <- cluster_factor
    split_factor <- cluster_factor
    
    if (length(cluster_levels) > length(my_colors)) {
      set.seed(123)
      cluster_colors <- rainbow(length(cluster_levels))
    } else {
      cluster_colors <- my_colors[seq_along(cluster_levels)]
    }
    color_list$Cluster <- setNames(cluster_colors, cluster_levels)
  }
  
  # 3. Create annotation
  anno_obj <- NULL
  if (length(annotations_list) > 0) {
    anno_obj <- ComplexHeatmap::rowAnnotation(
      df = as.data.frame(annotations_list),
      col = color_list,
      show_annotation_name = TRUE
    )
  }
  
  # 4. Generate heatmap
  Heatmap(
    mat,
    name = "Correlation",
    col = color_mapping,
    right_annotation = anno_obj,
    cluster_rows = as.dendrogram(hc_rows),
    cluster_columns = as.dendrogram(hc_rows),
    show_row_names = FALSE,
    show_column_names = FALSE,
    use_raster = TRUE
  )
}