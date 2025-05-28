# Matrix processing and filtering functions -------------------------------------

#' Process matrix based on user input
#' 
#' @param df Input data frame
#' @param input Shiny input object
#' @return Processed data frame
process_matrix <- function(df, input) {
  if (input$matrix_type == "Reduced Matrix") {
    df <- reduce_matrix(df, input$matrix_size)
  } else if (input$matrix_type == "Top N Genes") {
    df <- filter_top_genes(df, input$top_n_genes)
  }
  if (input$filter_mode != "None") {
    df <- filter_matrix_by_mode(df, input)
  }
  return(df)
}

#' Reduce matrix size by sampling
reduce_matrix <- function(df, target_size) {
  if (nrow(df) <= target_size) return(df)
  step <- ceiling(nrow(df) / target_size)
  idx <- seq(1, nrow(df), by = step)
  df_reduced <- df[idx, , drop = FALSE]
  numeric_cols <- !(colnames(df) %in% c("GeneID", "category", "GeneName"))
  mat_reduced <- df_reduced[, numeric_cols, drop = FALSE][, idx, drop = FALSE]
  cbind(df_reduced[, intersect(c("GeneID", "category", "GeneName"), colnames(df)), drop = FALSE], mat_reduced)
}

#' Filter matrix based on selected mode
filter_matrix_by_mode <- function(df, input) {
  numeric_cols <- !(colnames(df) %in% c("GeneID", "category", "GeneName"))
  mat <- as.matrix(df[, numeric_cols, drop = FALSE])
  filtered_mat <- mat
  range <- if (input$filter_mode == "Inside Range") input$inside_range else input$outside_range
  if (input$filter_mode == "Inside Range") {
    filtered_mat[mat < range[1] | mat > range[2]] <- 0
  } else if (input$filter_mode == "Outside Range") {
    filtered_mat[mat >= range[1] & mat <= range[2]] <- 0
  }
  keep <- which(rowSums(abs(filtered_mat)) > 0)
  if (length(keep) < 2) return(NULL)
  df_filtered <- df[keep, , drop = FALSE]
  filtered_mat <- filtered_mat[keep, keep, drop = FALSE]
  cbind(df_filtered[, intersect(c("GeneID", "category", "GeneName"), colnames(df)), drop = FALSE], filtered_mat)
}

#' Filter top N genes by score
filter_top_genes <- function(df, top_n) {
  numeric_cols <- !(colnames(df) %in% c("GeneID", "category", "GeneName"))
  mat <- as.matrix(df[, numeric_cols, drop = FALSE])
  gene_scores <- rowSums(abs(mat))
  top_indices <- order(gene_scores, decreasing = TRUE)[1:min(top_n, nrow(mat))]
  df_top <- df[top_indices, , drop = FALSE]
  mat_top <- mat[top_indices, top_indices, drop = FALSE]
  cbind(df_top[, intersect(c("GeneID", "category", "GeneName"), colnames(df)), drop = FALSE], mat_top)
}