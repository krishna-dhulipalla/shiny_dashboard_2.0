# Utility functions for data formatting and normalization ----------------------

#' Get display labels for genes
#' 
#' @param meta DataFrame containing gene metadata
#' @param use_gene_name Logical indicating whether to use gene names
#' @return Named vector of display labels
get_display_labels <- function(meta, use_gene_name = FALSE) {
  if (use_gene_name && "GeneName" %in% colnames(meta)) {
    gene_name <- trimws(meta$GeneName)
    dupes <- duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE)
    gene_name[dupes] <- paste0(gene_name[dupes], "|", meta$GeneID[dupes])
    return(setNames(gene_name, meta$GeneID))
  } else {
    return(setNames(meta$GeneID, meta$GeneID))
  }
}

#' Normalize category labels
#' 
#' @param vec Input vector of categories
#' @return Normalized vector with NA/empty values replaced with "Other"
normalize_categories <- function(vec) {
  vec_names <- names(vec)
  vec <- as.character(vec)
  vec[vec == ""] <- NA
  vec[is.na(vec)] <- "Other"
  names(vec) <- vec_names
  return(vec)
}