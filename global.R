library(shiny)
library(shinydashboard)
library(data.table)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(circlize)
library(shinycssloaders)
library(Cairo)
library(DT)
library(ggplot2)
library(shinyjs)
library(heatmaply)
library(plotly)
library(ggnewscale)
library(Seurat)

file_paths <- list(
  "Kelsey coexpression_bulk" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/cor_ordered_bulk.rds",
  "Kelsey coexpression_all" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/cor_ordered_all.rds",
  "Tylor coexpression_bulk" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/tylor_cor_ordered_bulk.rds",
  "Tylor coexpression_all" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/tylor_cor_ordered_all.rds",
  "Coexpression_all" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/dataset1.rds",
  "Coexpression_bulk" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/dataset2.rds",
  "ATH Leaf 2k Pearson" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/2k_pearson.rds",
  "ATH Leaf 2k Spearman" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/2k_spearman.rds",
  "Cauline" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/cauline.rds",
  "flower" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/flower.rds",
  "leaf" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/leaf.rds",
  "root" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/root.rds",
  "shoot" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/shoot.rds",
  "silique" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/silique.rds",
  "stem" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/stem.rds"
)


seurat_paths <- list(
  "Cauline" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/seurat/cauline_ne_seurat.rds",
  "flower" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/seurat/flower_ne_seurat.rds",
  "leaf" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/seurat/leaf_ne_seurat.rds",
  "root" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/seurat/root_ne_seurat.rds",
  "shoot" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/seurat/shoot_ne_seurat.rds",
  "silique" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/seurat/silique_ne_seurat.rds",
  "stem" = "https://raw.githubusercontent.com/krishna-dhulipalla/gene_co-expression_tool/master/data/seurat/stem_ne_seurat.rds"
)
