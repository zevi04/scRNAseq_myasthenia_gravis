###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 15_Find_Neighbors.R
# Author  : Zeno Vimalan A.
# Date    : 28 July 2026
# Purpose : Construct Shared Nearest Neighbor (SNN) Graph
###########################################################

#----------------------------------------------------------
# Check that UMAP has been completed
#----------------------------------------------------------

if (!"umap" %in% Reductions(scaled_seurat)) {
  stop("Run Script 14 (UMAP) first.")
}

#----------------------------------------------------------
# Find Nearest Neighbors
#----------------------------------------------------------

cat("Constructing Shared Nearest Neighbor (SNN) graph...\n")

scaled_seurat <- FindNeighbors(
  object = scaled_seurat,
  reduction = "pca",
  dims = 1:25,
  verbose = TRUE
)

gc()

#----------------------------------------------------------
# Summary
#----------------------------------------------------------

cat("\n---------------------------------------\n")
cat("Neighbor Graph Construction Complete\n")
cat("---------------------------------------\n\n")

cat("Using Principal Components : 1-25\n")
cat("Neighbor graph successfully constructed.\n")
cat("Proceed directly to clustering.\n")