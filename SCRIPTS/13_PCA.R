###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 13_PCA_analysis.R
# Author  : Zeno Vimalan A.
# Date    : 28 July 2026
# Purpose : Perform Principal Component Analysis (PCA)
###########################################################

#----------------------------------------------------------
# 1. Clear Workspace
#----------------------------------------------------------

rm(list = ls())
gc()

#----------------------------------------------------------
# 2. Load Libraries
#----------------------------------------------------------

library(Seurat)
library(ggplot2)

#----------------------------------------------------------
# 3. Load Scaled Seurat Object
#----------------------------------------------------------

cat("Loading scaled Seurat object...\n")

scaled_seurat <- readRDS(
  "RESULTS/objects/scaled_seurat_object.rds"
)

gc()

#----------------------------------------------------------
# 4. Create Output Directory
#----------------------------------------------------------

dir.create(
  "FIGURES/PCA",
  recursive = TRUE,
  showWarnings = FALSE
)

#----------------------------------------------------------
# 5. Run PCA
#----------------------------------------------------------

cat("Running PCA...\n")

scaled_seurat <- RunPCA(
  object = scaled_seurat,
  features = VariableFeatures(scaled_seurat),
  npcs = 50,
  approx = TRUE,
  verbose = TRUE
)

cat("PCA completed successfully.\n")

gc()

#----------------------------------------------------------
# 6. PCA Scatter Plot
#----------------------------------------------------------

cat("Generating PCA scatter plot...\n")

pca_plot <- DimPlot(
  scaled_seurat,
  reduction = "pca"
) +
  ggtitle("Principal Component Analysis") +
  theme_classic(base_size = 14)

ggsave(
  filename = "FIGURES/PCA/PCA_Scatter_Plot.tiff",
  plot = pca_plot,
  width = 7,
  height = 6,
  dpi = 600,
  bg = "white"
)

rm(pca_plot)
gc()

#----------------------------------------------------------
# 7. Elbow Plot
#----------------------------------------------------------

cat("Generating Elbow Plot...\n")

elbow_plot <- ElbowPlot(
  scaled_seurat,
  ndims = 50
) +
  ggtitle("Elbow Plot") +
  theme_classic(base_size = 14)

ggsave(
  filename = "FIGURES/PCA/Elbow_Plot.tiff",
  plot = elbow_plot,
  width = 7,
  height = 6,
  dpi = 600,
  bg = "white"
)

rm(elbow_plot)
gc()

#----------------------------------------------------------
# 8. Summary
#----------------------------------------------------------

cat("\n---------------------------------------\n")
cat("PCA Analysis Complete\n")
cat("---------------------------------------\n\n")

cat("Cells :", ncol(scaled_seurat), "\n")
cat("Genes :", nrow(scaled_seurat), "\n")
cat("Principal Components : 50\n\n")

cat("Generated Files:\n")
cat(" - FIGURES/PCA/PCA_Scatter_Plot.tiff\n")
cat(" - FIGURES/PCA/Elbow_Plot.tiff\n")

cat("\nPCA completed successfully.\n")
cat("Proceed directly to UMAP.\n")