###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 14_UMAP_analysis.R
# Author  : Zeno Vimalan A.
# Date    : 28 July 2026
# Purpose : Perform UMAP dimensionality reduction
###########################################################

#----------------------------------------------------------
# Check that PCA has been completed
#----------------------------------------------------------

if (!exists("scaled_seurat")) {
  stop("Run Script 13 (PCA) first.")
}

if (!"pca" %in% Reductions(scaled_seurat)) {
  stop("PCA reduction not found.")
}

#----------------------------------------------------------
# Create Output Directory
#----------------------------------------------------------

dir.create(
  "FIGURES/UMAP",
  recursive = TRUE,
  showWarnings = FALSE
)

#----------------------------------------------------------
# Run UMAP
#----------------------------------------------------------

cat("Running UMAP...\n")

scaled_seurat <- RunUMAP(
  object = scaled_seurat,
  dims = 1:25,
  reduction = "pca",
  metric = "cosine",
  verbose = TRUE
)

cat("UMAP completed successfully.\n")

gc()

#----------------------------------------------------------
# Generate UMAP Plot
#----------------------------------------------------------

cat("Generating UMAP plot...\n")

umap_plot <- DimPlot(
  object = scaled_seurat,
  reduction = "umap",
  pt.size = 0.2
) +
  ggtitle("UMAP Projection of Single Cells") +
  theme_classic(base_size = 14)

ggsave(
  filename = "FIGURES/UMAP/UMAP.tiff",
  plot = umap_plot,
  width = 8,
  height = 6,
  dpi = 600,
  bg = "white"
)

rm(umap_plot)
gc()

#----------------------------------------------------------
# Summary
#----------------------------------------------------------

cat("\n---------------------------------------\n")
cat("UMAP Analysis Complete\n")
cat("---------------------------------------\n\n")

cat("UMAP successfully generated.\n")
cat("Figure saved to:\n")
cat("FIGURES/UMAP/UMAP.tiff\n")