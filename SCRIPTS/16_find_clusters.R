###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 16_Find_Clusters.R
# Author  : Zeno Vimalan A.
# Date    : 28 July 2026
# Purpose : Identify cell clusters
###########################################################

#----------------------------------------------------------
# Check Neighbor Graph
#----------------------------------------------------------

if (!"RNA_snn" %in% Graphs(scaled_seurat)) {
  stop("Run Script 15 (FindNeighbors) first.")
}

#----------------------------------------------------------
# Create Output Directory
#----------------------------------------------------------

dir.create(
  "FIGURES/CLUSTERING",
  recursive = TRUE,
  showWarnings = FALSE
)

#----------------------------------------------------------
# Identify Cell Clusters
#----------------------------------------------------------

cat("Running graph-based clustering...\n")

scaled_seurat <- FindClusters(
  object = scaled_seurat,
  resolution = 0.8,
  algorithm = 4,
  verbose = TRUE
)

cat("Clustering completed successfully.\n")

gc()

#----------------------------------------------------------
# Generate Cluster UMAP
#----------------------------------------------------------

cluster_plot <- DimPlot(
  object = scaled_seurat,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.1,
  raster = TRUE
) +
  ggtitle("UMAP of Cell Clusters") +
  theme_classic(base_size = 14)

ggsave(
  filename = "FIGURES/CLUSTERING/UMAP_Clusters.tiff",
  plot = cluster_plot,
  width = 8,
  height = 6,
  dpi = 600,
  bg = "white"
)

rm(cluster_plot)
gc()

#----------------------------------------------------------
# Display Cluster Information
#----------------------------------------------------------

cat("\n---------------------------------------\n")
cat("Cluster Summary\n")
cat("---------------------------------------\n\n")

print(table(Idents(scaled_seurat)))

cat("\nNumber of Clusters :",
    length(levels(Idents(scaled_seurat))),
    "\n")

cat("\nCluster UMAP saved to:\n")
cat("FIGURES/CLUSTERING/UMAP_Clusters.tiff\n")

cat("\nProceed to marker gene identification.\n")