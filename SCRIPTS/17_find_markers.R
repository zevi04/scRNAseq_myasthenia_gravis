###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 17_Find_Markers.R
# Author  : Zeno Vimalan A.
# Date    : 28 July 2026
# Purpose : Identify marker genes for each cell cluster
###########################################################

#----------------------------------------------------------
# 1. Load Libraries
#----------------------------------------------------------

library(Seurat)
library(dplyr)
library(future)

#----------------------------------------------------------
# 2. Enable Parallel Processing
#----------------------------------------------------------

plan(multisession, workers = 12)

#----------------------------------------------------------
# 3. Join Normalized Data Layers (IMPORTANT)
#----------------------------------------------------------
#
#cat("Joining normalized data layers...\n")
#
#scaled_seurat[["RNA"]] <- JoinLayers(
#  object = scaled_seurat[["RNA"]],
#  layers = "data",
 # new = "data"
#)
#
#cat("Normalized data layers joined successfully.\n\n")
#
#----------------------------------------------------------
# 4. Create Output Directory
#----------------------------------------------------------

dir.create(
  "RESULTS/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

#----------------------------------------------------------
# 5. Identify Marker Genes
#----------------------------------------------------------

cat("---------------------------------------\n")
cat("Finding Marker Genes\n")
cat("---------------------------------------\n\n")

markers <- FindAllMarkers(
  object = scaled_seurat,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  return.thresh = 0.05,
  verbose = TRUE
)

cat("\nMarker gene identification completed.\n")

#----------------------------------------------------------
# 6. Save Complete Marker Table
#----------------------------------------------------------

write.csv(
  markers,
  file = "RESULTS/tables/All_Cluster_Markers.csv",
  row.names = FALSE
)

#----------------------------------------------------------
# 7. Save Top 20 Marker Genes Per Cluster
#----------------------------------------------------------

top20_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 20,
    with_ties = FALSE
  )

write.csv(
  top20_markers,
  file = "RESULTS/tables/Top20_Markers_Per_Cluster.csv",
  row.names = FALSE
)

#----------------------------------------------------------
# 8. Save Highly Significant Markers
#----------------------------------------------------------

significant_markers <- markers %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC > 1
  )

write.csv(
  significant_markers,
  file = "RESULTS/tables/Significant_Markers.csv",
  row.names = FALSE
)

#----------------------------------------------------------
# 9. Restore Sequential Processing
#----------------------------------------------------------

plan(sequential)

gc()

#----------------------------------------------------------
# 10. Summary
#----------------------------------------------------------

cat("\n---------------------------------------\n")
cat("Marker Gene Identification Complete\n")
cat("---------------------------------------\n\n")

cat("Total marker genes identified :", nrow(markers), "\n")
cat("Number of clusters           :", length(unique(markers$cluster)), "\n\n")

cat("Files generated:\n")
cat("1. RESULTS/tables/All_Cluster_Markers.csv\n")
cat("2. RESULTS/tables/Top20_Markers_Per_Cluster.csv\n")
cat("3. RESULTS/tables/Significant_Markers.csv\n")

cat("\nProceed to cell type annotation.\n")