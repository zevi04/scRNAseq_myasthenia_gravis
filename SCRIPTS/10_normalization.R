###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 10_normalization.R
# Author  : Zeno Vimalan A.
# Date    : 27 July 2026
# Purpose : Normalize filtered scRNA-seq data
###########################################################

#----------------------------------------------------------
# 1. Clear Workspace
#----------------------------------------------------------

rm(list = ls())

#----------------------------------------------------------
# 2. Load Libraries
#----------------------------------------------------------

library(Seurat)

#----------------------------------------------------------
# 3. Load Filtered Seurat Object
#----------------------------------------------------------

filtered_seurat <- readRDS(
  "RESULTS/objects/filtered_seurat_object.rds"
)

#----------------------------------------------------------
# 4. Normalize Data
#----------------------------------------------------------

cat("Starting normalization...\n")

filtered_seurat <- NormalizeData(
  object = filtered_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

cat("Normalization completed successfully.\n")

#----------------------------------------------------------
# 5. Save Normalized Seurat Object
#----------------------------------------------------------

saveRDS(
  filtered_seurat,
  file = "RESULTS/objects/normalized_seurat_object.rds"
)

#----------------------------------------------------------
# 6. Verify Output
#----------------------------------------------------------

cat("\nNormalized Seurat object saved successfully.\n")
cat("Location : RESULTS/objects/normalized_seurat_object.rds\n")