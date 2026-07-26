###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 06_merge_seurat_objects.R
# Author  : Zeno Vimalan A.
# Date    : 26 July 2026
# Purpose : Merge all Seurat objects into one dataset
###########################################################

#----------------------------------------------------------
# 1. Clear workspace
#----------------------------------------------------------

rm(list = ls())

#----------------------------------------------------------
# 2. Load library
#----------------------------------------------------------

library(Seurat)

#----------------------------------------------------------
# 3. Set working directory
#----------------------------------------------------------

setwd("C:/Users/zvima/OneDrive/Documents/GitHub/scRNAseq_myasthenia_gravis")

#----------------------------------------------------------
# 4. Load Seurat object list
#----------------------------------------------------------

seurat_list <- readRDS(
  "RESULTS/seurat_object_list.rds"
)

#----------------------------------------------------------
# 5. Merge all Seurat objects
#----------------------------------------------------------

merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "Myasthenia_Gravis_scRNAseq"
)

#----------------------------------------------------------
# 6. Display summary
#----------------------------------------------------------

print(merged_seurat)

cat("\nNumber of cells:\n")
print(ncol(merged_seurat))

cat("\nNumber of genes:\n")
print(nrow(merged_seurat))

#----------------------------------------------------------
# 7. Save merged object
#----------------------------------------------------------

saveRDS(
  merged_seurat,
  file = "RESULTS/merged_seurat_object.rds"
)

cat("\nMerged Seurat object saved successfully!\n")
