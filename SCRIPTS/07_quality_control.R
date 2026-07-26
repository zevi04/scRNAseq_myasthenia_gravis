###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 07_quality_control.R
# Author  : Zeno Vimalan A.
# Purpose : Calculate quality control metrics
###########################################################

#----------------------------------------------------------
# 1. Clear workspace
#----------------------------------------------------------

rm(list = ls())

#----------------------------------------------------------
# 2. Load libraries
#----------------------------------------------------------

library(Seurat)
library(ggplot2)

#----------------------------------------------------------
# 3. Load merged Seurat object
#----------------------------------------------------------

merged_seurat <- readRDS(
  "RESULTS/objects/merged_seurat_object.rds"
)

#----------------------------------------------------------
# 4. Calculate mitochondrial percentage
#----------------------------------------------------------

merged_seurat[["percent.mt"]] <- PercentageFeatureSet(
  merged_seurat,
  pattern = "^MT-"
)

#----------------------------------------------------------
# 5. Display summaries
#----------------------------------------------------------

cat("\n===== Number of Genes per Cell =====\n")
print(summary(merged_seurat$nFeature_RNA))

cat("\n===== Number of UMIs per Cell =====\n")
print(summary(merged_seurat$nCount_RNA))

cat("\n===== Percent Mitochondrial Reads =====\n")
print(summary(merged_seurat$percent.mt))

#----------------------------------------------------------
# 6. Save updated Seurat object
#----------------------------------------------------------

saveRDS(
  merged_seurat,
  "RESULTS/objects/merged_seurat_object.rds"
)

cat("\nQC metrics calculated and saved successfully.\n")