###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 12_scaling.R
# Author  : Zeno Vimalan A.
# Date    : 27 July 2026
# Purpose : Scale highly variable genes
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
# 3. Load Seurat Object
#----------------------------------------------------------

variable_features_seurat <- readRDS(
  "RESULTS/objects/variable_features_seurat_object.rds"
)

#----------------------------------------------------------
# 4. Scale Data
#----------------------------------------------------------

cat("Scaling data...\n")

variable_features_seurat <- ScaleData(
  object = variable_features_seurat,
  features = VariableFeatures(variable_features_seurat),
  verbose = TRUE
)

cat("Scaling completed successfully.\n")

#----------------------------------------------------------
# 5. Save Scaled Object
#----------------------------------------------------------

saveRDS(
  variable_features_seurat,
  file = "RESULTS/objects/scaled_seurat_object.rds"
)

#----------------------------------------------------------
# 6. Summary
#----------------------------------------------------------

cat("\n---------------------------------------\n")
cat("Scaling Complete\n")
cat("---------------------------------------\n\n")

cat("Variable genes scaled :",
    length(VariableFeatures(variable_features_seurat)),
    "\n")

cat("Scaled object saved to:\n")
cat("RESULTS/objects/scaled_seurat_object.rds\n")