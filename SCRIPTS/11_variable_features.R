###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 11_variable_features.R
# Author  : Zeno Vimalan A.
# Date    : 27 July 2026
# Purpose : Identify highly variable genes
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
# 3. Load Normalized Seurat Object
#----------------------------------------------------------

normalized_seurat <- readRDS(
  "RESULTS/objects/normalized_seurat_object.rds"
)

#----------------------------------------------------------
# 4. Identify Highly Variable Features
#----------------------------------------------------------

cat("Identifying highly variable genes...\n")

normalized_seurat <- FindVariableFeatures(
  object = normalized_seurat,
  selection.method = "vst",
  nfeatures = 2000
)

cat("Variable feature selection completed.\n")

#----------------------------------------------------------
# 5. Extract Variable Genes
#----------------------------------------------------------

variable_genes <- VariableFeatures(normalized_seurat)

#----------------------------------------------------------
# 6. Summary
#----------------------------------------------------------

cat("\n---------------------------------------\n")
cat("Variable Feature Selection Complete\n")
cat("---------------------------------------\n\n")

cat("Total genes             :", nrow(normalized_seurat), "\n")
cat("Variable genes selected :", length(variable_genes), "\n")

#----------------------------------------------------------
# 7. Save Variable Gene List
#----------------------------------------------------------

write.csv(
  data.frame(Gene = variable_genes),
  file = "RESULTS/tables/highly_variable_genes.csv",
  row.names = FALSE
)

#----------------------------------------------------------
# 8. Save Updated Seurat Object
#----------------------------------------------------------

saveRDS(
  normalized_seurat,
  file = "RESULTS/objects/variable_features_seurat_object.rds"
)

#----------------------------------------------------------
# 9. Completion Message
#----------------------------------------------------------

cat("\nFiles saved successfully.\n")

cat("Variable genes : RESULTS/tables/highly_variable_genes.csv\n")

cat("Seurat object  : RESULTS/objects/variable_features_seurat_object.rds\n")