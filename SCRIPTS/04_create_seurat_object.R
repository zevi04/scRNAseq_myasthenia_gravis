###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 04_create_seurat_objects.R
# Author  : Zeno Vimalan A.
# Date    : 26 July 2026
# Purpose : Create one Seurat object for each sample
###########################################################

#----------------------------------------------------------
# 1. Clear workspace
#----------------------------------------------------------

rm(list = ls())

#----------------------------------------------------------
# 2. Load libraries
#----------------------------------------------------------

library(data.table)
library(Seurat)

#----------------------------------------------------------
# 3. Set working directory
#----------------------------------------------------------

setwd("C:/Users/zvima/OneDrive/Documents/GitHub/scRNAseq_myasthenia_gravis")

#----------------------------------------------------------
# 4. Read metadata
#----------------------------------------------------------

metadata <- fread("RESULTS/sample_metadata.csv")

#----------------------------------------------------------
# 5. Create an empty list
#----------------------------------------------------------

seurat_list <- list()

#----------------------------------------------------------
# 6. Read only the FIRST sample (for now)
#----------------------------------------------------------

sample_name <- metadata$sample_id[1]

sample_file <- metadata$filename[1]

cat("Processing:", sample_name, "\n")

counts <- fread(
  file.path("DATA/PROCESSED", sample_file)
)

#----------------------------------------------------------
# 7. Convert first column into row names
#----------------------------------------------------------

gene_names <- counts$ID

counts <- as.data.frame(counts)

rownames(counts) <- gene_names

counts$ID <- NULL

#----------------------------------------------------------
# 8. Convert to matrix
#----------------------------------------------------------

count_matrix <- as.matrix(counts)

#----------------------------------------------------------
# 9. Create Seurat object
#----------------------------------------------------------

seurat_object <- CreateSeuratObject(
  counts = count_matrix,
  project = sample_name
)

#----------------------------------------------------------
# 10. Store inside list
#----------------------------------------------------------

seurat_list[[sample_name]] <- seurat_object

#----------------------------------------------------------
# 11. Show summary
#----------------------------------------------------------

print(seurat_object)

