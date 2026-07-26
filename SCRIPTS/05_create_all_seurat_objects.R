###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 05_create_all_seurat_objects.R
# Author  : Zeno Vimalan A.
# Date    : 26 July 2026
# Purpose : Create Seurat objects for all samples
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

setwd("C:/Bioinformatics/scRNAseq_myasthenia_gravis")

#----------------------------------------------------------
# 4. Read metadata
#----------------------------------------------------------

metadata <- fread("RESULTS/metadata/sample_metadata.csv")

#----------------------------------------------------------
# 5. Create empty list
#----------------------------------------------------------

seurat_list <- list()

#----------------------------------------------------------
# 6. Loop through every sample
#----------------------------------------------------------

for(i in seq_len(nrow(metadata))) {
  
  sample_name <- metadata$sample_id[i]
  sample_file <- metadata$filename[i]
  
  cat("Processing:", sample_name, "\n")
  
  counts <- fread(file.path("DATA/PROCESSED", sample_file))
  
  gene_names <- counts$ID
  
  counts <- as.data.frame(counts)
  
  rownames(counts) <- gene_names
  
  counts$ID <- NULL
  
  count_matrix <- as.matrix(counts)
  
  seurat_list[[sample_name]] <- CreateSeuratObject(
    counts = count_matrix,
    project = sample_name
  )
  
}

cat("\n")
cat(length(seurat_list), "Seurat objects successfully created.\n")
