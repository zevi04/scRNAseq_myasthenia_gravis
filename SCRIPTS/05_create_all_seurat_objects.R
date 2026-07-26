###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 05_create_all_seurat_objects.R
# Author  : Zeno Vimalan A.
# Date    : 26 July 2026
# Purpose : Create one Seurat object for every sample
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
# 6. Loop through every sample
#----------------------------------------------------------

for(i in 1:nrow(metadata))
{
  
  sample_name <- metadata$sample_id[i]
  sample_file <- metadata$filename[i]
  
  cat("Processing:", sample_name, "\n")
  
  # Read expression matrix
  counts <- fread(
    file.path("DATA/PROCESSED", sample_file)
  )
  
  # Extract gene names
  gene_names <- counts$ID
  
  # Convert to data frame
  counts <- as.data.frame(counts)
  
  # Set row names
  rownames(counts) <- gene_names
  
  # Remove gene column
  counts$ID <- NULL
  
  # Convert to matrix
  count_matrix <- as.matrix(counts)
  
  # Create Seurat object
  seurat_object <- CreateSeuratObject(
    counts = count_matrix,
    project = sample_name
  )
  
  # Store object
  seurat_list[[sample_name]] <- seurat_object
  
}

#----------------------------------------------------------
# 7. Summary
#----------------------------------------------------------

cat("\n---------------------------------\n")
cat("Total Seurat Objects Created :", length(seurat_list), "\n")
cat("---------------------------------\n")

print(names(seurat_list))

#----------------------------------------------------------
# 8. Save objects
#----------------------------------------------------------

saveRDS(
  seurat_list,
  file = "RESULTS/seurat_object_list.rds"
)

cat("\nAll Seurat objects saved successfully!\n")