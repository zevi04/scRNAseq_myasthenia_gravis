###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 03_check_dataset_consistency.R
# Author  : Zeno Vimalan A.
# Date    : 26 July 2026
# Purpose : Check whether all samples have the same genes
###########################################################

#----------------------------------------------------------
# 1. Clear workspace
#----------------------------------------------------------

rm(list = ls())

#----------------------------------------------------------
# 2. Load package
#----------------------------------------------------------

library(data.table)

#----------------------------------------------------------
# 3. Set working directory
#----------------------------------------------------------

setwd("C:/Users/zvima/OneDrive/Documents/GitHub/scRNAseq_myasthenia_gravis")

#----------------------------------------------------------
# 4. Read sample metadata
#----------------------------------------------------------

metadata <- fread("RESULTS/sample_metadata.csv")

#----------------------------------------------------------
# 5. Read first column (gene IDs) from every sample
#----------------------------------------------------------

gene_lists <- lapply(metadata$filename, function(file){
  
  fread(
    file.path("DATA/PROCESSED", file),
    select = 1
  )$ID
  
})

#----------------------------------------------------------
# 6. Compare gene lists
#----------------------------------------------------------

reference <- gene_lists[[1]]

same_genes <- sapply(gene_lists, function(x){
  
  identical(reference, x)
  
})

#----------------------------------------------------------
# 7. Display results
#----------------------------------------------------------

print(same_genes)

cat("\n")

cat(sum(same_genes), "out of", length(same_genes),
    "samples have identical gene lists.\n")