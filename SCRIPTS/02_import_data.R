###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 01_import_data.R
# Author  : Zeno Vimalan A.
# Date    : 26 July 2026
# Purpose : Import and inspect one processed scRNA-seq sample
###########################################################

#----------------------------------------------------------
# 1. Clear workspace
#----------------------------------------------------------

rm(list = ls())

#----------------------------------------------------------
# 2. Load required package
#----------------------------------------------------------

library(data.table)

#----------------------------------------------------------
# 3. Set working directory
#----------------------------------------------------------

setwd("C:/Users/zvima/OneDrive/Documents/GitHub/scRNAseq_myasthenia_gravis")

#----------------------------------------------------------
# 4. Define the file to import
#----------------------------------------------------------

sample_file <- "DATA/PROCESSED/GSM7266236_A1.txt.gz"

#----------------------------------------------------------
# 5. Check whether the file exists
#----------------------------------------------------------

if (!file.exists(sample_file)) {
  stop("ERROR: Sample file not found!")
}

cat("Sample file found.\n\n")

#----------------------------------------------------------
# 6. Import a small portion of the dataset
#    (First 10 genes only)
#----------------------------------------------------------

sample_data <- fread(
  sample_file,
  nrows = 10
)

cat("Data imported successfully.\n\n")

#----------------------------------------------------------
# 7. Explore the imported data
#----------------------------------------------------------

cat("Class:\n")
print(class(sample_data))

cat("\nDimensions:\n")
print(dim(sample_data))

cat("\nFirst 10 column names:\n")
print(colnames(sample_data)[1:10])

cat("\nFirst 6 rows and first 5 columns:\n")
print(sample_data[, 1:5])

#----------------------------------------------------------
# End of Script
#----------------------------------------------------------