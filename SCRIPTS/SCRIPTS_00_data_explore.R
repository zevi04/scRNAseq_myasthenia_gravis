###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 00_dataset_exploration.R
# Author  : Zeno Vimalan A.
# Purpose : Explore the downloaded GEO dataset
###########################################################

# Clear workspace
rm(list = ls())

# Load package
library(data.table)

# Set working directory
setwd("C:/Users/zvima/OneDrive/Documents/GitHub/scRNAseq_myasthenia_gravis")

# Check current directory
getwd()

# Show project folders
list.files()

# List all processed files
files <- list.files(
  path = "DATA/PROCESSED",
  pattern = "\\.txt\\.gz$"
)

# Number of files
length(files)

# First 10 files
head(files, 10)

# Check that the first file exists
file.exists(file.path("DATA/PROCESSED", files[1]))


