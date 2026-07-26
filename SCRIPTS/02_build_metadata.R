###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 02_build_metadata.R
# Author  : Zeno Vimalan A.
# Purpose : Build sample metadata from filenames
###########################################################

rm(list = ls())

library(data.table)

setwd("C:/Bioinformatics/scRNAseq_myasthenia_gravis")

#----------------------------------------------------------
# List processed files
#----------------------------------------------------------

files <- list.files(
  "DATA/PROCESSED",
  pattern = "\\.txt\\.gz$"
)

#----------------------------------------------------------
# Create metadata table
#----------------------------------------------------------

metadata <- data.table(
  filename = files
)

# Sample IDs
metadata[, sample_id := sub("^GSM[0-9]+_", "", filename)]
metadata[, sample_id := sub("\\.txt\\.gz$", "", sample_id)]

#----------------------------------------------------------
# Assign biological groups
#----------------------------------------------------------

metadata[, group := fifelse(
  grepl("^A", sample_id),
  "AChR_Positive",
  fifelse(
    grepl("^H", sample_id),
    "Healthy_Control",
    fifelse(
      grepl("a$", sample_id),
      "Seronegative_Pre",
      "Seronegative_Post"
    )
  )
)]

#----------------------------------------------------------
# View metadata
#----------------------------------------------------------

print(metadata)

#----------------------------------------------------------
# Save metadata
#----------------------------------------------------------

fwrite(
  metadata,
  "RESULTS/metadata/sample_metadata.csv"
)

cat("\nMetadata saved successfully!\n")