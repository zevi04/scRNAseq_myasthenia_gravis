###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 09_qc_filtering.R
# Author  : Zeno Vimalan A.
# Date    : 27 July 2026
# Purpose : Filter low-quality cells based on QC metrics
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
# 3. Load Merged Seurat Object
#----------------------------------------------------------

merged_seurat <- readRDS(
  "RESULTS/objects/merged_seurat_object.rds"
)

#----------------------------------------------------------
# 4. Record Number of Cells Before Filtering
#----------------------------------------------------------

cells_before <- ncol(merged_seurat)

cat("Cells before filtering :", cells_before, "\n")

#----------------------------------------------------------
# 5. Apply Quality Control Filters
#----------------------------------------------------------

filtered_seurat <- subset(
  merged_seurat,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 15
)

#----------------------------------------------------------
# 6. Record Number of Cells After Filtering
#----------------------------------------------------------

cells_after <- ncol(filtered_seurat)

cells_removed <- cells_before - cells_after

percent_removed <- round(
  (cells_removed / cells_before) * 100,
  2
)

cat("Cells after filtering  :", cells_after, "\n")
cat("Cells removed          :", cells_removed, "\n")
cat("Percent removed (%)    :", percent_removed, "\n")

#----------------------------------------------------------
# 7. Create QC Summary Table
#----------------------------------------------------------

qc_summary <- data.frame(
  
  Metric = c(
    "Cells before filtering",
    "Cells after filtering",
    "Cells removed",
    "Percent removed"
  ),
  
  Value = c(
    cells_before,
    cells_after,
    cells_removed,
    percent_removed
  )
  
)

#----------------------------------------------------------
# 8. Save QC Summary
#----------------------------------------------------------

write.csv(
  qc_summary,
  file = "RESULTS/tables/QC_Filtering_Summary.csv",
  row.names = FALSE
)

#----------------------------------------------------------
# 9. Save Filtered Seurat Object
#----------------------------------------------------------

saveRDS(
  filtered_seurat,
  file = "RESULTS/objects/filtered_seurat_object.rds"
)

#----------------------------------------------------------
# 10. Verify Files
#----------------------------------------------------------

cat("\nFiles saved successfully.\n")
cat("Filtered object : RESULTS/objects/filtered_seurat_object.rds\n")
cat("QC summary      : RESULTS/tables/QC_Filtering_Summary.csv\n")