###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 06_merge_seurat_objects.R
# Author  : Zeno Vimalan A.
# Date    : 26 July 2026
# Purpose : Merge all Seurat objects into a single dataset
###########################################################

#----------------------------------------------------------
# 1. Check that Script 05 has been run
#----------------------------------------------------------

if (!exists("seurat_list")) {
  stop("Run 05_create_all_seurat_objects.R first.")
}

#----------------------------------------------------------
# 2. Merge all Seurat objects
#----------------------------------------------------------

cat("Merging Seurat objects...\n")

merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.ids = names(seurat_list),
  project = "Myasthenia_Gravis_scRNAseq"
)

cat("Merge completed successfully!\n")

#----------------------------------------------------------
# 3. Display summary
#----------------------------------------------------------

cat("\nMerged Seurat object created successfully!\n\n")

print(merged_seurat)

cat("\nGenes :", nrow(merged_seurat), "\n")
cat("Cells :", ncol(merged_seurat), "\n")

#----------------------------------------------------------
# 4. Save merged Seurat object
#----------------------------------------------------------

saveRDS(
  object = merged_seurat,
  file = "RESULTS/merged_seurat_object.rds"
)

cat("\nMerged Seurat object saved successfully!\n")

#----------------------------------------------------------
# 5. Verify that the file was saved
#----------------------------------------------------------

if (file.exists("RESULTS/merged_seurat_object.rds")) {
  cat("Verified: merged_seurat_object.rds exists.\n")
} else {
  cat("ERROR: merged_seurat_object.rds was not found.\n")
} 