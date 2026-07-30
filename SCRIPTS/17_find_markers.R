###############################################################
# SCRIPT 17: CELL TYPE ANNOTATION USING SingleR
###############################################################

rm(list = ls())
gc()

###############################################################
# Load Libraries
###############################################################

library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)

###############################################################
# Create Output Folders
###############################################################

dir.create(
  "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "D:/Research/scRNAseq/RESULTS/objects",
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################
# Load Seurat Object
###############################################################

scaled_seurat <- readRDS(
  "D:/Research/scRNAseq/RESULTS/objects/GSE227835_joined_layers.rds"
)

###############################################################
# Aggregate Expression by Cluster
###############################################################

avg.exp <- AggregateExpression(
  object = scaled_seurat,
  assays = "RNA",
  group.by = "seurat_clusters",
  return.seurat = FALSE
)

avg.mat <- avg.exp$RNA

cat("Aggregated matrix dimensions:\n")
print(dim(avg.mat))
print(object.size(avg.mat))

###############################################################
# Create SingleCellExperiment
###############################################################

sce.cluster <- SingleCellExperiment(
  assays = list(logcounts = avg.mat)
)

###############################################################
# Load Monaco Immune Reference
###############################################################

monaco.ref <- MonacoImmuneData()

###############################################################
# Run SingleR
###############################################################

pred <- SingleR(
  test = sce.cluster,
  ref = monaco.ref,
  labels = monaco.ref$label.main
)

###############################################################
# Create Annotation Table
###############################################################

annotation <- data.frame(
  Cluster = colnames(avg.mat),
  Predicted_Cell_Type = pred$labels,
  Delta_Score = pred$delta.next
)

###############################################################
# Remove 'g' Prefix from Cluster Names
###############################################################

annotation$Cluster <- sub("^g", "", annotation$Cluster)

print(annotation)

###############################################################
# Save Annotation Table
###############################################################

write.csv(
  annotation,
  "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS/tables/SingleR_cluster_annotation.csv",
  row.names = FALSE
)

###############################################################
# Map Cell Type Labels to Seurat Object
###############################################################

cluster.map <- setNames(
  annotation$Predicted_Cell_Type,
  annotation$Cluster
)

labels <- unname(
  cluster.map[as.character(Idents(scaled_seurat))]
)

scaled_seurat <- AddMetaData(
  object = scaled_seurat,
  metadata = labels,
  col.name = "SingleR_Label"
)

###############################################################
# Quality Control
###############################################################

cat("\nCell type distribution:\n")
print(table(scaled_seurat$SingleR_Label))

###############################################################
# UMAP Visualization
###############################################################

p1 <- DimPlot(
  scaled_seurat,
  group.by = "SingleR_Label",
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  raster = TRUE
) +
  ggtitle("SingleR Cell Type Annotation")

print(p1)

ggsave(
  filename = "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS/figures/UMAP_SingleR_Annotation.png",
  plot = p1,
  width = 12,
  height = 10,
  dpi = 600
)

###############################################################
# Save Annotated Seurat Object
###############################################################

saveRDS(
  scaled_seurat,
  "D:/Research/scRNAseq/RESULTS/objects/GSE227835_annotated.rds",
  compress = FALSE
)

cat("\nSingleR annotation completed successfully.\n")