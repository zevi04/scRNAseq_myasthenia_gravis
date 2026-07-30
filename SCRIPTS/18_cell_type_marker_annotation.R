###############################################################
# SCRIPT 18: CELL TYPE MARKER IDENTIFICATION
# (Memory-Optimized Version)
###############################################################

rm(list = ls())
gc()

###############################################################
# Load Libraries
###############################################################

library(Seurat)
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

###############################################################
# Load Annotated Seurat Object
###############################################################

scaled_seurat <- readRDS(
  "D:/Research/scRNAseq/RESULTS/objects/GSE227835_annotated.rds"
)

###############################################################
# Set Cell Type Identities
###############################################################

Idents(scaled_seurat) <- "SingleR_Label"

cat("\nCell Type Distribution:\n")
print(table(Idents(scaled_seurat)))

###############################################################
# Downsample Cells
###############################################################

set.seed(123)

marker_object <- subset(
  scaled_seurat,
  downsample = 3000
)

cat("\nCells after downsampling:\n")
print(table(Idents(marker_object)))

###############################################################
# Find Marker Genes
###############################################################

markers <- FindAllMarkers(
  marker_object,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = TRUE
)

###############################################################
# Save Complete Marker Table
###############################################################

write.csv(
  markers,
  "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS/tables/CellType_Markers_All.csv",
  row.names = FALSE
)

###############################################################
# Top 5 Marker Genes
###############################################################

top5 <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(
    avg_log2FC,
    n = 5,
    with_ties = FALSE
  )

write.csv(
  top5,
  "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS/tables/Top10_Markers_Per_CellType.csv",
  row.names = FALSE
)

###############################################################
# Dot Plot
###############################################################

top.genes <- unique(top5$gene)

p1 <- DotPlot(
  marker_object,
  features = top.genes
) +
  RotatedAxis() +
  ggtitle("Top Marker Genes")

print(p1)

ggsave(
  "C:/Bioinformatics/scRNAseq_myasthenia_gravis/FIGURES/DotPlot_TopMarkers.png",
  plot = p1,
  width = 16,
  height = 8,
  dpi = 600
)

###############################################################
# Scale Only Marker Genes
###############################################################

marker_object <- ScaleData(
  marker_object,
  features = top.genes,
  verbose = FALSE
)

###############################################################
# Heatmap
###############################################################

p2 <- DoHeatmap(
  marker_object,
  features = top.genes,
  size = 4,
  group.bar = TRUE,
  draw.lines = TRUE
) +
  ggtitle("Top Marker Genes")

print(p2)

ggsave(
  "C:/Bioinformatics/scRNAseq_myasthenia_gravis/FIGURES/Heatmap_TopMarkers.png",
  plot = p2,
  width = 14,
  height = 10,
  dpi = 600
)

###############################################################
# Summary
###############################################################

cat("\nMarker identification completed successfully.\n")

cat("\nNumber of marker genes:\n")
print(nrow(markers))

cat("\nTop markers per cell type:\n")
print(top5)