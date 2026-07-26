###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 08_qc_visualization.R
# Author  : Zeno Vimalan A.
# Date    : 26 July 2026
# Purpose : Visualize quality control metrics
###########################################################

#----------------------------------------------------------
# 1. Clear workspace
#----------------------------------------------------------

rm(list = ls())

#----------------------------------------------------------
# 2. Load libraries
#----------------------------------------------------------

library(Seurat)
library(ggplot2)

#----------------------------------------------------------
# 3. Load merged Seurat object
#----------------------------------------------------------

merged_seurat <- readRDS(
  "RESULTS/objects/merged_seurat_object.rds"
)

#----------------------------------------------------------
# 4. Create output folder
#----------------------------------------------------------

dir.create("FIGURES/QC", recursive = TRUE, showWarnings = FALSE)

#----------------------------------------------------------
# 5. Violin plots
#----------------------------------------------------------

vln_plot <- VlnPlot(
  merged_seurat,
  features = c("nFeature_RNA",
               "nCount_RNA",
               "percent.mt"),
  pt.size = 0,
  ncol = 3
)

vln_plot <- vln_plot &
  theme_classic(base_size = 14) &
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8
    ),
    axis.title.x = element_blank(),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    legend.position = "none"
  )

ggsave(
  filename = "FIGURES/QC/QC_Violin_Plot.png",
  plot = vln_plot,
  width = 16,
  height = 6,
  dpi = 300
)

#----------------------------------------------------------
# 6. Scatter plot: Counts vs Features
#----------------------------------------------------------

scatter1 <- FeatureScatter(
  merged_seurat,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

ggsave(
  filename = "FIGURES/QC/Counts_vs_Features.png",
  plot = scatter1,
  width = 6,
  height = 5,
  dpi = 300
)

#----------------------------------------------------------
# 7. Scatter plot: Counts vs Mitochondrial %
#----------------------------------------------------------

scatter2 <- FeatureScatter(
  merged_seurat,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

ggsave(
  filename = "FIGURES/QC/Counts_vs_Mito.png",
  plot = scatter2,
  width = 6,
  height = 5,
  dpi = 300
)

#----------------------------------------------------------
# 8. Scatter plot: Features vs Mitochondrial %
#----------------------------------------------------------

scatter3 <- FeatureScatter(
  merged_seurat,
  feature1 = "nFeature_RNA",
  feature2 = "percent.mt"
)

ggsave(
  filename = "FIGURES/QC/Features_vs_Mito.png",
  plot = scatter3,
  width = 6,
  height = 5,
  dpi = 300
)

#----------------------------------------------------------
# 9. Display plots in RStudio
#----------------------------------------------------------

print(vln_plot)
print(scatter1)
print(scatter2)
print(scatter3)

cat("\nQC plots generated successfully!\n")
cat("Saved in: FIGURES/QC/\n")