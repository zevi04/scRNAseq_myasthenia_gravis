###########################################################
# Project : scRNA-seq Analysis of Myasthenia Gravis
# Script  : 08_qc_visualization.R
# Author  : Zeno Vimalan A.
# Purpose : Generate Quality Control Figures
###########################################################

#----------------------------------------------------------
# 1. Clear Workspace
#----------------------------------------------------------

rm(list = ls())

#----------------------------------------------------------
# 2. Load Libraries
#----------------------------------------------------------

library(Seurat)
library(ggplot2)

#----------------------------------------------------------
# 3. Load Seurat Object
#----------------------------------------------------------

merged_seurat <- readRDS(
  "RESULTS/objects/merged_seurat_object.rds"
)

#----------------------------------------------------------
# 4. Create Output Folder
#----------------------------------------------------------

dir.create(
  "FIGURES/QC",
  recursive = TRUE,
  showWarnings = FALSE
)

#----------------------------------------------------------
# 5. Common Theme
#----------------------------------------------------------

qc_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    axis.text = element_text(
      colour = "black",
      size = 10
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(12,12,12,12)
  )

#----------------------------------------------------------
# 6. Save Function
#----------------------------------------------------------

save_plot <- function(plot, filename){
  
  ggsave(
    filename = filename,
    plot = plot,
    width = 8,
    height = 6,
    dpi = 600,
    bg = "white"
  )
  
}

#----------------------------------------------------------
# 7. Violin Plot Function
#----------------------------------------------------------

create_violin_plot <- function(feature,
                               title,
                               ylab,
                               filename){
  
  p <- VlnPlot(
    merged_seurat,
    features = feature,
    pt.size = 0
  ) +
    labs(
      title = title,
      x = "",
      y = ylab
    ) +
    qc_theme
  
  save_plot(p, filename)
  
  return(p)
  
}

#----------------------------------------------------------
# 8. Scatter Plot Function
#----------------------------------------------------------

create_scatter_plot <- function(feature1,
                                feature2,
                                title,
                                xlabel,
                                ylabel,
                                filename){
  
  p <- FeatureScatter(
    merged_seurat,
    feature1 = feature1,
    feature2 = feature2
  ) +
    labs(
      title = title,
      x = xlabel,
      y = ylabel
    ) +
    qc_theme
  
  save_plot(p, filename)
  
  return(p)
  
}

#----------------------------------------------------------
# 9. Figure 1
#----------------------------------------------------------

genes_plot <- create_violin_plot(
  
  feature = "nFeature_RNA",
  
  title = "Distribution of Detected Genes per Cell",
  
  ylab = "Detected Genes",
  
  filename = "FIGURES/QC/Figure_01_Distribution_of_Detected_Genes_per_Cell.tif"
  
)

#----------------------------------------------------------
# 10. Figure 2
#----------------------------------------------------------

umi_plot <- create_violin_plot(
  
  feature = "nCount_RNA",
  
  title = "Distribution of UMI Counts per Cell",
  
  ylab = "UMI Counts",
  
  filename = "FIGURES/QC/Figure_02_Distribution_of_UMI_Counts_per_Cell.tif"
  
)

#----------------------------------------------------------
# 11. Figure 3
#----------------------------------------------------------

mito_plot <- create_violin_plot(
  
  feature = "percent.mt",
  
  title = "Distribution of Mitochondrial Gene Percentage",
  
  ylab = "Mitochondrial Gene Percentage (%)",
  
  filename = "FIGURES/QC/Figure_03_Distribution_of_Mitochondrial_Gene_Percentage.tif"
  
)

#----------------------------------------------------------
# 12. Figure 4
#----------------------------------------------------------

counts_vs_genes <- create_scatter_plot(
  
  feature1 = "nCount_RNA",
  
  feature2 = "nFeature_RNA",
  
  title = "Relationship Between UMI Counts and Detected Genes",
  
  xlabel = "UMI Counts",
  
  ylabel = "Detected Genes",
  
  filename = "FIGURES/QC/Figure_04_UMI_Counts_vs_Detected_Genes.tif"
  
)

#----------------------------------------------------------
# 13. Figure 5
#----------------------------------------------------------

counts_vs_mito <- create_scatter_plot(
  
  feature1 = "nCount_RNA",
  
  feature2 = "percent.mt",
  
  title = "Relationship Between UMI Counts and Mitochondrial Gene Percentage",
  
  xlabel = "UMI Counts",
  
  ylabel = "Mitochondrial Gene Percentage (%)",
  
  filename = "FIGURES/QC/Figure_05_UMI_Counts_vs_Mitochondrial_Gene_Percentage.tif"
  
)

#----------------------------------------------------------
# 14. Figure 6
#----------------------------------------------------------

genes_vs_mito <- create_scatter_plot(
  
  feature1 = "nFeature_RNA",
  
  feature2 = "percent.mt",
  
  title = "Relationship Between Detected Genes and Mitochondrial Gene Percentage",
  
  xlabel = "Detected Genes",
  
  ylabel = "Mitochondrial Gene Percentage (%)",
  
  filename = "FIGURES/QC/Figure_06_Detected_Genes_vs_Mitochondrial_Gene_Percentage.tif"
  
)

#----------------------------------------------------------
# 15. Display Plots
#----------------------------------------------------------

print(genes_plot)
print(umi_plot)
print(mito_plot)

print(counts_vs_genes)
print(counts_vs_mito)
print(genes_vs_mito)

cat("\nQC figures generated successfully.\n")
cat("Location : FIGURES/QC/\n")