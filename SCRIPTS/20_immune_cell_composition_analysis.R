###############################################################
# SCRIPT 20: IMMUNE CELL COMPOSITION ANALYSIS
###############################################################

rm(list = ls())
gc()

###############################################################
# Load Libraries
###############################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

###############################################################
# Create Output Folders
###############################################################

dir.create(
  "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS/Cell_Composition",
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
# Create Disease Groups
###############################################################

scaled_seurat$Condition <- case_when(
  grepl("^H", scaled_seurat$orig.ident) ~ "Healthy",
  grepl("^A", scaled_seurat$orig.ident) ~ "AChR_MG",
  grepl("a$", scaled_seurat$orig.ident) ~ "SNMG_Pre",
  grepl("b$", scaled_seurat$orig.ident) ~ "SNMG_Post"
)

###############################################################
# Extract Metadata
###############################################################

metadata <- scaled_seurat@meta.data %>%
  select(
    Sample = orig.ident,
    Condition,
    Cell_Type = SingleR_Label
  )

###############################################################
# Cell Counts Per Sample
###############################################################

cell.counts <- metadata %>%
  count(
    Sample,
    Condition,
    Cell_Type
  )

###############################################################
# Cell Percentages
###############################################################

cell.percentages <- cell.counts %>%
  group_by(Sample) %>%
  mutate(
    Percentage = n / sum(n) * 100
  ) %>%
  ungroup()

###############################################################
# Save Tables
###############################################################

write.csv(
  cell.counts,
  "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS/Cell_Composition/Cell_Counts_Per_Sample.csv",
  row.names = FALSE
)

write.csv(
  cell.percentages,
  "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS/Cell_Composition/Cell_Percentages_Per_Sample.csv",
  row.names = FALSE
)

###############################################################
# Stacked Bar Plot
###############################################################

p1 <- ggplot(
  cell.percentages,
  aes(
    x = Sample,
    y = Percentage,
    fill = Cell_Type
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.9
  ) +
  facet_grid(
    ~Condition,
    scales = "free_x",
    space = "free_x"
  ) +
  labs(
    title = "Peripheral Immune Cell Composition",
    x = "Sample",
    y = "Percentage (%)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

ggsave(
  filename = "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS/Cell_Composition/Cell_Composition_StackedBar.png",
  plot = p1,
  width = 14,
  height = 7,
  dpi = 300
)

###############################################################
# Boxplots
###############################################################

p2 <- ggplot(
  cell.percentages,
  aes(
    x = Condition,
    y = Percentage,
    fill = Condition
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    size = 1.5,
    alpha = 0.8
  ) +
  facet_wrap(
    ~Cell_Type,
    scales = "free_y"
  ) +
  labs(
    title = "Immune Cell Composition Across Disease Groups",
    y = "Percentage (%)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "none"
  )

ggsave(
  filename = "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS/Cell_Composition/Cell_Composition_Boxplots.png",
  plot = p2,
  width = 15,
  height = 10,
  dpi = 300
)

###############################################################
# Kruskal-Wallis Test
###############################################################

statistics <- cell.percentages %>%
  group_by(Cell_Type) %>%
  summarise(
    P_Value = kruskal.test(
      Percentage ~ Condition
    )$p.value
  ) %>%
  mutate(
    FDR = p.adjust(
      P_Value,
      method = "BH"
    )
  ) %>%
  arrange(FDR)

write.csv(
  statistics,
  "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS/Cell_Composition/Cell_Composition_Statistics.csv",
  row.names = FALSE
)

###############################################################
# Summary Table
###############################################################

summary.table <- cell.percentages %>%
  group_by(
    Condition,
    Cell_Type
  ) %>%
  summarise(
    Mean = mean(Percentage),
    SD = sd(Percentage),
    .groups = "drop"
  )

write.csv(
  summary.table,
  "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS/Cell_Composition/Cell_Composition_Summary.csv",
  row.names = FALSE
)

###############################################################
# Finished
###############################################################

cat("\n=========================================\n")
cat("Immune Cell Composition Analysis Complete\n")
cat("=========================================\n")