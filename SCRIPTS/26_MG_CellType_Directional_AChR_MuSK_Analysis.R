###############################################################################
# SCRIPT 26 : MG CELL-TYPE DIRECTIONAL ACTIVITY + AChR / MuSK TARGETED ANALYSIS
#
# PURPOSE
#   1. Identify which immune cell types show predominantly upregulated versus
#      downregulated significant transcriptional changes in each MG comparison.
#   2. Quantify DEG burden, direction balance, and effect-size magnitude.
#   3. Generate manuscript-ready comparison-specific figures and tables.
#   4. Examine AChR-related receptor genes and MUSK in the available DEG data
#      WITHOUT incorrectly equating antibody-positive MG with gene overexpression.
#   5. Produce data-driven Results text and interpretation metadata.
#
# IMPORTANT INTERPRETATION
#   "AChR-positive" and "MuSK-positive" MG refer to autoantibody status.
#   They do NOT automatically mean CHRNA1/MUSK gene overexpression.
#
# RUN AFTER
#   The integrated analysis has created `master_deg`.
#
# REQUIRED master_deg COLUMNS
#   Gene, CellType, Comparison, avg_log2FC, p_val_adj
###############################################################################

options(stringsAsFactors = FALSE)

###############################################################################
# 0. SETUP
###############################################################################

required_objects <- c(
  "master_deg",
  "TABLE_DIR",
  "FIGURE_DIR",
  "REPORT_DIR",
  "ALPHA",
  "LOGFC_THRESHOLD"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    "Script 26 missing required objects: ",
    paste(missing_objects, collapse = ", "),
    "\nRun the integrated DEG analysis first."
  )
}

required_columns <- c(
  "Gene",
  "CellType",
  "Comparison",
  "avg_log2FC",
  "p_val_adj"
)

missing_columns <- setdiff(required_columns, colnames(master_deg))

if (length(missing_columns) > 0) {
  stop(
    "master_deg is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

pkg_needed <- c("dplyr", "tidyr", "ggplot2")
pkg_missing <- pkg_needed[
  !vapply(pkg_needed, requireNamespace, logical(1), quietly = TRUE)
]

if (length(pkg_missing) > 0) {
  stop(
    "Install required packages first: ",
    paste(pkg_missing, collapse = ", ")
  )
}

SCRIPT26_TABLE_DIR <- file.path(TABLE_DIR, "Script26_MG_Directional_Analysis")
SCRIPT26_FIGURE_DIR <- file.path(FIGURE_DIR, "Script26_MG_Directional_Analysis")
SCRIPT26_REPORT_DIR <- file.path(REPORT_DIR, "Script26_MG_Directional_Analysis")

dir.create(SCRIPT26_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SCRIPT26_FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SCRIPT26_REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

save26 <- function(x, filename) {
  utils::write.csv(
    x,
    file.path(SCRIPT26_TABLE_DIR, filename),
    row.names = FALSE
  )
}

clean_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]+", "_", x)
  gsub("_+", "_", x)
}

message("\n==============================================================")
message(" SCRIPT 26 : MG CELL-TYPE DIRECTIONAL ACTIVITY ANALYSIS")
message("==============================================================\n")

###############################################################################
# 1. CREATE SIGNIFICANT DEG DATASET
###############################################################################

significant_deg <- master_deg |>
  dplyr::filter(
    !is.na(Gene),
    Gene != "",
    !is.na(CellType),
    CellType != "",
    !is.na(Comparison),
    Comparison != "",
    !is.na(avg_log2FC),
    !is.na(p_val_adj),
    p_val_adj < ALPHA,
    abs(avg_log2FC) >= LOGFC_THRESHOLD
  ) |>
  dplyr::mutate(
    Direction = dplyr::case_when(
      avg_log2FC > 0 ~ "Upregulated",
      avg_log2FC < 0 ~ "Downregulated",
      TRUE ~ "No_change"
    )
  )

if (nrow(significant_deg) == 0) {
  stop(
    "No significant DEGs remain after filtering at adjusted P < ",
    ALPHA,
    " and |avg_log2FC| >= ",
    LOGFC_THRESHOLD
  )
}

cell_types <- sort(unique(significant_deg$CellType))
comparisons <- sort(unique(significant_deg$Comparison))

message("Significant DEG entries : ", nrow(significant_deg))
message("Unique significant genes: ", dplyr::n_distinct(significant_deg$Gene))
message("Cell types              : ", length(cell_types))
message("Comparisons             : ", length(comparisons))
message("\nCell types:")
message(paste(cell_types, collapse = " | "))
message("\nComparisons:")
message(paste(comparisons, collapse = " | "))

save26(
  significant_deg,
  "26.1_All_Significant_DEG_Entries.csv"
)

###############################################################################
# 2. CELL TYPE × COMPARISON DIRECTIONAL SUMMARY
###############################################################################

direction_counts <- significant_deg |>
  dplyr::count(
    CellType,
    Comparison,
    Direction,
    name = "DEG_Count"
  ) |>
  tidyr::complete(
    CellType = cell_types,
    Comparison = comparisons,
    Direction = c("Upregulated", "Downregulated"),
    fill = list(DEG_Count = 0L)
  ) |>
  tidyr::pivot_wider(
    names_from = Direction,
    values_from = DEG_Count,
    values_fill = 0
  ) |>
  dplyr::mutate(
    Total_Significant_DEGs = Upregulated + Downregulated,
    Percent_Up = dplyr::if_else(
      Total_Significant_DEGs > 0,
      100 * Upregulated / Total_Significant_DEGs,
      0
    ),
    Percent_Down = dplyr::if_else(
      Total_Significant_DEGs > 0,
      100 * Downregulated / Total_Significant_DEGs,
      0
    ),
    Net_Direction_Score = dplyr::if_else(
      Total_Significant_DEGs > 0,
      (Upregulated - Downregulated) / Total_Significant_DEGs,
      0
    ),
    Direction_Bias = dplyr::case_when(
      Net_Direction_Score >= 0.20 ~ "Predominantly upregulated",
      Net_Direction_Score <= -0.20 ~ "Predominantly downregulated",
      TRUE ~ "Mixed / balanced"
    )
  )

###############################################################################
# 3. EFFECT-SIZE SUMMARY
###############################################################################

effect_summary <- significant_deg |>
  dplyr::group_by(CellType, Comparison) |>
  dplyr::summarise(
    Mean_Log2FC = mean(avg_log2FC, na.rm = TRUE),
    Median_Log2FC = median(avg_log2FC, na.rm = TRUE),
    Mean_Abs_Log2FC = mean(abs(avg_log2FC), na.rm = TRUE),

    Mean_Up_Log2FC = if (
      any(avg_log2FC > 0, na.rm = TRUE)
    ) {
      mean(avg_log2FC[avg_log2FC > 0], na.rm = TRUE)
    } else {
      NA_real_
    },

    Mean_Down_Log2FC = if (
      any(avg_log2FC < 0, na.rm = TRUE)
    ) {
      mean(avg_log2FC[avg_log2FC < 0], na.rm = TRUE)
    } else {
      NA_real_
    },

    Max_Up_Log2FC = if (
      any(avg_log2FC > 0, na.rm = TRUE)
    ) {
      max(avg_log2FC[avg_log2FC > 0], na.rm = TRUE)
    } else {
      NA_real_
    },

    Strongest_Down_Log2FC = if (
      any(avg_log2FC < 0, na.rm = TRUE)
    ) {
      min(avg_log2FC[avg_log2FC < 0], na.rm = TRUE)
    } else {
      NA_real_
    },

    .groups = "drop"
  )

cell_state_summary <- direction_counts |>
  dplyr::left_join(
    effect_summary,
    by = c("CellType", "Comparison")
  ) |>
  dplyr::mutate(
    # This is a descriptive transcriptomic perturbation score.
    # It is NOT a direct measurement of physiological cell activation.
    Transcriptional_Perturbation_Score =
      log1p(Total_Significant_DEGs) * Mean_Abs_Log2FC,

    Upregulation_Score =
      log1p(Upregulated) *
      dplyr::if_else(
        is.na(Mean_Up_Log2FC),
        0,
        Mean_Up_Log2FC
      ),

    Downregulation_Score =
      log1p(Downregulated) *
      abs(
        dplyr::if_else(
          is.na(Mean_Down_Log2FC),
          0,
          Mean_Down_Log2FC
        )
      )
  )

save26(
  cell_state_summary,
  "26.2_CellType_Comparison_Directional_Summary.csv"
)

###############################################################################
# 4. RANK CELL TYPES WITHIN EACH COMPARISON
###############################################################################

ranked_cells <- cell_state_summary |>
  dplyr::group_by(Comparison) |>
  dplyr::mutate(
    Rank_Total_Perturbation =
      dplyr::min_rank(dplyr::desc(Transcriptional_Perturbation_Score)),
    Rank_Upregulation =
      dplyr::min_rank(dplyr::desc(Upregulation_Score)),
    Rank_Downregulation =
      dplyr::min_rank(dplyr::desc(Downregulation_Score))
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(
    Comparison,
    Rank_Total_Perturbation
  )

save26(
  ranked_cells,
  "26.3_CellType_Rankings_by_Comparison.csv"
)

###############################################################################
# 5. TOP UPREGULATED AND DOWNREGULATED GENES PER CELL TYPE
###############################################################################

top_up_genes <- significant_deg |>
  dplyr::filter(Direction == "Upregulated") |>
  dplyr::group_by(Comparison, CellType) |>
  dplyr::arrange(
    dplyr::desc(avg_log2FC),
    p_val_adj,
    .by_group = TRUE
  ) |>
  dplyr::slice_head(n = 15) |>
  dplyr::ungroup()

top_down_genes <- significant_deg |>
  dplyr::filter(Direction == "Downregulated") |>
  dplyr::group_by(Comparison, CellType) |>
  dplyr::arrange(
    avg_log2FC,
    p_val_adj,
    .by_group = TRUE
  ) |>
  dplyr::slice_head(n = 15) |>
  dplyr::ungroup()

save26(
  top_up_genes,
  "26.4_Top15_Upregulated_Genes_per_CellType_Comparison.csv"
)

save26(
  top_down_genes,
  "26.4_Top15_Downregulated_Genes_per_CellType_Comparison.csv"
)

###############################################################################
# 6. FIGURE A : UPREGULATED VS DOWNREGULATED DEG COUNTS
###############################################################################

plot_direction <- significant_deg |>
  dplyr::count(
    Comparison,
    CellType,
    Direction,
    name = "DEG_Count"
  ) |>
  dplyr::mutate(
    Plot_Count = dplyr::if_else(
      Direction == "Downregulated",
      -DEG_Count,
      DEG_Count
    )
  )

p_direction <- ggplot2::ggplot(
  plot_direction,
  ggplot2::aes(
    x = CellType,
    y = Plot_Count,
    fill = Direction
  )
) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.5
  ) +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(
    ~ Comparison,
    scales = "free_y"
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Upregulated" = "#C43C39",
      "Downregulated" = "#3B6FB6"
    )
  ) +
  ggplot2::scale_y_continuous(
    labels = abs
  ) +
  ggplot2::labs(
    title =
      "Direction of Significant Transcriptional Changes Across Immune Cell Types",
    subtitle =
      "Bars to the right represent upregulated DEGs; bars to the left represent downregulated DEGs",
    x = NULL,
    y = "Number of significant DEGs",
    fill = "Direction",
    caption = paste0(
      "Adjusted P < ", ALPHA,
      "; |average log2FC| >= ", LOGFC_THRESHOLD,
      ". DEG direction is not equivalent to direct physiological activation."
    )
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 14
      ),
    legend.position = "bottom",
    strip.text =
      ggplot2::element_text(
        face = "bold"
      )
  )

ggplot2::ggsave(
  file.path(
    SCRIPT26_FIGURE_DIR,
    "26_Figure_A_Up_vs_Down_by_CellType.png"
  ),
  p_direction,
  width = 14,
  height = 8.5,
  dpi = 320,
  bg = "white"
)

ggplot2::ggsave(
  file.path(
    SCRIPT26_FIGURE_DIR,
    "26_Figure_A_Up_vs_Down_by_CellType.pdf"
  ),
  p_direction,
  width = 14,
  height = 8.5
)

###############################################################################
# 7. FIGURE B : NET DIRECTION HEATMAP
###############################################################################

direction_heat <- cell_state_summary |>
  dplyr::select(
    CellType,
    Comparison,
    Net_Direction_Score
  ) |>
  tidyr::pivot_wider(
    names_from = Comparison,
    values_from = Net_Direction_Score,
    values_fill = 0
  )

direction_mat <- direction_heat |>
  dplyr::select(-CellType) |>
  as.matrix()

rownames(direction_mat) <- direction_heat$CellType
storage.mode(direction_mat) <- "numeric"

if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("circlize", quietly = TRUE)) {

  direction_col <- circlize::colorRamp2(
    c(-1, 0, 1),
    c("#2166AC", "#F7F7F7", "#B2182B")
  )

  direction_ht <- ComplexHeatmap::Heatmap(
    direction_mat,
    name = "Net direction",
    col = direction_col,
    cluster_rows = nrow(direction_mat) > 1,
    cluster_columns = ncol(direction_mat) > 1,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp =
      grid::gpar(
        fontsize = 9,
        fontface = "bold"
      ),
    column_names_gp =
      grid::gpar(
        fontsize = 9,
        fontface = "bold"
      ),
    column_names_rot = 35,
    row_title = "Immune cell types",
    column_title =
      "MG-associated transcriptional direction by comparison",
    heatmap_legend_param = list(
      title = "Net direction",
      at = c(-1, 0, 1),
      labels = c(
        "Downregulated",
        "Balanced",
        "Upregulated"
      )
    ),
    border = TRUE
  )

  grDevices::png(
    file.path(
      SCRIPT26_FIGURE_DIR,
      "26_Figure_B_Net_Direction_Heatmap.png"
    ),
    width = 3000,
    height = 2500,
    res = 320
  )

  ComplexHeatmap::draw(
    direction_ht,
    heatmap_legend_side = "right"
  )

  grDevices::dev.off()

  grDevices::pdf(
    file.path(
      SCRIPT26_FIGURE_DIR,
      "26_Figure_B_Net_Direction_Heatmap.pdf"
    ),
    width = 10,
    height = 8
  )

  ComplexHeatmap::draw(
    direction_ht,
    heatmap_legend_side = "right"
  )

  grDevices::dev.off()
}

###############################################################################
# 8. FIGURE C : TRANSCRIPTIONAL PERTURBATION
###############################################################################

p_perturbation <- ggplot2::ggplot(
  cell_state_summary,
  ggplot2::aes(
    x = stats::reorder(
      CellType,
      Transcriptional_Perturbation_Score
    ),
    y = Transcriptional_Perturbation_Score,
    fill = Direction_Bias
  )
) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(
    ~ Comparison,
    scales = "free_y"
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Predominantly upregulated" = "#C43C39",
      "Mixed / balanced" = "#A7A7A7",
      "Predominantly downregulated" = "#3B6FB6"
    )
  ) +
  ggplot2::labs(
    title =
      "Relative Transcriptional Perturbation Across Immune Cell Types",
    subtitle =
      "Composite descriptive score combining significant-DEG burden and mean absolute effect size",
    x = NULL,
    y = "Transcriptional perturbation score",
    fill = "Direction bias",
    caption =
      "This score ranks transcriptional perturbation; it is not a direct assay of cellular activation."
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 14
      ),
    legend.position = "bottom",
    strip.text =
      ggplot2::element_text(
        face = "bold"
      )
  )

ggplot2::ggsave(
  file.path(
    SCRIPT26_FIGURE_DIR,
    "26_Figure_C_Transcriptional_Perturbation.png"
  ),
  p_perturbation,
  width = 14,
  height = 8.5,
  dpi = 320,
  bg = "white"
)

ggplot2::ggsave(
  file.path(
    SCRIPT26_FIGURE_DIR,
    "26_Figure_C_Transcriptional_Perturbation.pdf"
  ),
  p_perturbation,
  width = 14,
  height = 8.5
)

###############################################################################
# 9. TARGETED AChR / MuSK GENE ANALYSIS
###############################################################################
#
# Autoantibody subtype != receptor gene overexpression.
#
# AChR at the neuromuscular junction is composed of muscle nicotinic receptor
# subunits. We therefore search a defined receptor-related panel rather than
# treating the text "AChR+" as expression of one gene.
###############################################################################

AChR_related_genes <- c(
  "CHRNA1",
  "CHRNB1",
  "CHRND",
  "CHRNE",
  "CHRNG",
  "RAPSN",
  "AGRN",
  "LRP4",
  "DOK7"
)

MuSK_related_genes <- c(
  "MUSK",
  "LRP4",
  "AGRN",
  "DOK7",
  "RAPSN"
)

nmj_target_panel <- unique(
  c(
    AChR_related_genes,
    MuSK_related_genes
  )
)

available_master_genes <- unique(
  toupper(master_deg$Gene)
)

target_availability <- data.frame(
  Gene = nmj_target_panel,
  Target_Group = dplyr::case_when(
    nmj_target_panel == "MUSK" ~ "MuSK",
    nmj_target_panel %in% c(
      "CHRNA1",
      "CHRNB1",
      "CHRND",
      "CHRNE",
      "CHRNG"
    ) ~ "AChR receptor subunit",
    TRUE ~ "NMJ-associated"
  ),
  Present_in_Master_DEG_Table =
    nmj_target_panel %in% available_master_genes,
  stringsAsFactors = FALSE
)

save26(
  target_availability,
  "26.5_AChR_MuSK_Target_Gene_Availability.csv"
)

target_deg <- master_deg |>
  dplyr::mutate(
    Gene_Upper = toupper(Gene)
  ) |>
  dplyr::filter(
    Gene_Upper %in% nmj_target_panel
  ) |>
  dplyr::mutate(
    Target_Group = dplyr::case_when(
      Gene_Upper == "MUSK" ~ "MuSK",
      Gene_Upper %in% c(
        "CHRNA1",
        "CHRNB1",
        "CHRND",
        "CHRNE",
        "CHRNG"
      ) ~ "AChR receptor subunit",
      TRUE ~ "NMJ-associated"
    ),
    Significant = !is.na(p_val_adj) &
      p_val_adj < ALPHA &
      !is.na(avg_log2FC) &
      abs(avg_log2FC) >= LOGFC_THRESHOLD,
    Expression_Direction = dplyr::case_when(
      Significant & avg_log2FC > 0 ~ "Significantly upregulated",
      Significant & avg_log2FC < 0 ~ "Significantly downregulated",
      TRUE ~ "Not significant by Script 26 thresholds"
    )
  )

save26(
  target_deg,
  "26.6_AChR_MuSK_Target_DEG_Results.csv"
)

###############################################################################
# 10. TARGET-GENE FIGURE, ONLY IF TARGETS ARE ACTUALLY PRESENT
###############################################################################

if (nrow(target_deg) > 0) {

  target_plot_df <- target_deg |>
    dplyr::filter(
      !is.na(avg_log2FC)
    )

  if (nrow(target_plot_df) > 0) {

    p_targets <- ggplot2::ggplot(
      target_plot_df,
      ggplot2::aes(
        x = CellType,
        y = avg_log2FC,
        fill = Expression_Direction
      )
    ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linewidth = 0.45
      ) +
      ggplot2::geom_point(
        shape = 21,
        size = 3,
        alpha = 0.85
      ) +
      ggplot2::facet_grid(
        Gene_Upper ~ Comparison,
        scales = "free_y"
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(
        values = c(
          "Significantly upregulated" = "#C43C39",
          "Significantly downregulated" = "#3B6FB6",
          "Not significant by Script 26 thresholds" = "#C8C8C8"
        )
      ) +
      ggplot2::labs(
        title =
          "AChR / MuSK / Neuromuscular-Junction Target Genes in the DEG Dataset",
        subtitle =
          "This tests transcript-level evidence; antibody-positive status is not interpreted as gene overexpression",
        x = NULL,
        y = "Average log2 fold-change",
        fill = "DEG status"
      ) +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(
        plot.title =
          ggplot2::element_text(
            face = "bold"
          ),
        legend.position = "bottom",
        strip.text =
          ggplot2::element_text(
            face = "bold",
            size = 8
          )
      )

    ggplot2::ggsave(
      file.path(
        SCRIPT26_FIGURE_DIR,
        "26_Figure_D_AChR_MuSK_Target_Genes.png"
      ),
      p_targets,
      width = 15,
      height = max(
        7,
        2.2 * dplyr::n_distinct(target_plot_df$Gene_Upper)
      ),
      dpi = 320,
      bg = "white",
      limitsize = FALSE
    )

    ggplot2::ggsave(
      file.path(
        SCRIPT26_FIGURE_DIR,
        "26_Figure_D_AChR_MuSK_Target_Genes.pdf"
      ),
      p_targets,
      width = 15,
      height = max(
        7,
        2.2 * dplyr::n_distinct(target_plot_df$Gene_Upper)
      ),
      limitsize = FALSE
    )
  }
}

###############################################################################
# 11. AUTOMATIC DATA-DRIVEN INTERPRETATION
###############################################################################

top_total <- ranked_cells |>
  dplyr::filter(Rank_Total_Perturbation == 1) |>
  dplyr::select(
    Comparison,
    CellType,
    Total_Significant_DEGs,
    Transcriptional_Perturbation_Score,
    Direction_Bias
  )

top_up <- ranked_cells |>
  dplyr::filter(Rank_Upregulation == 1) |>
  dplyr::select(
    Comparison,
    CellType,
    Upregulated,
    Percent_Up,
    Mean_Up_Log2FC,
    Upregulation_Score
  )

top_down <- ranked_cells |>
  dplyr::filter(Rank_Downregulation == 1) |>
  dplyr::select(
    Comparison,
    CellType,
    Downregulated,
    Percent_Down,
    Mean_Down_Log2FC,
    Downregulation_Score
  )

save26(
  top_total,
  "26.7_Most_Transcriptionally_Perturbed_CellType.csv"
)

save26(
  top_up,
  "26.7_Strongest_Upregulated_CellType.csv"
)

save26(
  top_down,
  "26.7_Strongest_Downregulated_CellType.csv"
)

interpretation_lines <- c(
  "SCRIPT 26 : DATA-DRIVEN INTERPRETATION",
  "======================================",
  "",
  paste0(
    "Significant DEG definition: adjusted P < ",
    ALPHA,
    " and |average log2FC| >= ",
    LOGFC_THRESHOLD,
    "."
  ),
  "",
  "IMPORTANT:",
  paste0(
    "The terms 'more active' and 'less active' are not used as direct ",
    "physiological claims. This script measures transcriptional perturbation ",
    "and directional DEG bias. Functional activation requires pathway, ",
    "protein, signalling or experimental evidence."
  ),
  "",
  paste0(
    "Analysed cell types: ",
    paste(cell_types, collapse = ", "),
    "."
  ),
  "",
  paste0(
    "Analysed comparisons: ",
    paste(comparisons, collapse = ", "),
    "."
  ),
  "",
  "MOST TRANSCRIPTIONALLY PERTURBED CELL TYPE BY COMPARISON",
  "--------------------------------------------------------"
)

for (i in seq_len(nrow(top_total))) {
  interpretation_lines <- c(
    interpretation_lines,
    paste0(
      top_total$Comparison[i],
      ": ",
      top_total$CellType[i],
      " ranked highest for the composite transcriptional perturbation score, ",
      "with ",
      top_total$Total_Significant_DEGs[i],
      " significant DEG entries and a ",
      tolower(top_total$Direction_Bias[i]),
      " profile."
    )
  )
}

interpretation_lines <- c(
  interpretation_lines,
  "",
  "STRONGEST UPREGULATED PROFILE BY COMPARISON",
  "-------------------------------------------"
)

for (i in seq_len(nrow(top_up))) {
  interpretation_lines <- c(
    interpretation_lines,
    paste0(
      top_up$Comparison[i],
      ": ",
      top_up$CellType[i],
      " ranked highest for the upregulation score, with ",
      top_up$Upregulated[i],
      " upregulated DEG entries (",
      round(top_up$Percent_Up[i], 1),
      "% of its significant DEG entries)."
    )
  )
}

interpretation_lines <- c(
  interpretation_lines,
  "",
  "STRONGEST DOWNREGULATED PROFILE BY COMPARISON",
  "---------------------------------------------"
)

for (i in seq_len(nrow(top_down))) {
  interpretation_lines <- c(
    interpretation_lines,
    paste0(
      top_down$Comparison[i],
      ": ",
      top_down$CellType[i],
      " ranked highest for the downregulation score, with ",
      top_down$Downregulated[i],
      " downregulated DEG entries (",
      round(top_down$Percent_Down[i], 1),
      "% of its significant DEG entries)."
    )
  )
}

###############################################################################
# 12. AChR / MuSK INTERPRETATION WITHOUT OVERCLAIMING
###############################################################################

achr_present <- target_availability |>
  dplyr::filter(
    Target_Group == "AChR receptor subunit",
    Present_in_Master_DEG_Table
  ) |>
  dplyr::pull(Gene)

musk_present <- target_availability |>
  dplyr::filter(
    Gene == "MUSK",
    Present_in_Master_DEG_Table
  ) |>
  nrow() > 0

sig_achr <- target_deg |>
  dplyr::filter(
    Target_Group == "AChR receptor subunit",
    Significant
  )

sig_musk <- target_deg |>
  dplyr::filter(
    Gene_Upper == "MUSK",
    Significant
  )

interpretation_lines <- c(
  interpretation_lines,
  "",
  "AChR / MuSK TARGETED ANALYSIS",
  "-----------------------------",
  paste0(
    "AChR-related receptor subunits detected anywhere in the master DEG table: ",
    ifelse(
      length(achr_present) > 0,
      paste(achr_present, collapse = ", "),
      "none"
    ),
    "."
  ),
  paste0(
    "MUSK detected anywhere in the master DEG table: ",
    ifelse(musk_present, "yes", "no"),
    "."
  ),
  paste0(
    "Significant AChR receptor-subunit DEG records passing Script 26 thresholds: ",
    nrow(sig_achr),
    "."
  ),
  paste0(
    "Significant MUSK DEG records passing Script 26 thresholds: ",
    nrow(sig_musk),
    "."
  ),
  "",
  paste0(
    "AChR-positive and MuSK-positive myasthenia gravis are antibody-defined ",
    "disease subtypes. Detection of an AChR-related or MUSK transcript as ",
    "upregulated would represent transcript-level differential expression in ",
    "the sampled cells; it should not be used as evidence that antibody-positive ",
    "MG is caused by receptor overexpression."
  )
)

writeLines(
  interpretation_lines,
  file.path(
    SCRIPT26_REPORT_DIR,
    "26_Data_Driven_Interpretation.txt"
  )
)

###############################################################################
# 13. FIGURE LEGENDS
###############################################################################

legend_A <- paste0(
  "Figure 26A. Direction of significant transcriptional changes across immune ",
  "cell types and MG-related comparisons. Significant DEGs were defined using ",
  "adjusted P < ", ALPHA, " and |average log2 fold-change| >= ",
  LOGFC_THRESHOLD, ". Positive/rightward bars indicate the number of ",
  "upregulated DEGs and negative/leftward bars indicate the number of ",
  "downregulated DEGs. Panels correspond to individual comparisons. DEG ",
  "direction represents transcript-level change and should not alone be ",
  "interpreted as direct physiological activation or inhibition."
)

legend_B <- paste0(
  "Figure 26B. Net transcriptional direction across immune cell types. Each ",
  "cell represents (number of upregulated DEGs - number of downregulated DEGs) ",
  "/ total significant DEGs for a given cell type and comparison. Values ",
  "approaching +1 indicate predominantly upregulated significant changes, ",
  "values approaching -1 indicate predominantly downregulated changes, and ",
  "values near zero indicate a mixed directional profile."
)

legend_C <- paste0(
  "Figure 26C. Relative transcriptional perturbation across immune cell types. ",
  "The descriptive perturbation score combines significant-DEG burden with ",
  "mean absolute log2 fold-change [log(1 + number of significant DEGs) x mean ",
  "|log2FC|]. Bar colour denotes directional bias. This metric is intended for ",
  "within-analysis ranking and is not a direct experimental measurement of ",
  "immune-cell activation."
)

legend_D <- paste0(
  "Figure 26D. Targeted examination of acetylcholine-receptor (AChR)-related, ",
  "MuSK and neuromuscular-junction-associated genes represented in the DEG ",
  "dataset. Points show average log2 fold-change by immune cell type and ",
  "comparison, with colour indicating whether the transcript passes the ",
  "Script 26 DEG thresholds. AChR-positive or MuSK-positive clinical status ",
  "denotes autoantibody status and is not interpreted as receptor-gene ",
  "overexpression."
)

figure_legends <- data.frame(
  Figure = c(
    "26A",
    "26B",
    "26C",
    "26D"
  ),
  Legend = c(
    legend_A,
    legend_B,
    legend_C,
    legend_D
  ),
  stringsAsFactors = FALSE
)

save26(
  figure_legends,
  "26.8_Figure_Legends.csv"
)

writeLines(
  c(
    "SCRIPT 26 FIGURE LEGENDS",
    "========================",
    "",
    "Figure 26A",
    legend_A,
    "",
    "Figure 26B",
    legend_B,
    "",
    "Figure 26C",
    legend_C,
    "",
    "Figure 26D",
    legend_D
  ),
  file.path(
    SCRIPT26_REPORT_DIR,
    "26_Figure_Legends.txt"
  )
)

###############################################################################
# 14. ANALYSIS METADATA
###############################################################################

metadata26 <- data.frame(
  Field = c(
    "Analysis",
    "Adjusted P threshold",
    "Absolute average log2FC threshold",
    "Significant DEG entries",
    "Unique significant genes",
    "Number of cell types",
    "Cell types",
    "Number of comparisons",
    "Comparisons",
    "Net direction score",
    "Perturbation score",
    "AChR interpretation",
    "MuSK interpretation"
  ),
  Value = c(
    "MG cell-type directional transcriptional analysis",
    ALPHA,
    LOGFC_THRESHOLD,
    nrow(significant_deg),
    dplyr::n_distinct(significant_deg$Gene),
    length(cell_types),
    paste(cell_types, collapse = "; "),
    length(comparisons),
    paste(comparisons, collapse = "; "),
    "(Upregulated - Downregulated) / Total significant DEGs",
    "log(1 + significant DEG count) x mean absolute log2FC",
    "AChR-positive status is treated as autoantibody status, not receptor overexpression",
    "MuSK-positive status is treated as autoantibody status, not MUSK overexpression"
  ),
  stringsAsFactors = FALSE
)

save26(
  metadata26,
  "26.9_Analysis_Metadata.csv"
)

###############################################################################
# 15. FINAL MANIFEST
###############################################################################

manifest26 <- data.frame(
  Output = c(
    "26.1_All_Significant_DEG_Entries.csv",
    "26.2_CellType_Comparison_Directional_Summary.csv",
    "26.3_CellType_Rankings_by_Comparison.csv",
    "26.4_Top15_Upregulated_Genes_per_CellType_Comparison.csv",
    "26.4_Top15_Downregulated_Genes_per_CellType_Comparison.csv",
    "26.5_AChR_MuSK_Target_Gene_Availability.csv",
    "26.6_AChR_MuSK_Target_DEG_Results.csv",
    "26.7_Most_Transcriptionally_Perturbed_CellType.csv",
    "26.7_Strongest_Upregulated_CellType.csv",
    "26.7_Strongest_Downregulated_CellType.csv",
    "26.8_Figure_Legends.csv",
    "26.9_Analysis_Metadata.csv"
  ),
  Directory = SCRIPT26_TABLE_DIR,
  stringsAsFactors = FALSE
)

save26(
  manifest26,
  "26_Output_Manifest.csv"
)

###############################################################################
# 16. FINISH
###############################################################################

message("\n==============================================================")
message(" SCRIPT 26 COMPLETED")
message("==============================================================")
message("Tables : ", SCRIPT26_TABLE_DIR)
message("Figures: ", SCRIPT26_FIGURE_DIR)
message("Reports: ", SCRIPT26_REPORT_DIR)
message("")
message("Key files to inspect first:")
message("  26.2_CellType_Comparison_Directional_Summary.csv")
message("  26.3_CellType_Rankings_by_Comparison.csv")
message("  26_Figure_A_Up_vs_Down_by_CellType.png")
message("  26_Figure_B_Net_Direction_Heatmap.png")
message("  26_Figure_C_Transcriptional_Perturbation.png")
message("  26.6_AChR_MuSK_Target_DEG_Results.csv")
message("  26_Data_Driven_Interpretation.txt")
message("==============================================================\n")
