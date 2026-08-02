###############################################################################
# 24_5_REPORTING_REPAIR.R
#
# Purpose:
#   Repair and replace the downstream reporting/visualisation layer of the
#   integrated scRNA-seq DEG analysis (Sections 5.4 onward).
#
# IMPORTANT:
#   Run AFTER Sections 5.1-5.3 have successfully created/loaded `master_deg`
#   and defined TABLE_DIR, FIGURE_DIR, REPORT_DIR, ALPHA, LOGFC_THRESHOLD.
#
# Repairs:
#   - recurrence analysis uses SIGNIFICANT DEGs, not every tested gene
#   - readable recurrence heatmap and recurrence-class distribution
#   - robust conserved-signature scoring (no .x/.y join bug)
#   - detailed metadata tables
#   - figure-specific, data-driven legends
#   - quantitative Results text with figure references
#   - publication-oriented methodology paragraph
#   - coloured methodology flowchart
#   - supplementary tables
###############################################################################

options(stringsAsFactors = FALSE)

###############################################################################
# 0. SAFETY / COMPATIBILITY SETUP
###############################################################################

required_objects <- c("master_deg", "TABLE_DIR", "FIGURE_DIR", "REPORT_DIR",
                      "ALPHA", "LOGFC_THRESHOLD")
missing_objects <- required_objects[!vapply(required_objects, exists,
                                            logical(1), inherits = TRUE)]
if (length(missing_objects) > 0) {
  stop(
    "Missing required objects: ",
    paste(missing_objects, collapse = ", "),
    "\nRun Sections 5.1-5.3 first, then run this repair script."
  )
}

required_columns <- c("Gene", "CellType", "Comparison", "avg_log2FC", "p_val_adj")
missing_columns <- setdiff(required_columns, colnames(master_deg))
if (length(missing_columns) > 0) {
  stop("master_deg is missing: ", paste(missing_columns, collapse = ", "))
}

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("LEGEND_DIR")) {
  LEGEND_DIR <- file.path(REPORT_DIR, "Figure_Legends")
}
dir.create(LEGEND_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("save_csv")) {
  save_csv <- function(x, filename) {
    utils::write.csv(x, file.path(TABLE_DIR, filename), row.names = FALSE)
  }
}

if (!exists("log_message")) {
  log_message <- function(x) message("[INFO] ", x)
}

if (!exists("section_timer")) {
  section_timer <- function(label, expr) {
    message("\n>>> ", label)
    t0 <- Sys.time()
    force(expr)
    message("<<< Completed in ",
            round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
            " s")
  }
}

pkg_needed <- c("dplyr", "tidyr", "ggplot2", "ComplexHeatmap", "circlize")
pkg_missing <- pkg_needed[
  !vapply(pkg_needed, requireNamespace, logical(1), quietly = TRUE)
]
if (length(pkg_missing) > 0) {
  stop("Install required packages first: ", paste(pkg_missing, collapse = ", "))
}

###############################################################################
# 1. COMMON SIGNIFICANT-DEG DATASET
###############################################################################

significant_deg <- master_deg |>
  dplyr::filter(
    !is.na(Gene), Gene != "",
    !is.na(CellType), CellType != "",
    !is.na(Comparison), Comparison != "",
    !is.na(avg_log2FC),
    !is.na(p_val_adj),
    p_val_adj < ALPHA,
    abs(avg_log2FC) >= LOGFC_THRESHOLD
  )

if (nrow(significant_deg) == 0) {
  stop(
    "No significant DEGs pass p_val_adj < ", ALPHA,
    " and |avg_log2FC| >= ", LOGFC_THRESHOLD, "."
  )
}

cell_types <- sort(unique(significant_deg$CellType))
comparisons <- sort(unique(significant_deg$Comparison))
n_cell_types <- length(cell_types)
n_comparisons <- length(comparisons)

log_message(paste(
  "Repair analysis uses", nrow(significant_deg), "significant DEG entries,",
  dplyr::n_distinct(significant_deg$Gene), "genes,",
  n_cell_types, "cell types and", n_comparisons, "comparisons."
))

###############################################################################
# SECTION 5.4 : SHARED / UNIQUE DEGs + RECURRENCE VISUALISATION
###############################################################################

section_timer("SECTION 5.4 REPAIR : DEG Recurrence", {

  gene_celltype <- significant_deg |>
    dplyr::distinct(Gene, CellType) |>
    dplyr::mutate(Present = 1L)

  occurrence_matrix <- gene_celltype |>
    tidyr::pivot_wider(
      names_from = CellType,
      values_from = Present,
      values_fill = 0
    )

  occurrence_numeric <- occurrence_matrix |>
    dplyr::select(-Gene) |>
    as.matrix()

  storage.mode(occurrence_numeric) <- "numeric"
  rownames(occurrence_numeric) <- occurrence_matrix$Gene

  comparison_matrix <- significant_deg |>
    dplyr::distinct(Gene, Comparison) |>
    dplyr::mutate(Present = 1L) |>
    tidyr::pivot_wider(
      names_from = Comparison,
      values_from = Present,
      values_fill = 0
    )

  comparison_numeric <- comparison_matrix |>
    dplyr::select(-Gene) |>
    as.matrix()

  storage.mode(comparison_numeric) <- "numeric"
  rownames(comparison_numeric) <- comparison_matrix$Gene

  effect_stats <- significant_deg |>
    dplyr::group_by(Gene) |>
    dplyr::summarise(
      Mean_Log2FC = mean(avg_log2FC, na.rm = TRUE),
      Median_Log2FC = median(avg_log2FC, na.rm = TRUE),
      Mean_Abs_Log2FC = mean(abs(avg_log2FC), na.rm = TRUE),
      Max_Abs_Log2FC = max(abs(avg_log2FC), na.rm = TRUE),
      Min_Adjusted_P = min(p_val_adj, na.rm = TRUE),
      Median_Adjusted_P = median(p_val_adj, na.rm = TRUE),
      .groups = "drop"
    )

  gene_summary <- data.frame(
    Gene = rownames(occurrence_numeric),
    CellType_Frequency = rowSums(occurrence_numeric),
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(
      data.frame(
        Gene = rownames(comparison_numeric),
        Comparison_Frequency = rowSums(comparison_numeric),
        stringsAsFactors = FALSE
      ),
      by = "Gene"
    ) |>
    dplyr::left_join(effect_stats, by = "Gene") |>
    dplyr::mutate(
      Category = dplyr::case_when(
        CellType_Frequency == 1 ~ "Cell-type-specific",
        CellType_Frequency == 2 ~ "Shared across 2 cell types",
        CellType_Frequency %in% 3:4 ~ "Moderately recurrent (3-4)",
        CellType_Frequency >= 5 & CellType_Frequency < n_cell_types ~
          "Highly recurrent (5+)",
        CellType_Frequency == n_cell_types ~ "Conserved across all cell types",
        TRUE ~ "Other"
      )
    )

  save_csv(occurrence_matrix, "5.4_Gene_CellType_Occurrence_Matrix.csv")
  save_csv(comparison_matrix, "5.4_Gene_Comparison_Occurrence_Matrix.csv")
  save_csv(gene_summary, "5.4_Gene_Occurrence_Summary.csv")

  shared_genes <- gene_summary |>
    dplyr::filter(CellType_Frequency >= 2) |>
    dplyr::arrange(
      dplyr::desc(CellType_Frequency),
      dplyr::desc(Comparison_Frequency),
      dplyr::desc(Mean_Abs_Log2FC)
    )

  unique_genes <- gene_summary |>
    dplyr::filter(CellType_Frequency == 1) |>
    dplyr::arrange(dplyr::desc(Mean_Abs_Log2FC))

  conserved_genes <- gene_summary |>
    dplyr::filter(CellType_Frequency == n_cell_types) |>
    dplyr::arrange(dplyr::desc(Mean_Abs_Log2FC))

  save_csv(shared_genes, "5.4_Shared_Genes.csv")
  save_csv(unique_genes, "5.4_Unique_Genes.csv")
  save_csv(conserved_genes, "5.4_Conserved_Genes.csv")

  recurrence_distribution <- gene_summary |>
    dplyr::count(CellType_Frequency, name = "Number_of_Genes") |>
    tidyr::complete(
      CellType_Frequency = seq_len(n_cell_types),
      fill = list(Number_of_Genes = 0)
    ) |>
    dplyr::mutate(
      Recurrence_Class = dplyr::case_when(
        CellType_Frequency == 1 ~ "Cell-type-specific",
        CellType_Frequency == 2 ~ "Shared in 2",
        CellType_Frequency %in% 3:4 ~ "Moderately recurrent (3-4)",
        CellType_Frequency >= 5 & CellType_Frequency < n_cell_types ~
          "Highly recurrent (5+)",
        CellType_Frequency == n_cell_types ~ "All cell types",
        TRUE ~ "Other"
      )
    )

  save_csv(recurrence_distribution, "5.4_Gene_Recurrence_Distribution.csv")

  recurrence_colors <- c(
    "Cell-type-specific" = "#4E79A7",
    "Shared in 2" = "#59A14F",
    "Moderately recurrent (3-4)" = "#F28E2B",
    "Highly recurrent (5+)" = "#E15759",
    "All cell types" = "#7B2CBF",
    "Other" = "#9D9D9D"
  )

  p_recurrence <- ggplot2::ggplot(
    recurrence_distribution,
    ggplot2::aes(
      x = factor(CellType_Frequency, levels = seq_len(n_cell_types)),
      y = Number_of_Genes,
      fill = Recurrence_Class
    )
  ) +
    ggplot2::geom_col(width = 0.76) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(Number_of_Genes == 0, "", Number_of_Genes)),
      vjust = -0.35,
      size = 3.8
    ) +
    ggplot2::scale_fill_manual(values = recurrence_colors, drop = FALSE) +
    ggplot2::labs(
      title = "Recurrence of Significant DEGs Across Immune Cell Types",
      subtitle = paste0(
        "A recurrence value of k indicates that a gene is significantly ",
        "differentially expressed in k of ", n_cell_types, " analysed cell types"
      ),
      x = paste0("Number of cell types containing the significant DEG (1-", n_cell_types, ")"),
      y = "Number of genes",
      fill = "Recurrence class",
      caption = paste0(
        "Significance: adjusted P < ", ALPHA,
        "; |average log2 fold-change| >= ", LOGFC_THRESHOLD, "."
      )
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(size = 10.5),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(
    file.path(FIGURE_DIR, "5.4_Gene_Recurrence_Distribution.png"),
    p_recurrence, width = 11, height = 7.5, dpi = 320, bg = "white"
  )
  ggplot2::ggsave(
    file.path(FIGURE_DIR, "5.4_Gene_Recurrence_Distribution.pdf"),
    p_recurrence, width = 11, height = 7.5, device = cairo_pdf
  )

  # Select top recurrent genes for a readable binary heatmap.
  top_recurrence <- gene_summary |>
    dplyr::filter(CellType_Frequency >= 2) |>
    dplyr::arrange(
      dplyr::desc(CellType_Frequency),
      dplyr::desc(Comparison_Frequency),
      dplyr::desc(Mean_Abs_Log2FC),
      Min_Adjusted_P
    ) |>
    dplyr::slice_head(n = 50)

  save_csv(top_recurrence, "5.4_Top_Recurrence_Genes.csv")

  recurrence_heat_matrix <- occurrence_numeric[
    top_recurrence$Gene, , drop = FALSE
  ]

  # Constant rows remain scientifically reported in tables but are removed
  # from the clustering display because they contain no pattern information.
  variable_rows <- apply(
    recurrence_heat_matrix, 1,
    function(x) length(unique(x)) > 1
  )
  recurrence_heat_matrix <- recurrence_heat_matrix[
    variable_rows, , drop = FALSE
  ]

  if (nrow(recurrence_heat_matrix) > 0) {
    recurrence_ht <- ComplexHeatmap::Heatmap(
      recurrence_heat_matrix,
      name = "Significant DEG",
      col = circlize::colorRamp2(c(0, 1), c("#F7F7F7", "#2166AC")),
      cluster_rows = nrow(recurrence_heat_matrix) > 1,
      cluster_columns = ncol(recurrence_heat_matrix) > 1,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 7),
      column_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
      column_names_rot = 45,
      row_title = "Top recurrent significant DEGs",
      column_title = "Immune cell types",
      heatmap_legend_param = list(
        title = "Significant\nDEG",
        at = c(0, 1),
        labels = c("Absent", "Present")
      ),
      border = TRUE,
      use_raster = FALSE
    )

    grDevices::pdf(
      file.path(FIGURE_DIR, "5.4_Gene_Recurrence_Heatmap.pdf"),
      width = 9.5, height = 12
    )
    ComplexHeatmap::draw(recurrence_ht, heatmap_legend_side = "right")
    grDevices::dev.off()

    grDevices::png(
      file.path(FIGURE_DIR, "5.4_Gene_Recurrence_Heatmap.png"),
      width = 3000, height = 3800, res = 320
    )
    ComplexHeatmap::draw(recurrence_ht, heatmap_legend_side = "right")
    grDevices::dev.off()
  } else {
    warning("No variable recurrent-gene rows available for recurrence heatmap.")
  }

  category_summary <- gene_summary |>
    dplyr::count(Category, name = "Number_of_Genes") |>
    dplyr::arrange(dplyr::desc(Number_of_Genes))
  save_csv(category_summary, "5.4_Gene_Category_Summary.csv")

  report_54 <- data.frame(
    Metric = c(
      "Significant DEG entries",
      "Unique significant genes",
      "Cell types",
      "Comparisons",
      "Cell-type-specific genes",
      "Genes recurrent in >=2 cell types",
      "Genes recurrent in >=5 cell types",
      "Genes conserved across all cell types"
    ),
    Value = c(
      nrow(significant_deg),
      nrow(gene_summary),
      n_cell_types,
      n_comparisons,
      sum(gene_summary$CellType_Frequency == 1),
      sum(gene_summary$CellType_Frequency >= 2),
      sum(gene_summary$CellType_Frequency >= 5),
      sum(gene_summary$CellType_Frequency == n_cell_types)
    )
  )
  save_csv(report_54, "5.4_Report.csv")
})

###############################################################################
# SECTION 5.5 : CONSERVED DISEASE SIGNATURES -- REPAIRED
###############################################################################

section_timer("SECTION 5.5 REPAIR : Conserved Disease Signatures", {

  # IMPORTANT: use the significant-DEG recurrence summary from repaired 5.4.
  # No duplicate frequency columns are joined, preventing .x/.y failures.
  signature_table <- gene_summary |>
    dplyr::mutate(
      Safe_Min_P = pmax(Min_Adjusted_P, .Machine$double.xmin),
      Significance_Strength = -log10(Safe_Min_P),
      Signature_Score =
        Mean_Abs_Log2FC *
        CellType_Frequency *
        Comparison_Frequency *
        Significance_Strength
    ) |>
    dplyr::arrange(dplyr::desc(Signature_Score))

  # A conserved candidate must recur in >=3 cell types and be significant
  # because gene_summary itself was constructed only from significant DEGs.
  conserved_signature <- signature_table |>
    dplyr::filter(CellType_Frequency >= 3) |>
    dplyr::arrange(dplyr::desc(Signature_Score))

  top_signature <- conserved_signature |>
    dplyr::slice_head(n = 40)

  save_csv(signature_table, "5.5_All_Gene_Signatures.csv")
  save_csv(conserved_signature, "5.5_Conserved_Disease_Signatures.csv")
  save_csv(top_signature, "5.5_Top_Conserved_Signatures.csv")

  if (nrow(top_signature) > 0) {
    heat_df <- significant_deg |>
      dplyr::filter(Gene %in% top_signature$Gene) |>
      dplyr::group_by(Gene, CellType) |>
      dplyr::summarise(
        Mean_Log2FC = mean(avg_log2FC, na.rm = TRUE),
        .groups = "drop"
      )

    heat_wide <- heat_df |>
      tidyr::pivot_wider(
        names_from = CellType,
        values_from = Mean_Log2FC,
        values_fill = 0
      )

    signature_mat <- heat_wide |>
      dplyr::select(-Gene) |>
      as.matrix()
    rownames(signature_mat) <- heat_wide$Gene
    storage.mode(signature_mat) <- "numeric"

    # Robust symmetric colour limit: avoids a single extreme value washing
    # out the remainder of the heatmap.
    lim <- stats::quantile(
      abs(signature_mat[is.finite(signature_mat)]),
      probs = 0.95, na.rm = TRUE, names = FALSE
    )
    if (!is.finite(lim) || lim <= 0) lim <- 1

    sig_col_fun <- circlize::colorRamp2(
      c(-lim, 0, lim),
      c("#2166AC", "#F7F7F7", "#B2182B")
    )

    signature_ht <- ComplexHeatmap::Heatmap(
      signature_mat,
      name = "Mean log2FC",
      col = sig_col_fun,
      cluster_rows = nrow(signature_mat) > 1,
      cluster_columns = ncol(signature_mat) > 1,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 7),
      column_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
      column_names_rot = 45,
      row_title = "Top recurrent disease-associated genes",
      column_title = "Immune cell types",
      border = TRUE,
      use_raster = FALSE
    )

    grDevices::pdf(
      file.path(FIGURE_DIR, "5.5_Conserved_Signature_Heatmap.pdf"),
      width = 10, height = 11
    )
    ComplexHeatmap::draw(signature_ht, heatmap_legend_side = "right")
    grDevices::dev.off()

    grDevices::png(
      file.path(FIGURE_DIR, "5.5_Conserved_Signature_Heatmap.png"),
      width = 3200, height = 3500, res = 320
    )
    ComplexHeatmap::draw(signature_ht, heatmap_legend_side = "right")
    grDevices::dev.off()
  }

  report_55 <- data.frame(
    Metric = c(
      "Significant genes evaluated",
      "Genes recurrent in >=3 cell types",
      "Genes shown in conserved-signature heatmap"
    ),
    Value = c(
      nrow(signature_table),
      nrow(conserved_signature),
      nrow(top_signature)
    )
  )
  save_csv(report_55, "5.5_Conserved_Signature_Report.csv")
})

###############################################################################
# SECTION 5.6 : ANALYSIS / FIGURE METADATA
###############################################################################

section_timer("SECTION 5.6 REPAIR : Metadata", {

  dataset_accession <- "GSE227835"

  metadata <- data.frame(
    Field = c(
      "Dataset accession",
      "Analysis level",
      "Number of analysed cell types",
      "Analysed cell types",
      "Number of disease comparisons",
      "Disease comparisons",
      "Adjusted P-value threshold",
      "Absolute log2FC threshold",
      "Significant DEG entries",
      "Unique significant genes",
      "Recurrence definition",
      "Recurrence heatmap selection",
      "Conserved-signature definition",
      "Conserved-signature ranking",
      "Heatmap colour interpretation"
    ),
    Value = c(
      dataset_accession,
      "Cell-type-resolved differential-expression integration",
      n_cell_types,
      paste(cell_types, collapse = "; "),
      n_comparisons,
      paste(comparisons, collapse = "; "),
      as.character(ALPHA),
      as.character(LOGFC_THRESHOLD),
      nrow(significant_deg),
      dplyr::n_distinct(significant_deg$Gene),
      paste0(
        "Number of analysed cell types in which a gene passes adjusted P < ",
        ALPHA, " and |avg_log2FC| >= ", LOGFC_THRESHOLD
      ),
      "Top recurrent genes ranked by cell-type frequency, comparison frequency, effect size and adjusted P value; maximum 50 displayed",
      "Significant DEG recurrent in at least 3 analysed cell types",
      "Mean absolute log2FC x cell-type frequency x comparison frequency x -log10(minimum adjusted P)",
      "Binary recurrence: white = absent, blue = present; effect-size heatmap: blue = lower/negative log2FC, white = zero, red = higher/positive log2FC"
    ),
    stringsAsFactors = FALSE
  )

  save_csv(metadata, "5.6_Analysis_Metadata.csv")

  figure_metadata <- data.frame(
    Figure_File = c(
      "5.4_Gene_Recurrence_Distribution",
      "5.4_Gene_Recurrence_Heatmap",
      "5.5_Conserved_Signature_Heatmap"
    ),
    Purpose = c(
      "Quantifies how broadly significant DEGs recur across analysed immune cell types",
      "Shows cell-type presence/absence patterns for the strongest recurrent significant DEGs",
      "Shows direction and magnitude of mean differential expression for recurrent disease-associated genes"
    ),
    Data_Unit = c(
      "Genes",
      "Gene x cell type",
      "Gene x cell type"
    ),
    stringsAsFactors = FALSE
  )
  save_csv(figure_metadata, "5.6_Figure_Metadata.csv")
})

###############################################################################
# SECTION 5.7 : DETAILED, DATA-DRIVEN FIGURE LEGENDS
###############################################################################

section_timer("SECTION 5.7 REPAIR : Detailed Figure Legends", {

  total_sig_genes <- nrow(gene_summary)
  specific_n <- sum(gene_summary$CellType_Frequency == 1)
  shared2plus_n <- sum(gene_summary$CellType_Frequency >= 2)
  recurrent5_n <- sum(gene_summary$CellType_Frequency >= 5)
  all_n <- sum(gene_summary$CellType_Frequency == n_cell_types)

  specific_pct <- round(100 * specific_n / total_sig_genes, 1)
  shared_pct <- round(100 * shared2plus_n / total_sig_genes, 1)
  recurrent5_pct <- round(100 * recurrent5_n / total_sig_genes, 1)

  dist_text <- paste(
    paste0(
      recurrence_distribution$CellType_Frequency,
      " cell type(s): ",
      recurrence_distribution$Number_of_Genes,
      " gene(s)"
    ),
    collapse = "; "
  )

  legend_recurrence_distribution <- paste0(
    "Figure 1. Recurrence distribution of significant differentially expressed ",
    "genes (DEGs) across immune cell types. The x-axis represents the number ",
    "of analysed immune cell types in which an individual gene met the DEG ",
    "criteria (adjusted P < ", ALPHA, " and |average log2 fold-change| >= ",
    LOGFC_THRESHOLD, "), rather than the identity of a specific cell type. ",
    "The y-axis gives the number of unique genes at each recurrence level. ",
    "Bars are coloured by recurrence class: cell-type-specific (1 cell type), ",
    "shared in 2 cell types, moderately recurrent (3-4), highly recurrent ",
    "(5 or more), and conserved across all analysed cell types when present. ",
    "The analysed cell types were: ", paste(cell_types, collapse = ", "), ". ",
    "Observed recurrence counts were: ", dist_text, ". ",
    specific_n, " of ", total_sig_genes, " significant genes (", specific_pct,
    "%) were cell-type-specific, whereas ", shared2plus_n, " (", shared_pct,
    "%) recurred in at least two cell types. ", recurrent5_n, " genes (",
    recurrent5_pct, "%) recurred in five or more cell types. This distribution ",
    "therefore separates predominantly cell-type-restricted transcriptional ",
    "responses from a smaller recurrent component shared across immune lineages."
  )

  legend_recurrence_heatmap <- paste0(
    "Figure 2. Cell-type recurrence map of significant DEGs. Rows represent ",
    "the highest-ranking recurrent genes and columns represent the analysed ",
    "immune cell types (", paste(cell_types, collapse = ", "), "). A blue cell ",
    "indicates that the gene satisfied the DEG criteria in that cell type; ",
    "white indicates that it did not satisfy those criteria. Genes were ranked ",
    "by cell-type recurrence, comparison recurrence, mean absolute log2 fold-change, ",
    "and adjusted P value, and up to 50 genes were selected for visualization. ",
    "Rows with identical presence across every displayed cell type were excluded ",
    "from clustering because they provide no pattern variation; such genes remain ",
    "retained in the exported analytical tables. Hierarchical clustering groups ",
    "genes and cell types with similar recurrence patterns."
  )

  legend_signature_heatmap <- paste0(
    "Figure 3. Conserved disease-associated transcriptional signatures across ",
    "immune cell types. Rows represent recurrent significant genes detected in ",
    "at least three analysed cell types and prioritized by a composite score ",
    "integrating mean absolute log2 fold-change, cell-type recurrence, comparison ",
    "recurrence and -log10 of the minimum adjusted P value. Columns represent ",
    "the analysed immune cell types. Values are mean log2 fold-changes calculated ",
    "from significant DEG entries for each gene-cell-type combination. Red denotes ",
    "positive mean log2 fold-change, blue denotes negative mean log2 fold-change, ",
    "and white denotes values near zero. Both rows and columns are hierarchically ",
    "clustered to reveal shared versus cell-type-divergent response patterns."
  )

  legends <- data.frame(
    Figure = c("Figure 1", "Figure 2", "Figure 3"),
    Source_File = c(
      "5.4_Gene_Recurrence_Distribution",
      "5.4_Gene_Recurrence_Heatmap",
      "5.5_Conserved_Signature_Heatmap"
    ),
    Legend = c(
      legend_recurrence_distribution,
      legend_recurrence_heatmap,
      legend_signature_heatmap
    ),
    stringsAsFactors = FALSE
  )

  utils::write.csv(
    legends,
    file.path(LEGEND_DIR, "5.7_Detailed_Figure_Legends.csv"),
    row.names = FALSE
  )

  txt <- c(
    "DETAILED FIGURE LEGENDS",
    "=======================",
    "",
    unlist(lapply(seq_len(nrow(legends)), function(i) {
      c(legends$Figure[i], legends$Legend[i], "")
    }))
  )
  writeLines(txt, file.path(LEGEND_DIR, "5.7_Detailed_Figure_Legends.txt"))

  if (requireNamespace("officer", quietly = TRUE)) {
    doc <- officer::read_docx()
    doc <- officer::body_add_par(doc, "Detailed Figure Legends",
                                 style = "heading 1")
    for (i in seq_len(nrow(legends))) {
      doc <- officer::body_add_par(doc, legends$Figure[i], style = "heading 2")
      doc <- officer::body_add_par(doc, legends$Legend[i], style = "Normal")
    }
    print(doc, target = file.path(LEGEND_DIR,
                                  "5.7_Detailed_Figure_Legends.docx"))
  }
})

###############################################################################
# SECTION 6 : QUANTITATIVE RESULTS + BIOLOGICAL INTERPRETATION
###############################################################################

section_timer("SECTION 6 REPAIR : Results Writing", {

  total_sig_genes <- nrow(gene_summary)
  specific_n <- sum(gene_summary$CellType_Frequency == 1)
  shared2_n <- sum(gene_summary$CellType_Frequency == 2)
  moderate_n <- sum(gene_summary$CellType_Frequency %in% 3:4)
  high_n <- sum(gene_summary$CellType_Frequency >= 5 &
                  gene_summary$CellType_Frequency < n_cell_types)
  all_n <- sum(gene_summary$CellType_Frequency == n_cell_types)

  pct <- function(x) round(100 * x / total_sig_genes, 1)

  top_recurrent_text <- if (nrow(shared_genes) > 0) {
    top_n <- min(10L, nrow(shared_genes))
    paste(shared_genes$Gene[seq_len(top_n)], collapse = ", ")
  } else {
    "none"
  }

  top_signature_text <- if (nrow(top_signature) > 0) {
    top_n <- min(10L, nrow(top_signature))
    paste(top_signature$Gene[seq_len(top_n)], collapse = ", ")
  } else {
    "none"
  }

  results_paragraphs <- c(
    "6. BIOLOGICAL INTERPRETATION / RESULTS",
    "",
    "6.1 Integrated significant-DEG landscape",
    paste0(
      "Applying the predefined thresholds (adjusted P < ", ALPHA,
      " and |average log2 fold-change| >= ", LOGFC_THRESHOLD, ") identified ",
      nrow(significant_deg), " significant DEG entries representing ",
      total_sig_genes, " unique genes across ", n_cell_types,
      " immune cell types and ", n_comparisons, " disease comparisons. ",
      "The analysed cell types were ", paste(cell_types, collapse = ", "),
      ". The comparisons represented in the integrated dataset were ",
      paste(comparisons, collapse = ", "), "."
    ),
    "",
    "6.2 Recurrence across immune cell types",
    paste0(
      "The recurrence analysis demonstrated that the transcriptional response ",
      "was dominated by cell-type-restricted changes (Figure 1). Of ",
      total_sig_genes, " significant genes, ", specific_n, " (", pct(specific_n),
      "%) were significant in only one immune cell type. A further ", shared2_n,
      " (", pct(shared2_n), "%) were shared by exactly two cell types, while ",
      moderate_n, " (", pct(moderate_n), "%) recurred across three to four ",
      "cell types. Only ", high_n + all_n, " genes (",
      pct(high_n + all_n), "%) were detected in five or more cell types. ",
      "This pattern indicates substantial cell-type specificity alongside a ",
      "smaller core of recurrent disease-associated transcriptional changes."
    ),
    "",
    "6.3 Recurrent gene patterns",
    paste0(
      "The binary recurrence heatmap (Figure 2) resolves which immune cell ",
      "types contribute to the recurrent component. Unlike the recurrence ",
      "distribution, where the x-axis denotes recurrence frequency, the heatmap ",
      "columns correspond to individual cell-type identities. The highest-ranking ",
      "recurrent genes included ", top_recurrent_text, ". These candidates were ",
      "prioritized using recurrence across cell types and comparisons together ",
      "with effect-size information, allowing broadly shared changes to be ",
      "distinguished from lineage-restricted signals."
    ),
    "",
    "6.4 Conserved disease-associated signatures",
    paste0(
      nrow(conserved_signature), " genes were recurrent in at least three cell ",
      "types and were therefore carried forward as conserved/recurrent candidate ",
      "signatures. Figure 3 displays the highest-ranking candidates using mean ",
      "log2 fold-change rather than binary presence, thereby showing both the ",
      "direction and magnitude of the transcriptional response. The leading ",
      "ranked candidates included ", top_signature_text, ". A gene may recur ",
      "across several cell types while differing in direction or magnitude; ",
      "therefore recurrence and effect-size heatmaps should be interpreted ",
      "together rather than as interchangeable measures."
    ),
    "",
    "6.5 Overall interpretation",
    paste0(
      "Taken together, the integrated analysis supports a model in which most ",
      "significant transcriptional alterations are cell-type-specific, while a ",
      "smaller set of genes recur across multiple immune populations. The latter ",
      "group represents candidate shared disease-associated mechanisms, whereas ",
      "the cell-type-specific component may reflect lineage-dependent responses. ",
      "These findings are computational associations from the present dataset and ",
      "should be treated as candidates for pathway-level interpretation and ",
      "independent experimental or external-dataset validation rather than as ",
      "validated biomarkers."
    )
  )

  writeLines(
    results_paragraphs,
    file.path(REPORT_DIR, "Section6_Quantitative_Results_and_Interpretation.txt")
  )

  if (requireNamespace("officer", quietly = TRUE)) {
    doc <- officer::read_docx()
    for (line in results_paragraphs) {
      if (grepl("^6\\.", line)) {
        doc <- officer::body_add_par(doc, line, style = "heading 2")
      } else if (line != "") {
        doc <- officer::body_add_par(doc, line, style = "Normal")
      }
    }
    print(
      doc,
      target = file.path(
        REPORT_DIR,
        "Section6_Quantitative_Results_and_Interpretation.docx"
      )
    )
  }

  result_metrics <- data.frame(
    Metric = c(
      "Unique significant genes",
      "Cell-type-specific genes",
      "Shared in exactly 2 cell types",
      "Recurrent in 3-4 cell types",
      "Recurrent in >=5 cell types",
      "Conserved/recurrent candidates (>=3 cell types)"
    ),
    Count = c(
      total_sig_genes, specific_n, shared2_n, moderate_n,
      high_n + all_n, nrow(conserved_signature)
    ),
    Percent_of_significant_genes = c(
      100, pct(specific_n), pct(shared2_n), pct(moderate_n),
      pct(high_n + all_n), pct(nrow(conserved_signature))
    )
  )
  save_csv(result_metrics, "Section6_Key_Result_Metrics.csv")
})

###############################################################################
# SECTION 7 : MANUSCRIPT-READY SUPPLEMENTARY TABLES
###############################################################################

section_timer("SECTION 7 REPAIR : Supplementary Tables", {

  supp_dir <- file.path(REPORT_DIR, "Supplementary")
  dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(
    significant_deg,
    file.path(supp_dir, "Table_S1_All_Significant_DEG_Entries.csv"),
    row.names = FALSE
  )
  write.csv(
    gene_summary,
    file.path(supp_dir, "Table_S2_Gene_Recurrence_and_Effect_Summary.csv"),
    row.names = FALSE
  )
  write.csv(
    shared_genes,
    file.path(supp_dir, "Table_S3_Shared_Recurrent_Genes.csv"),
    row.names = FALSE
  )
  write.csv(
    unique_genes,
    file.path(supp_dir, "Table_S4_CellType_Specific_Genes.csv"),
    row.names = FALSE
  )
  write.csv(
    conserved_signature,
    file.path(supp_dir, "Table_S5_Conserved_Disease_Signatures.csv"),
    row.names = FALSE
  )
  write.csv(
    recurrence_distribution,
    file.path(supp_dir, "Table_S6_Recurrence_Distribution.csv"),
    row.names = FALSE
  )

  supp_index <- data.frame(
    Table = paste0("S", 1:6),
    Description = c(
      "All significant DEG entries passing the integrated thresholds",
      "Gene-level recurrence, comparison frequency and effect-size statistics",
      "Genes significant in at least two analysed cell types",
      "Genes significant in exactly one analysed cell type",
      "Recurrent disease-associated candidate signatures",
      "Number of genes observed at each cell-type recurrence frequency"
    )
  )
  write.csv(
    supp_index,
    file.path(supp_dir, "Supplementary_Table_Index.csv"),
    row.names = FALSE
  )
})

###############################################################################
# SECTION 8 : METHODOLOGY PARAGRAPH + COLOURED FLOWCHART
###############################################################################

section_timer("SECTION 8 REPAIR : Methodology", {

  methodology_paragraph <- paste0(
    "Integrated differential-expression analysis was performed using the ",
    "cell-type-resolved scRNA-seq outputs from GEO accession GSE227835. ",
    "Differential-expression results generated for individual immune cell types ",
    "and disease comparisons were consolidated into a master DEG table. Entries ",
    "were considered significant when the adjusted P value was < ", ALPHA,
    " and the absolute average log2 fold-change was >= ", LOGFC_THRESHOLD,
    ". Significant genes were then mapped across ", n_cell_types,
    " analysed immune cell types (", paste(cell_types, collapse = ", "),
    ") and ", n_comparisons, " comparisons (", paste(comparisons, collapse = ", "),
    "). For each gene, cell-type recurrence was defined as the number of immune ",
    "cell types in which the gene met the significance criteria, while comparison ",
    "recurrence represented the number of disease contrasts in which it was ",
    "detected. Genes were classified as cell-type-specific or recurrent according ",
    "to these frequencies. Recurrent genes were visualized using a binary ",
    "presence/absence heatmap, while conserved disease-associated candidates were ",
    "prioritized using a composite score integrating mean absolute log2 fold-change, ",
    "cell-type recurrence, comparison recurrence and adjusted-P-value strength. ",
    "For the conserved-signature heatmap, mean log2 fold-change values were ",
    "summarized for each gene-cell-type combination to retain both response ",
    "direction and magnitude. Analytical metadata, detailed figure legends, ",
    "quantitative Results text and manuscript-ready supplementary tables were ",
    "generated from the same filtered dataset to maintain consistency between ",
    "statistical outputs, visualizations and reporting."
  )

  writeLines(
    methodology_paragraph,
    file.path(REPORT_DIR, "Section8_Methodology_Paragraph.txt")
  )

  if (requireNamespace("officer", quietly = TRUE)) {
    doc <- officer::read_docx()
    doc <- officer::body_add_par(doc, "Methodology", style = "heading 1")
    doc <- officer::body_add_par(doc, methodology_paragraph, style = "Normal")
    print(
      doc,
      target = file.path(REPORT_DIR, "Section8_Methodology_Paragraph.docx")
    )
  }

  # DiagrammeR flowchart is optional: text/report outputs still complete if
  # DiagrammeR-related packages are unavailable.
  diagram_pkgs <- c("DiagrammeR", "DiagrammeRsvg", "rsvg")
  diagram_ok <- all(vapply(
    diagram_pkgs, requireNamespace, logical(1), quietly = TRUE
  ))

  if (diagram_ok) {

    cell_label <- paste(cell_types, collapse = ", ")
    comparison_label <- paste(comparisons, collapse = ", ")

    # Keep node labels concise; complete metadata are in 5.6_Analysis_Metadata.csv.
    graph_code <- paste0(
'digraph methodology {
graph [
  layout = dot,
  rankdir = LR,
  bgcolor = "white",
  nodesep = 0.45,
  ranksep = 0.65,
  pad = 0.35
]

node [
  shape = rectangle,
  style = "rounded,filled",
  fontname = "Helvetica",
  fontsize = 11,
  color = "#334155",
  penwidth = 1.2,
  margin = "0.18,0.12"
]

edge [
  color = "#64748B",
  penwidth = 1.5,
  arrowsize = 0.75
]

A [label = "GSE227835\\nscRNA-seq dataset", fillcolor = "#DBEAFE"]
B [label = "Cell-type-resolved\\nDEG results", fillcolor = "#DBEAFE"]

C [label = "Integrated master\\nDEG table", fillcolor = "#E0E7FF"]
D [label = "Significance filtering\\nadj. P < ', ALPHA,
'\\n|avg log2FC| >= ', LOGFC_THRESHOLD, '", fillcolor = "#E0E7FF"]

E [label = "Gene x cell-type\\nrecurrence matrix", fillcolor = "#DCFCE7"]
F [label = "Gene x comparison\\nrecurrence matrix", fillcolor = "#DCFCE7"]

G [label = "Cell-type-specific\\nand shared DEGs", fillcolor = "#FEF3C7"]
H [label = "Conserved / recurrent\\ndisease signatures", fillcolor = "#FEF3C7"]

I [label = "Recurrence\\nvisualisations", fillcolor = "#FCE7F3"]
J [label = "Effect-size\\nheatmaps", fillcolor = "#FCE7F3"]

K [label = "Metadata + legends\\n+ quantitative results", fillcolor = "#F3E8FF"]
L [label = "Supplementary tables\\n+ manuscript reporting", fillcolor = "#F3E8FF"]

A -> B -> C -> D
D -> E
D -> F
E -> G
F -> G
G -> H
G -> I
H -> J
I -> K
J -> K
K -> L
}'
    )

    graph <- DiagrammeR::grViz(graph_code)
    svg_txt <- DiagrammeRsvg::export_svg(graph)

    svg_file <- file.path(FIGURE_DIR, "Section8_Methodology_Workflow.svg")
    writeLines(svg_txt, svg_file)

    rsvg::rsvg_png(
      svg_file,
      file.path(FIGURE_DIR, "Section8_Methodology_Workflow.png"),
      width = 3600
    )
    rsvg::rsvg_pdf(
      svg_file,
      file.path(FIGURE_DIR, "Section8_Methodology_Workflow.pdf")
    )

  } else {
    warning(
      "Flowchart skipped because one or more packages are missing: ",
      paste(diagram_pkgs[!vapply(
        diagram_pkgs, requireNamespace, logical(1), quietly = TRUE
      )], collapse = ", "),
      ". Methodology paragraph was still generated."
    )
  }
})

###############################################################################
# FINAL OUTPUT MANIFEST
###############################################################################

manifest_files <- c(
  "5.4_Gene_CellType_Occurrence_Matrix.csv",
  "5.4_Gene_Comparison_Occurrence_Matrix.csv",
  "5.4_Gene_Occurrence_Summary.csv",
  "5.4_Gene_Recurrence_Distribution.csv",
  "5.4_Top_Recurrence_Genes.csv",
  "5.5_All_Gene_Signatures.csv",
  "5.5_Conserved_Disease_Signatures.csv",
  "5.6_Analysis_Metadata.csv",
  "5.6_Figure_Metadata.csv",
  "Section6_Key_Result_Metrics.csv"
)

manifest <- data.frame(
  File = manifest_files,
  Exists = file.exists(file.path(TABLE_DIR, manifest_files)),
  stringsAsFactors = FALSE
)
save_csv(manifest, "Reporting_Repair_Output_Manifest.csv")

cat("\n============================================================\n")
cat(" DOWNSTREAM REPORTING / VISUALISATION REPAIR COMPLETED\n")
cat("============================================================\n")
cat("Significant DEG entries :", nrow(significant_deg), "\n")
cat("Unique significant genes:", nrow(gene_summary), "\n")
cat("Cell types              :", n_cell_types, "\n")
cat("Comparisons             :", n_comparisons, "\n")
cat("Figures                 :", FIGURE_DIR, "\n")
cat("Tables                  :", TABLE_DIR, "\n")
cat("Reports                 :", REPORT_DIR, "\n")
cat("============================================================\n")
