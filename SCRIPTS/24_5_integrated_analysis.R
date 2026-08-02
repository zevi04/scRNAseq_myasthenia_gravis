###############################################################################
# SECTION 5.1 : LOAD & VALIDATE MASTER DEG TABLE
###############################################################################

section_timer("SECTION 5.1 : Load & Validate Master DEG Table", {
  
  cat("\n=================================================\n")
  cat(" GLOBAL DIFFERENTIAL EXPRESSION LANDSCAPE\n")
  cat(" SECTION 5.1 : MASTER DEG INITIALIZATION\n")
  cat("=================================================\n\n")
  
  log_message("Loading Master_DEG_Table.csv")
  
  #########################################################################
  # File existence
  #########################################################################
  
  master_file <-
    file.path(
      PATHS$tables,
      "Master_DEG_Table.csv"
    )
  
  if(!file.exists(master_file)){
    
    stop(
      paste(
        "Master DEG table not found:",
        master_file
      )
    )
    
  }
  
  #########################################################################
  # Load table
  #########################################################################
  
  master_deg <-
    read.csv(
      master_file,
      stringsAsFactors = FALSE
    )
  
  log_message(
    paste(
      "Master DEG table loaded:",
      nrow(master_deg),
      "rows"
    )
  )
  
  #########################################################################
  # Required columns
  #########################################################################
  
  required_columns <- c(
    
    "Gene",
    "p_val",
    "avg_log2FC",
    "pct.1",
    "pct.2",
    "p_val_adj",
    "CellType",
    "Comparison",
    "Direction"
    
  )
  
  missing_columns <-
    setdiff(
      required_columns,
      colnames(master_deg)
    )
  
  if(length(missing_columns) > 0){
    
    stop(
      
      paste(
        
        "Missing columns:",
        
        paste(
          missing_columns,
          collapse = ", "
        )
        
      )
      
    )
    
  }
  
  #########################################################################
  # Data type conversion
  #########################################################################
  
  master_deg$Gene <-
    as.character(master_deg$Gene)
  
  master_deg$CellType <-
    factor(master_deg$CellType)
  
  master_deg$Comparison <-
    factor(
      master_deg$Comparison,
      levels = c(
        "Healthy_vs_AChR_MG",
        "Healthy_vs_SNMG_Pre",
        "SNMG_Pre_vs_SNMG_Post"
      )
    )
  
  master_deg$Direction <-
    factor(
      master_deg$Direction,
      levels = c(
        "Upregulated",
        "Downregulated"
      )
    )
  
  numeric_columns <- c(
    
    "p_val",
    "avg_log2FC",
    "pct.1",
    "pct.2",
    "p_val_adj"
    
  )
  
  master_deg[numeric_columns] <-
    lapply(
      master_deg[numeric_columns],
      as.numeric
    )
  
  #########################################################################
  # Missing values
  #########################################################################
  
  missing_summary <-
    
    data.frame(
      
      Variable = names(master_deg),
      
      MissingValues =
        sapply(
          master_deg,
          function(x)
            sum(is.na(x))
        ),
      
      stringsAsFactors = FALSE
      
    )
  
  #########################################################################
  # Duplicate entries
  #########################################################################
  
  duplicate_rows <-
    duplicated(master_deg)
  
  n_duplicates <-
    sum(duplicate_rows)
  
  #########################################################################
  # Dataset statistics
  #########################################################################
  
  dataset_summary <-
    
    data.frame(
      
      Metric = c(
        
        "Total DEG Entries",
        
        "Unique Genes",
        
        "Cell Types",
        
        "Comparisons",
        
        "Upregulated Genes",
        
        "Downregulated Genes",
        
        "Duplicate Rows"
        
      ),
      
      Value = c(
        
        nrow(master_deg),
        
        length(unique(master_deg$Gene)),
        
        length(unique(master_deg$CellType)),
        
        length(unique(master_deg$Comparison)),
        
        sum(master_deg$Direction=="Upregulated"),
        
        sum(master_deg$Direction=="Downregulated"),
        
        n_duplicates
        
      ),
      
      stringsAsFactors = FALSE
      
    )
  
  #########################################################################
  # Save reports
  #########################################################################
  
  save_csv(
    dataset_summary,
    "5.1_Master_DEG_Dataset_Summary.csv"
  )
  
  save_csv(
    missing_summary,
    "5.1_Missing_Value_Report.csv"
  )
  
  #########################################################################
  # Critical validation
  #########################################################################
  
  if(any(missing_summary$MissingValues > 0)){
    
    warning(
      "Missing values detected in Master DEG table."
    )
    
    log_message(
      "WARNING : Missing values detected."
    )
    
  }
  
  if(n_duplicates > 0){
    
    warning(
      paste(
        n_duplicates,
        "duplicate rows detected."
      )
    )
    
    log_message(
      paste(
        "WARNING:",
        n_duplicates,
        "duplicate rows detected."
      )
    )
    
  }
  
  #########################################################################
  # Console summary
  #########################################################################
  
  cat("\n=========================================\n")
  cat(" MASTER DEG DATASET SUMMARY\n")
  cat("=========================================\n\n")
  
  print(dataset_summary)
  
  cat("\nCell Types\n")
  print(sort(unique(master_deg$CellType)))
  
  cat("\nComparisons\n")
  print(levels(master_deg$Comparison))
  
  cat("\nDirection Counts\n")
  
  print(
    table(master_deg$Direction)
  )
  
  cat("\n=========================================\n")
  cat(" Section 5.1 Completed Successfully\n")
  cat("=========================================\n")
  
  log_message("Section 5.1 completed successfully.")
  
})

###############################################################################
# SECTION 5.2 : GLOBAL DEG OVERVIEW
###############################################################################

section_timer("SECTION 5.2 : Global DEG Overview", {
  
  cat("\n=================================================\n")
  cat(" SECTION 5.2 : GLOBAL DEG OVERVIEW\n")
  cat("=================================================\n\n")
  
  log_message("Generating global DEG overview.")
  
  #########################################################################
  # Overall DEG statistics
  #########################################################################
  
  overall_summary <- data.frame(
    
    Metric = c(
      
      "Total DEG Entries",
      "Unique Genes",
      "Unique Cell Types",
      "Unique Comparisons",
      "Upregulated",
      "Downregulated"
      
    ),
    
    Value = c(
      
      nrow(master_deg),
      
      dplyr::n_distinct(master_deg$Gene),
      
      dplyr::n_distinct(master_deg$CellType),
      
      dplyr::n_distinct(master_deg$Comparison),
      
      sum(master_deg$Direction == "Upregulated"),
      
      sum(master_deg$Direction == "Downregulated")
      
    ),
    
    stringsAsFactors = FALSE
    
  )
  
  save_csv(
    overall_summary,
    "5.2_Overall_DEG_Summary.csv"
  )
  
  #########################################################################
  # Cell-type DEG summary
  #########################################################################
  
  celltype_summary <-
    
    master_deg %>%
    
    group_by(CellType) %>%
    
    summarise(
      
      Total_DEGs = n(),
      
      Upregulated = sum(Direction == "Upregulated"),
      
      Downregulated = sum(Direction == "Downregulated"),
      
      Unique_Genes = n_distinct(Gene),
      
      Mean_Log2FC = mean(avg_log2FC),
      
      Median_Log2FC = median(avg_log2FC),
      
      Mean_Adjusted_P =
        mean(p_val_adj),
      
      .groups = "drop"
      
    ) %>%
    
    arrange(desc(Total_DEGs))
  
  save_csv(
    
    celltype_summary,
    
    "5.2_CellType_DEG_Summary.csv"
    
  )
  
  #########################################################################
  # Comparison summary
  #########################################################################
  
  comparison_summary <-
    
    master_deg %>%
    
    group_by(Comparison) %>%
    
    summarise(
      
      Total_DEGs = n(),
      
      Upregulated = sum(Direction=="Upregulated"),
      
      Downregulated = sum(Direction=="Downregulated"),
      
      Unique_Genes = n_distinct(Gene),
      
      Mean_Log2FC = mean(avg_log2FC),
      
      Median_Log2FC = median(avg_log2FC),
      
      Mean_Adjusted_P =
        mean(p_val_adj),
      
      .groups="drop"
      
    ) %>%
    
    arrange(desc(Total_DEGs))
  
  save_csv(
    
    comparison_summary,
    
    "5.2_Comparison_DEG_Summary.csv"
    
  )
  
  #########################################################################
  # Direction summary
  #########################################################################
  
  direction_summary <-
    
    master_deg %>%
    
    group_by(Direction) %>%
    
    summarise(
      
      Total = n(),
      
      Percentage =
        round(
          100*n()/nrow(master_deg),
          2
        ),
      
      .groups="drop"
      
    )
  
  save_csv(
    
    direction_summary,
    
    "5.2_Direction_Summary.csv"
    
  )
  
  #########################################################################
  # Significant DEG summary
  #########################################################################
  
  significant_deg <-
    
    master_deg %>%
    
    filter(
      
      p_val_adj < ALPHA,
      
      abs(avg_log2FC) >= LOGFC_THRESHOLD
      
    )
  
  significant_summary <-
    
    significant_deg %>%
    
    group_by(
      
      CellType,
      
      Comparison
      
    ) %>%
    
    summarise(
      
      Significant_DEGs = n(),
      
      Upregulated =
        sum(Direction=="Upregulated"),
      
      Downregulated =
        sum(Direction=="Downregulated"),
      
      .groups="drop"
      
    )
  
  save_csv(
    
    significant_summary,
    
    "5.2_Significant_DEG_Summary.csv"
    
  )
  
  #########################################################################
  # Console output
  #########################################################################
  
  cat("\n=========================================\n")
  cat(" GLOBAL DEG OVERVIEW COMPLETED\n")
  cat("=========================================\n\n")
  
  print(overall_summary)
  
  cat("\nTop Cell Types\n")
  
  print(celltype_summary)
  
  cat("\nComparison Summary\n")
  
  print(comparison_summary)
  
  log_message("Section 5.2 completed successfully.")
  
})

###############################################################################
# SECTION 5.3 : GLOBAL DEG VISUALIZATION
###############################################################################

section_timer("SECTION 5.3 : Global DEG Visualization", {
  
  cat("\n=================================================\n")
  cat(" SECTION 5.3 : GLOBAL DEG VISUALIZATION\n")
  cat("=================================================\n\n")
  
  log_message("Generating global DEG visualizations.")
  
  #########################################################################
  # Load summary tables
  #########################################################################
  
  overall_summary <-
    read.csv(
      file.path(TABLE_DIR,
                "5.2_Overall_DEG_Summary.csv")
    )
  
  celltype_summary <-
    read.csv(
      file.path(TABLE_DIR,
                "5.2_CellType_DEG_Summary.csv")
    )
  
  comparison_summary <-
    read.csv(
      file.path(TABLE_DIR,
                "5.2_Comparison_DEG_Summary.csv")
    )
  
  direction_summary <-
    read.csv(
      file.path(TABLE_DIR,
                "5.2_Direction_Summary.csv")
    )
  
  significant_summary <-
    read.csv(
      file.path(TABLE_DIR,
                "5.2_Significant_DEG_Summary.csv")
    )
  
  #########################################################################
  # Cell-type DEG count
  #########################################################################
  
  p_celltype <-
    
    ggplot(
      
      celltype_summary,
      
      aes(
        
        x = reorder(CellType, Total_DEGs),
        
        y = Total_DEGs
        
      )
      
    ) +
    
    geom_col(
      
      fill = "#2C7FB8",
      
      width = 0.8
      
    ) +
    
    coord_flip() +
    
    geom_text(
      
      aes(label = Total_DEGs),
      
      hjust = -0.2,
      
      size = 4
      
    ) +
    
    labs(
      
      title = "Differentially Expressed Genes by Cell Type",
      
      x = "",
      
      y = "Number of DEGs"
      
    )
  
  save_plot(
    
    p_celltype,
    
    "5.3_CellType_DEG_Count"
    
  )
  
  #########################################################################
  # Comparison DEG count
  #########################################################################
  
  p_comparison <-
    
    ggplot(
      
      comparison_summary,
      
      aes(
        
        x = reorder(Comparison, Total_DEGs),
        
        y = Total_DEGs
        
      )
      
    ) +
    
    geom_col(
      
      fill = "#D95F02",
      
      width = 0.75
      
    ) +
    
    geom_text(
      
      aes(label = Total_DEGs),
      
      vjust = -0.3,
      
      size = 4
      
    ) +
    
    labs(
      
      title = "Differentially Expressed Genes Across Comparisons",
      
      x = "",
      
      y = "Number of DEGs"
      
    ) +
    
    theme(
      
      axis.text.x =
        
        element_text(
          
          angle = 20,
          
          hjust = 1
          
        )
      
    )
  
  save_plot(
    
    p_comparison,
    
    "5.3_Comparison_DEG_Count"
    
  )
  
  #########################################################################
  # Upregulated vs Downregulated
  #########################################################################
  
  stacked_df <-
    
    celltype_summary %>%
    
    select(
      
      CellType,
      
      Upregulated,
      
      Downregulated
      
    ) %>%
    
    pivot_longer(
      
      cols =
        
        c(
          
          Upregulated,
          
          Downregulated
          
        ),
      
      names_to = "Direction",
      
      values_to = "Count"
      
    )
  
  p_direction <-
    
    ggplot(
      
      stacked_df,
      
      aes(
        
        x = reorder(CellType, Count),
        
        y = Count,
        
        fill = Direction
        
      )
      
    ) +
    
    geom_col() +
    
    coord_flip() +
    
    scale_fill_manual(
      
      values = c(
        
        Upregulated = "#D73027",
        
        Downregulated = "#4575B4"
        
      )
      
    ) +
    
    labs(
      
      title = "Upregulated vs Downregulated Genes",
      
      x = "",
      
      y = "Number of DEGs"
      
    )
  
  save_plot(
    
    p_direction,
    
    "5.3_Up_vs_Down"
    
  )
  
  #########################################################################
  # Direction proportion
  #########################################################################
  
  p_pie <-
    
    ggplot(
      
      direction_summary,
      
      aes(
        
        x = "",
        
        y = Percentage,
        
        fill = Direction
        
      )
      
    ) +
    
    geom_col(
      
      width = 1
      
    ) +
    
    coord_polar("y") +
    
    geom_text(
      
      aes(
        
        label =
          
          paste0(
            
            Percentage,
            
            "%"
            
          )
        
      ),
      
      position =
        
        position_stack(
          
          vjust = 0.5
          
        )
      
    ) +
    
    scale_fill_manual(
      
      values = c(
        
        Upregulated = "#D73027",
        
        Downregulated = "#4575B4"
        
      )
      
    ) +
    
    labs(
      
      title = "Overall DEG Direction"
      
    ) +
    
    theme_void()
  
  save_plot(
    
    p_pie,
    
    "5.3_DEG_Direction_Pie"
    
  )
  
  #########################################################################
  # Significant DEG heatmap
  #########################################################################
  
  heatmap_matrix <-
    
    significant_summary %>%
    
    select(
      
      CellType,
      
      Comparison,
      
      Significant_DEGs
      
    ) %>%
    
    pivot_wider(
      
      names_from = Comparison,
      
      values_from = Significant_DEGs,
      
      values_fill = 0
      
    )
  
  heat_matrix <-
    
    as.matrix(
      
      heatmap_matrix[,-1]
      
    )
  
  rownames(heat_matrix) <-
    
    heatmap_matrix$CellType
  
  col_fun <-
    
    circlize::colorRamp2(
      
      c(
        
        min(heat_matrix),
        
        median(heat_matrix),
        
        max(heat_matrix)
        
      ),
      
      c(
        
        "white",
        
        "#FDB863",
        
        "#B2182B"
        
      )
      
    )
  
  ht <-
    
    ComplexHeatmap::Heatmap(
      
      heat_matrix,
      
      name = "DEGs",
      
      col = col_fun,
      
      cluster_rows = TRUE,
      
      cluster_columns = TRUE,
      
      row_names_side = "left",
      
      column_title = "Comparison",
      
      row_title = "Cell Type"
      
    )
  
  pdf(
    
    file.path(
      
      FIGURE_DIR,
      
      "5.3_Significant_DEG_Heatmap.pdf"
      
    ),
    
    width = 8,
    
    height = 7
    
  )
  
  ComplexHeatmap::draw(ht)
  
  dev.off()
  
  png(
    
    file.path(
      
      FIGURE_DIR,
      
      "5.3_Significant_DEG_Heatmap.png"
      
    ),
    
    width = 2200,
    
    height = 1800,
    
    res = 300
    
  )
  
  ComplexHeatmap::draw(ht)
  
  dev.off()
  
  #########################################################################
  # Summary report
  #########################################################################
  
  figure_summary <-
    
    data.frame(
      
      Figure = c(
        
        "Cell Type DEG Count",
        
        "Comparison DEG Count",
        
        "Up vs Down",
        
        "Direction Pie",
        
        "Significant DEG Heatmap"
        
      ),
      
      Status = "Generated",
      
      stringsAsFactors = FALSE
      
    )
  
  save_csv(
    
    figure_summary,
    
    "5.3_Figure_Summary.csv"
    
  )
  
  #########################################################################
  # Console output
  #########################################################################
  
  cat("\n=========================================\n")
  cat(" GLOBAL VISUALIZATION COMPLETED\n")
  cat("=========================================\n\n")
  
  cat("Figures generated : 5\n")
  
  log_message("Section 5.3 completed successfully.")
  
})

###############################################################################
# SECTION 5.4 : SHARED & UNIQUE DIFFERENTIALLY EXPRESSED GENES
###############################################################################

section_timer("SECTION 5.4 : Shared & Unique DEGs", {
  
  cat("\n=================================================\n")
  cat(" SECTION 5.4 : SHARED & UNIQUE DEGs\n")
  cat("=================================================\n\n")
  
  log_message("Building significant DEG occurrence matrices.")
  
  #########################################################################
  # 1. Validate required columns
  #########################################################################
  
  required_cols <- c(
    "Gene",
    "CellType",
    "Comparison",
    "avg_log2FC",
    "p_val_adj"
  )
  
  missing_cols <- setdiff(
    required_cols,
    colnames(master_deg)
  )
  
  if (length(missing_cols) > 0) {
    
    stop(
      paste(
        "Section 5.4 missing required columns:",
        paste(missing_cols, collapse = ", ")
      )
    )
    
  }
  
  #########################################################################
  # 2. Select significant DEGs only
  #########################################################################
  
  significant_deg <-
    master_deg |>
    dplyr::filter(
      !is.na(Gene),
      Gene != "",
      !is.na(CellType),
      !is.na(Comparison),
      !is.na(avg_log2FC),
      !is.na(p_val_adj),
      p_val_adj < ALPHA,
      abs(avg_log2FC) >= LOGFC_THRESHOLD
    )
  
  if (nrow(significant_deg) == 0) {
    
    stop(
      paste0(
        "Section 5.4: No significant DEGs remain after filtering. ",
        "Check ALPHA, LOGFC_THRESHOLD, and master_deg."
      )
    )
    
  }
  
  cat(
    "Significant DEG entries used :",
    nrow(significant_deg),
    "\n"
  )
  
  cat(
    "Unique significant genes     :",
    dplyr::n_distinct(significant_deg$Gene),
    "\n"
  )
  
  #########################################################################
  # 3. Gene × Cell Type occurrence matrix
  #########################################################################
  
  gene_celltype <-
    significant_deg |>
    dplyr::distinct(
      Gene,
      CellType
    ) |>
    dplyr::mutate(
      Present = 1L
    )
  
  occurrence_matrix <-
    gene_celltype |>
    tidyr::pivot_wider(
      names_from = CellType,
      values_from = Present,
      values_fill = 0
    )
  
  occurrence_numeric <-
    occurrence_matrix |>
    dplyr::select(-Gene) |>
    as.matrix()
  
  storage.mode(occurrence_numeric) <- "numeric"
  
  rownames(occurrence_numeric) <-
    occurrence_matrix$Gene
  
  #########################################################################
  # 4. Cell-type recurrence frequency
  #########################################################################
  
  occurrence_summary <-
    data.frame(
      Gene = rownames(occurrence_numeric),
      CellType_Frequency = rowSums(
        occurrence_numeric,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  
  #########################################################################
  # 5. Gene × Comparison occurrence matrix
  #########################################################################
  
  comparison_matrix <-
    significant_deg |>
    dplyr::distinct(
      Gene,
      Comparison
    ) |>
    dplyr::mutate(
      Present = 1L
    ) |>
    tidyr::pivot_wider(
      names_from = Comparison,
      values_from = Present,
      values_fill = 0
    )
  
  comparison_numeric <-
    comparison_matrix |>
    dplyr::select(-Gene) |>
    as.matrix()
  
  storage.mode(comparison_numeric) <- "numeric"
  
  rownames(comparison_numeric) <-
    comparison_matrix$Gene
  
  comparison_summary <-
    data.frame(
      Gene = rownames(comparison_numeric),
      Comparison_Frequency = rowSums(
        comparison_numeric,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  
  #########################################################################
  # 6. Gene-level effect-size statistics
  #########################################################################
  
  gene_effect_summary <-
    significant_deg |>
    dplyr::group_by(Gene) |>
    dplyr::summarise(
      
      Mean_Log2FC =
        mean(
          avg_log2FC,
          na.rm = TRUE
        ),
      
      Mean_Abs_Log2FC =
        mean(
          abs(avg_log2FC),
          na.rm = TRUE
        ),
      
      Max_Abs_Log2FC =
        max(
          abs(avg_log2FC),
          na.rm = TRUE
        ),
      
      Min_Adjusted_P =
        min(
          p_val_adj,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  #########################################################################
  # 7. Merge recurrence statistics
  #########################################################################
  
  gene_summary <-
    occurrence_summary |>
    dplyr::left_join(
      comparison_summary,
      by = "Gene"
    ) |>
    dplyr::left_join(
      gene_effect_summary,
      by = "Gene"
    )
  
  #########################################################################
  # 8. Biological recurrence categories
  #########################################################################
  
  total_celltypes <-
    ncol(occurrence_numeric)
  
  gene_summary <-
    gene_summary |>
    dplyr::mutate(
      
      Category =
        dplyr::case_when(
          
          CellType_Frequency == 1 ~
            "CellType_Specific",
          
          CellType_Frequency == 2 ~
            "Shared_2_CellTypes",
          
          CellType_Frequency >= 3 &
            CellType_Frequency < total_celltypes ~
            "Highly_Shared",
          
          CellType_Frequency == total_celltypes ~
            "Conserved_All_CellTypes",
          
          TRUE ~
            "Other"
        )
    )
  
  #########################################################################
  # 9. Export complete analytical matrices
  #########################################################################
  
  save_csv(
    occurrence_matrix,
    "5.4_Gene_CellType_Occurrence_Matrix.csv"
  )
  
  save_csv(
    comparison_matrix,
    "5.4_Gene_Comparison_Occurrence_Matrix.csv"
  )
  
  save_csv(
    gene_summary,
    "5.4_Gene_Occurrence_Summary.csv"
  )
  
  #########################################################################
  # 10. Shared, unique and conserved genes
  #########################################################################
  
  shared_genes <-
    gene_summary |>
    dplyr::filter(
      CellType_Frequency >= 2
    ) |>
    dplyr::arrange(
      dplyr::desc(CellType_Frequency),
      dplyr::desc(Comparison_Frequency),
      dplyr::desc(Mean_Abs_Log2FC)
    )
  
  unique_genes <-
    gene_summary |>
    dplyr::filter(
      CellType_Frequency == 1
    ) |>
    dplyr::arrange(
      dplyr::desc(Mean_Abs_Log2FC)
    )
  
  conserved_genes <-
    gene_summary |>
    dplyr::filter(
      CellType_Frequency == total_celltypes
    ) |>
    dplyr::arrange(
      dplyr::desc(Mean_Abs_Log2FC)
    )
  
  save_csv(
    shared_genes,
    "5.4_Shared_Genes.csv"
  )
  
  save_csv(
    unique_genes,
    "5.4_Unique_Genes.csv"
  )
  
  save_csv(
    conserved_genes,
    "5.4_Conserved_Genes.csv"
  )
  
  #########################################################################
  # 11. Category statistics
  #########################################################################
  
  category_summary <-
    gene_summary |>
    dplyr::count(
      Category,
      name = "Number_of_Genes"
    ) |>
    dplyr::arrange(
      dplyr::desc(Number_of_Genes)
    )
  
  save_csv(
    category_summary,
    "5.4_Gene_Category_Summary.csv"
  )
  
  #########################################################################
  # 12. Select informative genes for recurrence heatmap
  #
  # Full matrix is retained above.
  # Heatmap intentionally displays only the strongest recurrent genes.
  #########################################################################
  
  TOP_RECURRENCE_GENES <- 50
  
  recurrence_rank <-
    gene_summary |>
    dplyr::filter(
      CellType_Frequency >= 2
    ) |>
    dplyr::arrange(
      dplyr::desc(CellType_Frequency),
      dplyr::desc(Comparison_Frequency),
      dplyr::desc(Mean_Abs_Log2FC),
      Min_Adjusted_P
    )
  
  if (nrow(recurrence_rank) == 0) {
    
    stop(
      "Section 5.4: No recurrent significant genes were identified."
    )
    
  }
  
  top_recurrence_genes <-
    recurrence_rank |>
    dplyr::slice_head(
      n = min(
        TOP_RECURRENCE_GENES,
        nrow(recurrence_rank)
      )
    )
  
  save_csv(
    top_recurrence_genes,
    "5.4_Top_Recurrence_Genes.csv"
  )
  
  #########################################################################
  # 13. Build heatmap matrix
  #########################################################################
  
  recurrence_heat_matrix <-
    occurrence_numeric[
      top_recurrence_genes$Gene,
      ,
      drop = FALSE
    ]
  
  # Remove rows without variation.
  #
  # A gene present in every cell type is biologically important,
  # but a completely constant row contributes no clustering information.
  # Conserved genes remain available in the exported tables.
  
  variable_rows <-
    apply(
      recurrence_heat_matrix,
      1,
      function(x) {
        length(unique(x)) > 1
      }
    )
  
  recurrence_heat_matrix <-
    recurrence_heat_matrix[
      variable_rows,
      ,
      drop = FALSE
    ]
  
  if (nrow(recurrence_heat_matrix) == 0) {
    
    stop(
      paste0(
        "Section 5.4: Selected recurrence genes are constant across ",
        "all cell types; no informative recurrence heatmap can be drawn."
      )
    )
    
  }
  
  cat(
    "Genes displayed in recurrence heatmap :",
    nrow(recurrence_heat_matrix),
    "\n"
  )
  
  cat(
    "Cell types displayed                  :",
    ncol(recurrence_heat_matrix),
    "\n"
  )
  
  #########################################################################
  # 14. Recurrence heatmap colours
  #########################################################################
  
  recurrence_col_fun <-
    circlize::colorRamp2(
      c(0, 1),
      c(
        "#F7F7F7",
        "#2166AC"
      )
    )
  
  #########################################################################
  # 15. Gene recurrence heatmap
  #########################################################################
  
  recurrence_ht <-
    ComplexHeatmap::Heatmap(
      
      recurrence_heat_matrix,
      
      name = "Significant\nDEG",
      
      col = recurrence_col_fun,
      
      cluster_rows = TRUE,
      
      cluster_columns = TRUE,
      
      show_row_names = TRUE,
      
      show_column_names = TRUE,
      
      row_names_gp =
        grid::gpar(
          fontsize = 7
        ),
      
      column_names_gp =
        grid::gpar(
          fontsize = 9,
          fontface = "bold"
        ),
      
      column_names_rot = 45,
      
      row_title =
        "Recurrently Differentially Expressed Genes",
      
      column_title =
        "Immune Cell Types",
      
      heatmap_legend_param =
        list(
          
          title =
            "Significant\nDEG",
          
          at = c(0, 1),
          
          labels = c(
            "Absent",
            "Present"
          )
          
        ),
      
      border = TRUE,
      
      use_raster = FALSE
    )
  
  #########################################################################
  # 16. Export recurrence heatmap - PDF
  #########################################################################
  
  pdf(
    file.path(
      FIGURE_DIR,
      "5.4_Gene_Recurrence_Heatmap.pdf"
    ),
    width = 9,
    height = 12,
    onefile = FALSE
  )
  
  ComplexHeatmap::draw(
    recurrence_ht,
    heatmap_legend_side = "right"
  )
  
  dev.off()
  
  #########################################################################
  # 17. Export recurrence heatmap - PNG
  #########################################################################
  
  png(
    file.path(
      FIGURE_DIR,
      "5.4_Gene_Recurrence_Heatmap.png"
    ),
    width = 2700,
    height = 3600,
    res = 300,
    type = ifelse(
      capabilities("cairo"),
      "cairo",
      "windows"
    )
  )
  
  ComplexHeatmap::draw(
    recurrence_ht,
    heatmap_legend_side = "right"
  )
  
  dev.off()
  
  #########################################################################
  # 18. Recurrence distribution plot
  #########################################################################
  
  recurrence_distribution <-
    gene_summary |>
    dplyr::count(
      CellType_Frequency,
      name = "Genes"
    )
  
  p_recurrence <-
    ggplot(
      recurrence_distribution,
      aes(
        x = CellType_Frequency,
        y = Genes
      )
    ) +
    geom_col(
      fill = "#2166AC",
      width = 0.75
    ) +
    geom_text(
      aes(
        label = Genes
      ),
      vjust = -0.3,
      size = 3.5
    ) +
    scale_x_continuous(
      breaks =
        seq_len(total_celltypes)
    ) +
    labs(
      title =
        "Recurrence of Significant DEGs Across Cell Types",
      x =
        "Number of Cell Types",
      y =
        "Number of Genes"
    ) +
    theme_classic(base_size = 12)
  
  save_plot(
    p_recurrence,
    "5.4_Gene_Recurrence_Distribution"
  )
  
  #########################################################################
  # 19. Report
  #########################################################################
  
  report <-
    data.frame(
      
      Metric = c(
        "Significant DEG Entries",
        "Total Significant Genes",
        "Shared Genes",
        "Cell-Type-Specific Genes",
        "Conserved Genes",
        "Genes Displayed in Recurrence Heatmap"
      ),
      
      Value = c(
        nrow(significant_deg),
        nrow(gene_summary),
        nrow(shared_genes),
        nrow(unique_genes),
        nrow(conserved_genes),
        nrow(recurrence_heat_matrix)
      ),
      
      stringsAsFactors = FALSE
    )
  
  save_csv(
    report,
    "5.4_Report.csv"
  )
  
  #########################################################################
  # 20. Console summary
  #########################################################################
  
  cat("\n=========================================\n")
  cat(" SHARED & UNIQUE DEG ANALYSIS COMPLETED\n")
  cat("=========================================\n\n")
  
  cat(
    "Significant genes       :",
    nrow(gene_summary),
    "\n"
  )
  
  cat(
    "Shared genes            :",
    nrow(shared_genes),
    "\n"
  )
  
  cat(
    "Cell-type-specific genes:",
    nrow(unique_genes),
    "\n"
  )
  
  cat(
    "Conserved genes         :",
    nrow(conserved_genes),
    "\n"
  )
  
  cat(
    "Heatmap genes           :",
    nrow(recurrence_heat_matrix),
    "\n"
  )
  
  log_message(
    "Section 5.4 completed successfully."
  )
  
})

###############################################################################
# SECTION 5.5 : CONSERVED DISEASE SIGNATURES
###############################################################################

section_timer("SECTION 5.5 : Conserved Disease Signatures", {
  
  cat("\n=================================================\n")
  cat(" SECTION 5.5 : CONSERVED DISEASE SIGNATURES\n")
  cat("=================================================\n\n")
  
  log_message("Identifying conserved disease signatures.")
  
  #########################################################################
  # Load occurrence summary
  #########################################################################
  
  gene_summary <- read.csv(
    file.path(TABLE_DIR, "5.4_Gene_Occurrence_Summary.csv"),
    stringsAsFactors = FALSE
  )
  
  #########################################################################
  # Gene-level statistics
  #########################################################################
  
  gene_stats <-
    master_deg |>
    dplyr::group_by(Gene) |>
    dplyr::summarise(
      Mean_Log2FC = mean(avg_log2FC, na.rm = TRUE),
      Mean_Abs_Log2FC = mean(abs(avg_log2FC), na.rm = TRUE),
      Median_Log2FC = median(avg_log2FC, na.rm = TRUE),
      Mean_Adjusted_P = mean(p_val_adj, na.rm = TRUE),
      CellType_Frequency = dplyr::n_distinct(CellType),
      Comparison_Frequency = dplyr::n_distinct(Comparison),
      .groups = "drop"
    )
  
  #########################################################################
  # Merge occurrence information
  #########################################################################
  
  signature_table <-
    dplyr::left_join(
      gene_stats,
      gene_summary |>
        dplyr::select(Gene, Category),
      by = "Gene"
    )
  
  #########################################################################
  # Composite score
  #########################################################################
  
  signature_table <-
    signature_table |>
    dplyr::mutate(
      Signature_Score =
        Mean_Abs_Log2FC *
        CellType_Frequency *
        (-log10(Mean_Adjusted_P + 1e-300))
    )
  
  #########################################################################
  # Conserved signatures
  #########################################################################
  
  conserved_signature <-
    signature_table |>
    dplyr::filter(
      Mean_Adjusted_P < ALPHA,
      CellType_Frequency >= 3
    ) |>
    dplyr::arrange(
      dplyr::desc(Signature_Score)
    )
  
  top_signature <-
    head(conserved_signature, 50)
  
  #########################################################################
  # Export tables
  #########################################################################
  
  save_csv(
    signature_table,
    "5.5_All_Gene_Signatures.csv"
  )
  
  save_csv(
    conserved_signature,
    "5.5_Conserved_Disease_Signatures.csv"
  )
  
  save_csv(
    top_signature,
    "5.5_Top50_Conserved_Signatures.csv"
  )
  
  #########################################################################
  # Heatmap
  #########################################################################
  
  heat_df <-
    master_deg |>
    dplyr::filter(
      Gene %in% top_signature$Gene
    ) |>
    dplyr::group_by(Gene, CellType) |>
    dplyr::summarise(
      Mean_Log2FC = mean(avg_log2FC, na.rm = TRUE),
      .groups = "drop"
    )
  
  heat_matrix <-
    reshape2::dcast(
      heat_df,
      Gene ~ CellType,
      value.var = "Mean_Log2FC",
      fill = 0
    )
  
  mat <- as.matrix(heat_matrix[, -1])
  rownames(mat) <- heat_matrix$Gene
  
  col_fun <-
    circlize::colorRamp2(
      c(-2, 0, 2),
      c("#4575B4", "white", "#D73027")
    )
  
  ht <-
    ComplexHeatmap::Heatmap(
      mat,
      name = "log2FC",
      col = col_fun,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 7),
      column_title = "Cell Types",
      row_title = "Top Conserved Signatures"
    )
  
  pdf(
    file.path(
      FIGURE_DIR,
      "5.5_Conserved_Signature_Heatmap.pdf"
    ),
    width = 9,
    height = 10
  )
  
  ComplexHeatmap::draw(ht)
  dev.off()
  
  png(
    file.path(
      FIGURE_DIR,
      "5.5_Conserved_Signature_Heatmap.png"
    ),
    width = 2600,
    height = 3000,
    res = 300
  )
  
  ComplexHeatmap::draw(ht)
  dev.off()
  
  #########################################################################
  # Report
  #########################################################################
  
  report <- data.frame(
    Metric = c(
      "Genes Evaluated",
      "Conserved Signatures",
      "Top Signature Genes"
    ),
    Value = c(
      nrow(signature_table),
      nrow(conserved_signature),
      nrow(top_signature)
    ),
    stringsAsFactors = FALSE
  )
  
  save_csv(
    report,
    "5.5_Conserved_Signature_Report.csv"
  )
  
  log_message("Section 5.5 completed successfully.")
  
})

###############################################################################
# SECTION 5.6 : PUBLICATION FIGURE ASSEMBLY
###############################################################################

section_timer("SECTION 5.6 : Publication Figure Assembly", {
  
  log_message("Assembling publication-ready composite figures.")
  
  suppressPackageStartupMessages({
    library(cowplot)
    library(magick)
  })
  
  img <- function(name){
    pngf <- file.path(FIGURE_DIR, paste0(name, ".png"))
    if(!file.exists(pngf)) stop("Missing figure: ", pngf)
    magick::image_read(pngf)
  }
  
  export_panel <- function(panel, base){
    magick::image_write(panel,
                        path=file.path(FIGURE_DIR,paste0(base,".png")),
                        format="png")
    magick::image_write(panel,
                        path=file.path(FIGURE_DIR,paste0(base,".tiff")),
                        format="tiff")
    magick::image_write(panel,
                        path=file.path(FIGURE_DIR,paste0(base,".pdf")),
                        format="pdf")
    magick::image_write(panel,
                        path=file.path(FIGURE_DIR,paste0(base,".svg")),
                        format="svg")
  }
  
  label_panel <- function(im, lab){
    magick::image_annotate(
      im,
      text = lab,
      size = 80,
      gravity = "northwest",
      location = "+35+25",
      font = "Arial",
      weight = 700
    )
  }
  
  #######################################################################
  # Figure 1
  #######################################################################
  
  f1a <- label_panel(img("5.3_CellType_DEG_Count"),"A")
  f1b <- label_panel(img("5.3_Comparison_DEG_Count"),"B")
  f1c <- label_panel(img("5.3_Up_vs_Down"),"C")
  f1d <- label_panel(img("5.3_Significant_DEG_Heatmap"),"D")
  
  row1 <- image_append(c(f1a,f1b))
  row2 <- image_append(c(f1c,f1d))
  fig1 <- image_append(c(row1,row2),stack=TRUE)
  
  export_panel(fig1,"Figure_1_Global_DEG_Landscape")
  
  #######################################################################
  # Figure 2
  #######################################################################
  
  f2a <- label_panel(img("5.4_Gene_Recurrence_Heatmap"),"A")
  f2b <- label_panel(img("5.5_Conserved_Signature_Heatmap"),"B")
  f2c <- magick::image_blank(1600,1200,"white") |>
    image_annotate("Shared & Conserved Gene Statistics",
                   gravity="center",
                   size=60)
  
  f2c <- label_panel(f2c,"C")
  
  fig2 <- image_append(
    c(
      image_append(c(f2a,f2b)),
      image_append(c(f2c,magick::image_blank(1600,1200,"white")))
    ),
    stack=TRUE
  )
  
  export_panel(fig2,"Figure_2_Shared_Disease_Signatures")
  
  #######################################################################
  # Figure Manifest
  #######################################################################
  
  manifest <- data.frame(
    Figure=c(
      "Figure 1",
      "Figure 2"
    ),
    Title=c(
      "Global DEG Landscape",
      "Shared Disease Signatures"
    ),
    Panels=c(
      "A-D",
      "A-C"
    ),
    Output=c(
      "PNG/PDF/TIFF/SVG",
      "PNG/PDF/TIFF/SVG"
    ),
    stringsAsFactors=FALSE
  )
  
  save_csv(
    manifest,
    "5.6_Publication_Figure_Index.csv"
  )
  
  log_message("Publication figures assembled successfully.")
  
})


###############################################################################
# SECTION 5.7 : AUTOMATIC FIGURE LEGENDS
###############################################################################

section_timer("SECTION 5.7 : Automatic Figure Legends", {
  
  log_message("Generating publication-ready figure legends.")
  
  figures <- data.frame(
    Figure=c("Figure 1","Figure 2"),
    File=c("Figure_1_Global_DEG_Landscape",
           "Figure_2_Shared_Disease_Signatures"),
    Legend=c(
      "Figure 1. Global differential expression landscape across all analyzed cell types. (A) Total number of differentially expressed genes identified for each cell type. (B) DEG distribution across disease comparisons. (C) Upregulated and downregulated DEG counts by cell type. (D) Heatmap showing the distribution of significant DEGs across cell types and comparisons. Differential expression significance was determined using the thresholds defined in the analysis pipeline.",
      "Figure 2. Shared disease signatures identified across cell populations. (A) Gene recurrence heatmap illustrating the presence of DEGs across cell types. (B) Heatmap of the top conserved disease signatures ranked by composite signature score. (C) Summary panel describing shared and conserved genes identified across the dataset."
    ),
    stringsAsFactors=FALSE
  )
  
  save_csv(figures,
           "5.7_Figure_Legends.csv")
  
  txt <- c(
    "###############################################################################",
    "# FIGURE LEGENDS",
    "###############################################################################",
    "",
    paste0(figures$Figure, ". ", figures$Legend, collapse="\n\n")
  )
  
  writeLines(
    txt,
    con=file.path(LEGEND_DIR,
                  "Figure_Legends.txt")
  )
  
  if(requireNamespace("officer", quietly=TRUE)){
    doc <- officer::read_docx()
    doc <- officer::body_add_par(
      doc,
      "Figure Legends",
      style="heading 1"
    )
    
    for(i in seq_len(nrow(figures))){
      doc <- officer::body_add_par(
        doc,
        figures$Figure[i],
        style="heading 2"
      )
      doc <- officer::body_add_par(
        doc,
        figures$Legend[i],
        style="Normal"
      )
    }
    
    print(
      doc,
      target=file.path(
        LEGEND_DIR,
        "Figure_Legends.docx"
      )
    )
  }
  
  manifest <- data.frame(
    Output=c(
      "5.7_Figure_Legends.csv",
      "Figure_Legends.txt",
      "Figure_Legends.docx"
    ),
    Status=c(
      "Generated",
      "Generated",
      ifelse(requireNamespace("officer", quietly=TRUE),
             "Generated","Skipped")
    ),
    stringsAsFactors=FALSE
  )
  
  save_csv(
    manifest,
    "5.7_Legend_Output_Index.csv"
  )
  
  log_message("Automatic figure legends completed.")
  
})
