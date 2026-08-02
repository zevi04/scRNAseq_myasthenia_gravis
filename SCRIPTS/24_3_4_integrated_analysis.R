###############################################################################
# CHUNK 3 : DATA INTEGRITY & BIOLOGICAL READINESS
###############################################################################

section_timer("CHUNK 3 : Data Integrity & Biological Readiness", {
  
  cat("\n==============================\n")
  cat(" Running Data Validation\n")
  cat("==============================\n")
  
  ###############################################################################
  # Required columns
  ###############################################################################
  
  required_deg <- c(
    "Gene",
    "p_val",
    "avg_log2FC",
    "pct.1",
    "pct.2",
    "p_val_adj"
  )
  
  required_enrichment <- c(
    "Comparison",
    "Direction",
    "source",
    "term_name",
    "p_value",
    "intersection_size"
  )
  
  ###############################################################################
  # Validation tables
  ###############################################################################
  
  deg_validation <- data.frame()
  
  for(nm in names(deg_tables)){
    
    df <- deg_tables[[nm]]
    
    missing_cols <- setdiff(required_deg, colnames(df))
    
    deg_validation <- rbind(
      deg_validation,
      data.frame(
        File = nm,
        Rows = nrow(df),
        Columns = ncol(df),
        MissingColumns = ifelse(length(missing_cols)==0,
                                "None",
                                paste(missing_cols,
                                      collapse=", ")),
        MissingValues = sum(is.na(df)),
        stringsAsFactors = FALSE
      )
    )
    
  }
  
  ###############################################################################
  # Enrichment validation
  ###############################################################################
  
  enrichment_validation <- data.frame()
  
  for(nm in names(enrichment_tables)){
    
    df <- enrichment_tables[[nm]]
    
    missing_cols <- setdiff(required_enrichment,
                            colnames(df))
    
    enrichment_validation <- rbind(
      
      enrichment_validation,
      
      data.frame(
        
        File = nm,
        
        Rows = nrow(df),
        
        Columns = ncol(df),
        
        MissingColumns =
          ifelse(length(missing_cols)==0,
                 "None",
                 paste(missing_cols,
                       collapse=", ")),
        
        MissingValues = sum(is.na(df)),
        
        stringsAsFactors = FALSE
        
      )
      
    )
    
  }
  
  ###############################################################################
  # Save validation reports
  ###############################################################################
  
  write.csv(
    deg_validation,
    file.path(VALIDATION_DIR,
              "DEG_Validation.csv"),
    row.names = FALSE
  )
  
  write.csv(
    enrichment_validation,
    file.path(VALIDATION_DIR,
              "Enrichment_Validation.csv"),
    row.names = FALSE
  )
  
  ###############################################################################
  # Critical failures
  ###############################################################################
  
  if(any(deg_validation$MissingColumns != "None")){
    
    stop("One or more DEG tables are missing required columns.")
    
  }
  
  if(any(enrichment_validation$MissingColumns != "None")){
    
    stop("One or more enrichment tables are missing required columns.")
    
  }
  
  
  ###############################################################################
  # Overall statistics
  ###############################################################################
  
  summary_report <- data.frame(
    
    Metric = c(
      "Master Tables",
      "DEG Tables",
      "Enrichment Tables"
    ),
    
    Value = c(
      length(analysis_data),
      length(deg_tables),
      length(enrichment_tables)
    ),
    
    stringsAsFactors = FALSE
    
  )
  
  write.csv(
    
    summary_report,
    
    file.path(
      VALIDATION_DIR,
      "Data_Integrity_Summary.csv"
    ),
    
    row.names = FALSE
    
  )
  
  ###############################################################################
  # Console output
  ###############################################################################
  
  cat("\n=====================================\n")
  cat(" DATA VALIDATION COMPLETED\n")
  cat("=====================================\n")
  
  cat("Master Tables       :", length(analysis_data), "\n")
  cat("DEG Tables          :", length(deg_tables), "\n")
  cat("Enrichment Tables   :", length(enrichment_tables), "\n")
  
  log_message("Chunk 3 completed successfully.")
  
})

###############################################################################
# CHUNK 4 : BUILD MASTER DEG TABLE
###############################################################################

section_timer("CHUNK 4 : DEG Metadata Integration", {
  
  cat("\n=========================================\n")
  cat(" Building Master DEG Table\n")
  cat("=========================================\n")
  
  ###############################################################################
  # Create metadata from filenames
  ###############################################################################
  
  master_deg <- list()
  
  for(i in seq_along(deg_tables)){
    
    df <- deg_tables[[i]]
    
    file_name <- names(deg_tables)[i]
    
    comparison <- sub(
      ".*?(Healthy_vs_AChR_MG|Healthy_vs_SNMG_Pre|SNMG_Pre_vs_SNMG_Post)$",
      "\\1",
      file_name
    )
    
    celltype <- sub(
      paste0("_", comparison, "$"),
      "",
      file_name
    )
    
    celltype <- gsub("__", "_", celltype)
    celltype <- gsub("_", " ", celltype)
    
    df$CellType <- celltype
    df$Comparison <- comparison
    
    df$Direction <- ifelse(
      df$avg_log2FC > 0,
      "Upregulated",
      "Downregulated"
    )
    
    master_deg[[i]] <- df
    
  }
  
  master_deg <- do.call(rbind, master_deg)
  
  rownames(master_deg) <- NULL
  
  ###############################################################################
  # Save master DEG table
  ###############################################################################
  
  write.csv(
    master_deg,
    file.path(PATHS$tables, "Master_DEG_Table.csv"),
    row.names = FALSE
  )
  
  ###############################################################################
  # Summary statistics
  ###############################################################################
  
  deg_summary <- aggregate(
    Gene ~ CellType + Comparison + Direction,
    data = master_deg,
    FUN = length
  )
  
  colnames(deg_summary)[4] <- "Number_of_DEGs"
  
  write.csv(
    deg_summary,
    file.path(PATHS$tables,
              "DEG_Summary.csv"),
    row.names = FALSE
  )
  
  ###############################################################################
  # Cell type ranking
  ###############################################################################
  
  celltype_rank <- aggregate(
    Gene ~ CellType,
    data = master_deg,
    FUN = length
  )
  
  colnames(celltype_rank)[2] <- "Total_DEGs"
  
  celltype_rank <- celltype_rank[
    order(-celltype_rank$Total_DEGs),
  ]
  
  write.csv(
    celltype_rank,
    file.path(PATHS$tables,
              "CellType_DEG_Ranking.csv"),
    row.names = FALSE
  )
  
  ###############################################################################
  # Comparison ranking
  ###############################################################################
  
  comparison_rank <- aggregate(
    Gene ~ Comparison,
    data = master_deg,
    FUN = length
  )
  
  colnames(comparison_rank)[2] <- "Total_DEGs"
  
  comparison_rank <- comparison_rank[
    order(-comparison_rank$Total_DEGs),
  ]
  
  write.csv(
    comparison_rank,
    file.path(PATHS$tables,
              "Comparison_DEG_Ranking.csv"),
    row.names = FALSE
  )
  
  ###############################################################################
  # Console output
  ###############################################################################
  
  cat("\n=========================================\n")
  cat(" DEG Metadata Integration Complete\n")
  cat("=========================================\n")
  
  cat("Total DEG entries :", nrow(master_deg), "\n")
  cat("Unique genes      :", length(unique(master_deg$Gene)), "\n")
  cat("Cell types        :", length(unique(master_deg$CellType)), "\n")
  cat("Comparisons       :", length(unique(master_deg$Comparison)), "\n")
  
  log_message("Chunk 4 completed successfully.")
  
})