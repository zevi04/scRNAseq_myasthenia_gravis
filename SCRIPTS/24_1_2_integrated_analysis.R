###############################################################################
# SCRIPT 24
# Integrated Biological Analysis
#
# Project  : Single-cell RNA-seq Analysis of Myasthenia Gravis
# Dataset  : GSE227835
# Platform : Seurat v5.5.1
# Author   : Zeno Vimalan A.
#
# Description:
# This script integrates all downstream analyses generated from the previous
# pipeline, performs biological interpretation, statistical summarization,
# cross-analysis, and generates publication-ready outputs.
#
# Version  : 1.0
# Date     : Sys.Date()
###############################################################################

rm(list = ls())
graphics.off()
gc()

options(stringsAsFactors = FALSE)
options(scipen = 999)

set.seed(1234)

cat("\n=========================================================\n")
cat(" Integrated Biological Analysis & Publication Engine\n")
cat(" Script 24\n")
cat("=========================================================\n\n")

###############################################################################
# SECTION 1 : LOAD REQUIRED PACKAGES
###############################################################################

required_cran <- c(
  "tidyverse",
  "data.table",
  "ggplot2",
  "ggpubr",
  "patchwork",
  "cowplot",
  "reshape2",
  "readr",
  "stringr",
  "forcats",
  "RColorBrewer",
  "viridis",
  "igraph",
  "ggraph",
  "scales",
  "rstatix",
  "broom",
  "glue",
  "fs"
)

required_bioc <- c(
  "ComplexHeatmap",
  "circlize"
)

install_if_missing <- function(pkg){
  if(!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

# Install CRAN packages
invisible(lapply(required_cran, install_if_missing))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

for(pkg in required_bioc){
  if(!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

# Combine package lists
all_packages <- c(required_cran, required_bioc)

# Load packages
invisible(lapply(all_packages, library, character.only = TRUE))

cat("✔ Packages loaded successfully.\n")

###############################################################################
# SECTION 2 : PROJECT DIRECTORIES
###############################################################################

PROJECT_DIR <- getwd()

RESULTS_DIR <- file.path(PROJECT_DIR, "RESULTS")

MASTER_DIR <- file.path(RESULTS_DIR, "Master_Summary")

SCRIPT24_DIR <- file.path(RESULTS_DIR,
                          "Script24_Integrated_Analysis")

FIGURE_DIR <- file.path(SCRIPT24_DIR,
                        "Publication_Figures")

TABLE_DIR <- file.path(SCRIPT24_DIR,
                       "Tables")

REPORT_DIR <- file.path(SCRIPT24_DIR,
                        "Reports")

LEGEND_DIR <- file.path(SCRIPT24_DIR,
                        "Figure_Legends")

METHOD_DIR <- file.path(SCRIPT24_DIR,
                        "Methods")

STATS_DIR <- file.path(SCRIPT24_DIR,
                       "Statistics")

META_DIR <- file.path(SCRIPT24_DIR,
                      "Metadata")

VALIDATION_DIR <- file.path(SCRIPT24_DIR,
                            "Validation")

dir_list <- c(
  
  SCRIPT24_DIR,
  
  FIGURE_DIR,
  
  TABLE_DIR,
  
  REPORT_DIR,
  
  LEGEND_DIR,
  
  METHOD_DIR,
  
  STATS_DIR,
  
  META_DIR,
  
  VALIDATION_DIR
  
)

invisible(lapply(dir_list,
                 dir.create,
                 recursive = TRUE,
                 showWarnings = FALSE))

cat("✔ Output folders created.\n")

###############################################################################
# SECTION 3 : GLOBAL SETTINGS
###############################################################################

FIGURE_WIDTH <- 10
FIGURE_HEIGHT <- 8

DPI <- 600

ALPHA <- 0.05

LOGFC_THRESHOLD <- 0.5

SIGNIFICANCE_LEVELS <- c(
  
  "ns",
  "*",
  "**",
  "***",
  "****"
  
)

theme_set(
  
  theme_bw(base_size = 14) +
    
    theme(
      
      plot.title = element_text(
        
        face = "bold",
        hjust = 0.5
        
      ),
      
      axis.title = element_text(face = "bold"),
      
      legend.title = element_text(face = "bold"),
      
      strip.text = element_text(face = "bold"),
      
      panel.grid.major = element_blank(),
      
      panel.grid.minor = element_blank()
      
    )
  
)

cat("✔ Global plotting theme initialized.\n")

###############################################################################
# SECTION 4 : COLOUR PALETTES
###############################################################################

GROUP_COLOURS <- c(
  
  Healthy = "#1B9E77",
  
  AChR_MG = "#D95F02",
  
  SNMG_Pre = "#7570B3",
  
  SNMG_Post = "#E7298A"
  
)

DIRECTION_COLOURS <- c(
  
  Up = "#D73027",
  
  Down = "#4575B4"
  
)

DATABASE_COLOURS <- c(
  
  GO_BP = "#1F78B4",
  
  GO_CC = "#33A02C",
  
  GO_MF = "#E31A1C",
  
  KEGG = "#FF7F00",
  
  Reactome = "#6A3D9A",
  
  WikiPathways = "#B15928"
  
)

cat("✔ Colour palettes loaded.\n")

###############################################################################
# SECTION 5 : LOGGING
###############################################################################

LOG_FILE <- file.path(
  
  REPORT_DIR,
  
  "Script24_Log.txt"
  
)

log_message <- function(message){
  
  timestamp <- format(Sys.time(),
                      "%Y-%m-%d %H:%M:%S")
  
  line <- paste0("[",
                 timestamp,
                 "] ",
                 message)
  
  cat(line, "\n")
  
  write(
    
    line,
    
    file = LOG_FILE,
    
    append = TRUE
    
  )
  
}

log_message("Script 24 started.")

###############################################################################
# SECTION 6 : TIMER FUNCTIONS
###############################################################################

section_timer <- function(section_name,
                          expression){
  
  cat("\n------------------------------------------\n")
  cat(section_name, "\n")
  cat("------------------------------------------\n")
  
  start <- Sys.time()
  
  result <- force(expression)
  
  end <- Sys.time()
  
  runtime <- round(
    
    as.numeric(
      
      difftime(end,
               start,
               units = "secs")
      
    ),
    
    2
    
  )
  
  cat("Completed in",
      runtime,
      "seconds.\n")
  
  log_message(
    
    paste(
      
      section_name,
      
      "-",
      
      runtime,
      
      "seconds"
      
    )
    
  )
  
  return(result)
  
}

###############################################################################
# SECTION 7 : HELPER FUNCTIONS
###############################################################################

save_csv <- function(df,
                     filename){
  
  write.csv(
    
    df,
    
    file.path(TABLE_DIR,
              filename),
    
    row.names = FALSE
    
  )
  
}

save_plot <- function(plot,
                      filename){
  
  ggsave(
    
    filename = file.path(
      
      FIGURE_DIR,
      
      paste0(filename,
             ".png")
      
    ),
    
    plot = plot,
    
    dpi = DPI,
    
    width = FIGURE_WIDTH,
    
    height = FIGURE_HEIGHT
    
  )
  
  ggsave(
    
    filename = file.path(
      
      FIGURE_DIR,
      
      paste0(filename,
             ".pdf")
      
    ),
    
    plot = plot,
    
    width = FIGURE_WIDTH,
    
    height = FIGURE_HEIGHT
    
  )
  
}

significance_label <- function(p){
  
  case_when(
    
    is.na(p) ~ "",
    
    p < 0.0001 ~ "****",
    
    p < 0.001 ~ "***",
    
    p < 0.01 ~ "**",
    
    p < 0.05 ~ "*",
    
    TRUE ~ "ns"
    
  )
  
}

###############################################################################
# SECTION 8 : STARTUP CHECK
###############################################################################

cat("\n=========================================\n")
cat(" Script 24 Foundation Loaded Successfully\n")
cat("=========================================\n")

log_message("Initialization complete.")

###############################################################################
# PROJECT ROOT
###############################################################################

PROJECT_DIR <- "C:/Bioinformatics/scRNAseq_myasthenia_gravis"

RESULTS_DIR <- file.path(PROJECT_DIR, "RESULTS")

cat("Project Directory :", PROJECT_DIR, "\n")
cat("Results Directory :", RESULTS_DIR, "\n")

###############################################################################
# SECTION 9 : PROJECT PATHS
###############################################################################

section_timer("SECTION 9 : Project Paths", {
  
  PATHS <- list(
    
    master       = file.path(RESULTS_DIR, "Master_Summary"),
    deg          = file.path(RESULTS_DIR, "DEG"),
    enrichment   = file.path(RESULTS_DIR, "Functional_Enrichment"),
    composition  = file.path(RESULTS_DIR, "Cell_Composition"),
    figures      = file.path(RESULTS_DIR, "figures"),
    qc           = file.path(RESULTS_DIR, "qc"),
    tables       = file.path(RESULTS_DIR, "tables"),
    metadata     = file.path(RESULTS_DIR, "metadata"),
    objects      = file.path(RESULTS_DIR, "objects"),
    logs         = file.path(RESULTS_DIR, "logs")
    
  )
  
  required_dirs <- names(PATHS)
  
  dir_validation <- data.frame(
    
    Directory = required_dirs,
    Path = unlist(PATHS),
    Exists = sapply(PATHS, dir.exists),
    stringsAsFactors = FALSE
    
  )
  
  write.csv(
    dir_validation,
    file.path(VALIDATION_DIR, "Directory_Validation.csv"),
    row.names = FALSE
  )
  
  if(any(!dir_validation$Exists)){
    
    stop(
      paste(
        "Missing directories:\n",
        paste(dir_validation$Directory[!dir_validation$Exists],
              collapse = "\n")
      )
    )
    
  }
  
  log_message("Directory validation completed.")
  
})

###############################################################################
# SECTION 10 : LOAD MASTER FILES
###############################################################################

section_timer("SECTION 10 : Loading Master Files", {
  
  MASTER_FILES <- c(
    
    Pipeline_Statistics = "Pipeline_Statistics.csv",
    DEG_Long = "DEG_Long.csv",
    DEG_Wide = "DEG_Wide.csv",
    CellType_Ranking = "CellType_Ranking.csv",
    Comparison_Ranking = "Comparison_Ranking.csv",
    Enrichment_Master = "Enrichment_Master.csv",
    Database_Ranking = "Database_Ranking.csv",
    CellType_Enrichment_Ranking = "CellType_Enrichment_Ranking.csv",
    Comparison_Enrichment_Ranking = "Comparison_Enrichment_Ranking.csv",
    Direction_Ranking = "Direction_Ranking.csv"
    
  )
  
  analysis_data <- list()
  
  validation_table <- data.frame()
  
  for(i in seq_along(MASTER_FILES)){
    
    file_name <- MASTER_FILES[i]
    
    full_path <- file.path(PATHS$master, file_name)
    
    if(!file.exists(full_path)){
      
      validation_table <- rbind(
        
        validation_table,
        
        data.frame(
          
          File=file_name,
          Status="NOT FOUND",
          Rows=NA,
          Columns=NA
          
        )
        
      )
      
      next
      
    }
    
    df <- read.csv(full_path)
    
    analysis_data[[names(MASTER_FILES)[i]]] <- df
    
    validation_table <- rbind(
      
      validation_table,
      
      data.frame(
        
        File=file_name,
        Status="PASS",
        Rows=nrow(df),
        Columns=ncol(df)
        
      )
      
    )
    
  }
  
  write.csv(
    validation_table,
    file.path(VALIDATION_DIR,
              "Master_File_Validation.csv"),
    row.names = FALSE
  )
  
})

###############################################################################
# SECTION 11 : LOAD ALL DEG TABLES
###############################################################################

section_timer("SECTION 11 : Loading DEG Tables", {
  
  deg_files <- list.files(
    
    PATHS$deg,
    
    pattern="\\.csv$",
    
    full.names=TRUE,
    
    recursive=TRUE
    
  )
  
  deg_tables <- lapply(deg_files, read.csv)
  
  names(deg_tables) <- tools::file_path_sans_ext(
    
    basename(deg_files)
    
  )
  
  log_message(
    
    paste(length(deg_tables),
          "DEG tables loaded.")
    
  )
  
})

###############################################################################
# SECTION 12 : LOAD ALL ENRICHMENT TABLES
###############################################################################

section_timer("SECTION 12 : Loading Enrichment Tables", {
  
  enrichment_files <- list.files(
    
    PATHS$enrichment,
    
    pattern="\\.csv$",
    
    recursive=TRUE,
    
    full.names=TRUE
    
  )
  
  enrichment_tables <- lapply(
    
    enrichment_files,
    
    read.csv
    
  )
  
  names(enrichment_tables) <-
    
    tools::file_path_sans_ext(
      
      basename(enrichment_files)
      
    )
  
  log_message(
    
    paste(length(enrichment_tables),
          "Enrichment tables loaded.")
    
  )
  
})

###############################################################################
# SECTION 13 : INVENTORY
###############################################################################

section_timer("SECTION 13 : Inventory", {
  
  inventory <- data.frame(
    
    File=c(
      
      basename(deg_files),
      
      basename(enrichment_files)
      
    ),
    
    Category=c(
      
      rep("DEG", length(deg_files)),
      
      rep("Enrichment",
          length(enrichment_files))
      
    ),
    
    Size_MB=round(
      
      c(
        
        file.info(deg_files)$size,
        
        file.info(enrichment_files)$size
        
      )/1024^2,
      
      3
      
    ),
    
    stringsAsFactors=FALSE
    
  )
  
  write.csv(
    
    inventory,
    
    file.path(
      
      VALIDATION_DIR,
      
      "Analysis_Inventory.csv"
      
    ),
    
    row.names=FALSE
    
  )
  
})

###############################################################################
# SECTION 14 : SUMMARY
###############################################################################

cat("\n========================================\n")
cat(" Script 24 Data Loaded Successfully\n")
cat("========================================\n")

cat("Master Tables      :", length(analysis_data), "\n")
cat("DEG Tables         :", length(deg_tables), "\n")
cat("Enrichment Tables  :", length(enrichment_tables), "\n")