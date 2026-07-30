###############################################################
# SCRIPT 21: FUNCTIONAL ENRICHMENT ANALYSIS
# PART 1 - INITIALIZATION & DATA PREPARATION
###############################################################

rm(list = ls())
gc()

###############################################################
# LOAD LIBRARIES
###############################################################

library(clusterProfiler)
library(enrichplot)
library(DOSE)

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

###############################################################
# CREATE OUTPUT DIRECTORIES
###############################################################

base.dir <- "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS/Functional_Enrichment"

dir.list <- c(
  
  "GO_BP/Tables",
  "GO_BP/Dotplots",
  "GO_BP/Barplots",
  
  "GO_CC/Tables",
  "GO_CC/Dotplots",
  "GO_CC/Barplots",
  
  "GO_MF/Tables",
  "GO_MF/Dotplots",
  "GO_MF/Barplots",
  
  "KEGG/Tables",
  "KEGG/Dotplots",
  "KEGG/Barplots",
  
  "Master_Summary"
  
)

for(i in dir.list){
  
  dir.create(
    file.path(base.dir, i),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
}

###############################################################
# DEG FILES
###############################################################

deg.folder <-
  "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS/DEG"

deg.files <- list.files(
  
  path = deg.folder,
  
  pattern = "\\.csv$",
  
  full.names = TRUE
  
)

cat("\n=========================================\n")
cat("DEG files detected :", length(deg.files), "\n")
cat("=========================================\n")

###############################################################
# MASTER SUMMARY TABLES
###############################################################

GO_BP_summary <- data.frame()

GO_CC_summary <- data.frame()

GO_MF_summary <- data.frame()

KEGG_summary <- data.frame()

###############################################################
# PREPARE DEG FUNCTION
###############################################################

prepare_deg <- function(file){
  
  cat("\n----------------------------------------\n")
  cat("Reading :", basename(file), "\n")
  
  deg <- read_csv(
    
    file,
    
    show_col_types = FALSE
    
  )
  
  #############################################################
  # CHECK REQUIRED COLUMNS
  #############################################################
  
  required.columns <- c(
    
    "Gene",
    
    "avg_log2FC",
    
    "p_val_adj"
    
  )
  
  missing.columns <-
    setdiff(required.columns, colnames(deg))
  
  if(length(missing.columns) > 0){
    
    stop(
      
      paste(
        
        "Missing columns:",
        
        paste(missing.columns,
              collapse = ", ")
        
      )
      
    )
    
  }
  
  #############################################################
  # FILTER SIGNIFICANT DEGs
  #############################################################
  
  deg.sig <-
    
    deg %>%
    
    filter(
      
      p_val_adj < 0.05,
      
      abs(avg_log2FC) >= 0.5
      
    )
  
  #############################################################
  # SPLIT UP & DOWN
  #############################################################
  
  up.genes <-
    
    deg.sig %>%
    
    filter(avg_log2FC > 0) %>%
    
    pull(Gene) %>%
    
    unique()
  
  down.genes <-
    
    deg.sig %>%
    
    filter(avg_log2FC < 0) %>%
    
    pull(Gene) %>%
    
    unique()
  
  #############################################################
  # PRINT SUMMARY
  #############################################################
  
  cat("Significant DEGs :", nrow(deg.sig), "\n")
  
  cat("Upregulated      :", length(up.genes), "\n")
  
  cat("Downregulated    :", length(down.genes), "\n")
  
  #############################################################
  # RETURN
  #############################################################
  
  return(
    
    list(
      
      filename = tools::file_path_sans_ext(
        
        basename(file)
        
      ),
      
      up = up.genes,
      
      down = down.genes,
      
      deg = deg.sig
      
    )
    
  )
  
}

###############################################################
# READ ALL DEG FILES
###############################################################

deg.list <-
  
  lapply(
    
    deg.files,
    
    prepare_deg
    
  )

names(deg.list) <-
  
  sapply(
    
    deg.list,
    
    function(x) x$filename
    
  )

###############################################################
# FINAL SUMMARY
###############################################################

cat("\n=========================================\n")
cat("Part 1 completed successfully.\n")
cat("All DEG files prepared.\n")
cat("Ready for functional enrichment.\n")
cat("=========================================\n")

###############################################################
# PART 2
# RUN g:PROFILER ENRICHMENT
###############################################################

###############################################################
# Create output folders
###############################################################

db.list <- c(
  "GO_BP",
  "GO_CC",
  "GO_MF",
  "KEGG",
  "Reactome",
  "WikiPathways"
)

for(db in db.list){
  
  dir.create(
    file.path(base.dir, db, "Tables"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
}

###############################################################
# Master result object
###############################################################

all.results <- list()

###############################################################
# Convert list columns before writing CSV
###############################################################

flatten_lists <- function(df){
  
  is.list.col <- sapply(df, is.list)
  
  if(any(is.list.col)){
    
    df[is.list.col] <- lapply(
      
      df[is.list.col],
      
      function(col){
        
        sapply(col, function(x){
          
          if(length(x) == 0){
            
            NA_character_
            
          } else{
            
            paste(x, collapse = ";")
            
          }
          
        })
        
      }
      
    )
    
  }
  
  df
  
}

###############################################################
# Run enrichment
###############################################################

for(comp in names(deg.list)){
  
  cat("\n=====================================\n")
  cat("Processing:", comp, "\n")
  cat("=====================================\n")
  
  current <- deg.list[[comp]]
  
  gene.list <- list(
    
    Up = current$up,
    
    Down = current$down
    
  )
  
  for(direction in names(gene.list)){
    
    genes <- unique(gene.list[[direction]])
    
    cat(direction, "genes:", length(genes), "\n")
    
    if(length(genes) < 5){
      
      cat("Skipped (<5 genes)\n")
      
      next
      
    }
    
    enrich <- tryCatch(
      
      gost(
        
        query = genes,
        
        organism = "hsapiens",
        
        sources = c(
          
          "GO:BP",
          "GO:CC",
          "GO:MF",
          "KEGG",
          "REAC",
          "WP"
          
        ),
        
        correction_method = "fdr",
        
        significant = TRUE
        
      ),
      
      error=function(e){
        
        cat("Error:", e$message, "\n")
        
        return(NULL)
        
      }
      
    )
    
    if(is.null(enrich))
      next
    
    if(is.null(enrich$result))
      next
    
    if(nrow(enrich$result)==0){
      
      cat("No enriched terms\n")
      
      next
      
    }
    
    res <- enrich$result
    
    res$Comparison <- comp
    
    res$Direction <- direction
    
    all.results[[paste(comp, direction, sep = "_")]] <- res
    
    cat("Enriched terms:", nrow(res), "\n")
    cat("Enriched terms:", nrow(res), "\n")
    
    ###########################################################
    # Split databases
    ###########################################################
    
    go.bp <- res %>%
      filter(source == "GO:BP")
    
    go.cc <- res %>%
      filter(source == "GO:CC")
    
    go.mf <- res %>%
      filter(source == "GO:MF")
    
    kegg <- res %>%
      filter(source == "KEGG")
    
    reactome <- res %>%
      filter(source == "REAC")
    
    wiki <- res %>%
      filter(source == "WP")
    
    ###########################################################
    # Save Tables
    ###########################################################
    
    if(nrow(go.bp)>0){
      
      go.bp <- flatten_lists(go.bp)
      
      write.csv(
        
        go.bp,
        
        file.path(
          
          base.dir,
          
          "GO_BP",
          
          "Tables",
          
          paste0(comp,"_",direction,"_GO_BP.csv")
          
        ),
        
        row.names=FALSE
        
      )
      
    }
    
    if(nrow(go.cc)>0){
      
      go.cc <- flatten_lists(go.cc)
      
      write.csv(
        
        go.cc,
        
        file.path(
          
          base.dir,
          
          "GO_CC",
          
          "Tables",
          
          paste0(comp,"_",direction,"_GO_CC.csv")
          
        ),
        
        row.names=FALSE
        
      )
      
    }
    
    if(nrow(go.mf)>0){
      
      go.mf <- flatten_lists(go.mf)
      
       write.csv(
        
        go.mf,
        
        file.path(
          
          base.dir,
          
          "GO_MF",
          
          "Tables",
          
          paste0(comp,"_",direction,"_GO_MF.csv")
          
        ),
        
        row.names=FALSE
        
      )
      
    }
    
    if(nrow(kegg)>0){
      
      kegg <- flatten_lists(kegg)
      
       write.csv(
        
        kegg,
        
        file.path(
          
          base.dir,
          
          "KEGG",
          
          "Tables",
          
          paste0(comp,"_",direction,"_KEGG.csv")
          
        ),
        
        row.names=FALSE
        
      )
      
    }
    
    if(nrow(reactome)>0){
      
      reactome <- flatten_lists(reactome)
      
      write.csv(
        
        reactome,
        
        file.path(
          
          base.dir,
          
          "Reactome",
          
          "Tables",
          
          paste0(comp,"_",direction,"_Reactome.csv")
          
        ),
        
        row.names=FALSE
        
      )
      
    }
    
    if(nrow(wiki)>0){
      
      wiki <- flatten_lists(wiki)
      
      write.csv(
        
        wiki,
        
        file.path(
          
          base.dir,
          
          "WikiPathways",
          
          "Tables",
          
          paste0(comp,"_",direction,"_WikiPathways.csv")
          
        ),
        
        row.names=FALSE
        
      )
      
    }
    
  }
  
}

###############################################################
# PART 2 COMPLETE
###############################################################

cat("\n=====================================\n")
cat("Functional enrichment completed.\n")
cat("Results stored in:\n")
cat(base.dir,"\n")
cat("=====================================\n")