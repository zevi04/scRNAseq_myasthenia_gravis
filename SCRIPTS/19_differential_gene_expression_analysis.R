###############################################################
# SCRIPT 19: DIFFERENTIAL GENE EXPRESSION ANALYSIS
###############################################################

rm(list = ls())
gc()

###############################################################
# Load Libraries
###############################################################

library(Seurat)
library(dplyr)

###############################################################
# Create Output Folder
###############################################################

dir.create(
  "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS/DEG",
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
# Create Disease Metadata
###############################################################

scaled_seurat$Condition <- case_when(
  grepl("^H", scaled_seurat$orig.ident) ~ "Healthy",
  grepl("^A", scaled_seurat$orig.ident) ~ "AChR_MG",
  grepl("a$", scaled_seurat$orig.ident) ~ "SNMG_Pre",
  grepl("b$", scaled_seurat$orig.ident) ~ "SNMG_Post"
)

cat("\nCondition Distribution:\n")
print(table(scaled_seurat$Condition))

###############################################################
# Set Cell Type Identity
###############################################################

Idents(scaled_seurat) <- "SingleR_Label"

###############################################################
# Cell Types
###############################################################

cell.types <- levels(Idents(scaled_seurat))

###############################################################
# Comparisons
###############################################################

comparisons <- list(
  
  c("Healthy", "AChR_MG"),
  
  c("Healthy", "SNMG_Pre"),
  
  c("SNMG_Pre", "SNMG_Post")
  
)

###############################################################
# Differential Expression
###############################################################

for(cell in cell.types){
  
  cat("\n=========================================\n")
  cat("Cell Type:", cell, "\n")
  cat("=========================================\n")
  
  #############################################################
  # Subset One Cell Type
  #############################################################
  
  cell.object <- subset(
    scaled_seurat,
    idents = cell
  )
  
  #############################################################
  # Downsample Large Cell Types
  #############################################################
  
  set.seed(123)
  
  cell.object <- subset(
    cell.object,
    downsample = 3000
  )
  
  #############################################################
  # Set Identity to Disease Condition
  #############################################################
  
  Idents(cell.object) <- "Condition"
  
  #############################################################
  # Run Comparisons
  #############################################################
  
  for(comp in comparisons){
    
    group1 <- comp[1]
    group2 <- comp[2]
    
    n1 <- sum(Idents(cell.object) == group1)
    n2 <- sum(Idents(cell.object) == group2)
    
    cat("\n-----------------------------------------\n")
    cat(group1, "vs", group2, "\n")
    cat("Cells:", group1, "=", n1,
        "|", group2, "=", n2, "\n")
    
    if(n1 < 10 | n2 < 10){
      
      cat("Skipping (too few cells)\n")
      
      next
      
    }
    
    ###########################################################
    # Differential Expression
    ###########################################################
    
    deg <- FindMarkers(
      
      object = cell.object,
      
      ident.1 = group1,
      
      ident.2 = group2,
      
      test.use = "wilcox",
      
      min.pct = 0.10,
      
      logfc.threshold = 0.25
      
    )
    
    ###########################################################
    # Add Gene Names
    ###########################################################
    
    deg$Gene <- rownames(deg)
    
    deg <- deg %>%
      relocate(Gene)
    
    ###########################################################
    # Save Results
    ###########################################################
    
    filename <- paste0(
      
      "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS/DEG/",
      
      gsub("[ +]", "_", cell),
      
      "_",
      
      group1,
      
      "_vs_",
      
      group2,
      
      ".csv"
      
    )
    
    write.csv(
      deg,
      filename,
      row.names = FALSE
    )
    
    cat("DEGs Identified:", nrow(deg), "\n")
    
  }
  
}

###############################################################
# Finished
###############################################################

cat("\n=========================================\n")
cat("Differential Expression Analysis Complete\n")
cat("=========================================\n")