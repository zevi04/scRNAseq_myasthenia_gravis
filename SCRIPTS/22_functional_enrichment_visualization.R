###############################################################
# SCRIPT 22
# Functional Enrichment Visualization
###############################################################

library(ggplot2)
library(dplyr)

base.dir <- "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS/Functional_Enrichment"

databases <- c(
  "GO_BP",
  "GO_CC",
  "GO_MF",
  "KEGG",
  "Reactome",
  "WikiPathways"
)

###############################################################
# LOOP THROUGH DATABASES
###############################################################

for(db in databases){
  
  cat("\nProcessing:", db,"\n")
  
  table.dir <- file.path(base.dir,db,"Tables")
  dot.dir   <- file.path(base.dir,db,"Dotplots")
  bar.dir   <- file.path(base.dir,db,"Barplots")
  
  dir.create(dot.dir,showWarnings = FALSE,recursive = TRUE)
  dir.create(bar.dir,showWarnings = FALSE,recursive = TRUE)
  
  files <- list.files(
    table.dir,
    pattern="\\.csv$",
    full.names=TRUE
  )
  
  for(f in files){
    
    df <- read.csv(f)
    
    if(nrow(df)==0) next
    
    ###########################################################
    # PREPARE DATA
    ###########################################################
    
    df <- df %>%
      arrange(p_value)
    
    n.plot <- min(20,nrow(df))
    
    df <- df[1:n.plot,]
    
    df$GeneRatio <- df$intersection_size/df$query_size
    df$Score <- -log10(df$p_value)
    
    df$term_name <- factor(
      df$term_name,
      levels=rev(df$term_name)
    )
    
    name <- tools::file_path_sans_ext(basename(f))
    
    ###########################################################
    # DOTPLOT
    ###########################################################
    
    p1 <- ggplot(
      df,
      aes(
        x=GeneRatio,
        y=term_name,
        size=intersection_size,
        colour=Score
      )
    )+
      geom_point()+
      labs(
        title=name,
        x="Gene Ratio",
        y=NULL,
        colour="-log10(FDR)",
        size="Genes"
      )+
      theme_bw(base_size=12)
    
    ggsave(
      filename=file.path(dot.dir,paste0(name,"_Dotplot.png")),
      plot=p1,
      width=8,
      height=6,
      dpi=300
    )
    
    ###########################################################
    # BARPLOT
    ###########################################################
    
    p2 <- ggplot(
      df,
      aes(
        x=Score,
        y=term_name
      )
    )+
      geom_col()+
      labs(
        title=name,
        x="-log10(FDR)",
        y=NULL
      )+
      theme_bw(base_size=12)
    
    ggsave(
      filename=file.path(bar.dir,paste0(name,"_Barplot.png")),
      plot=p2,
      width=8,
      height=6,
      dpi=300
    )
    
    cat("✓",name,"\n")
    
  }
  
}

###############################################################
# COMPLETE
###############################################################

cat("\n=====================================\n")
cat("Visualization completed.\n")
cat("Plots saved successfully.\n")
cat("=====================================\n")