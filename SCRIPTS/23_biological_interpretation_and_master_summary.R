###############################################################
# SCRIPT 23
# PART 1
# PROJECT SCAN & PIPELINE STATISTICS
###############################################################

library(dplyr)
library(stringr)

###############################################################
# DEFINE DIRECTORIES
###############################################################

results.dir <- "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS"

summary.dir <- file.path(results.dir,"Master_Summary")

dir.create(summary.dir,
           recursive = TRUE,
           showWarnings = FALSE)

###############################################################
# LOCATE PROJECT FILES
###############################################################

cat("Scanning project...\n\n")

all.csv <- list.files(
  results.dir,
  pattern="\\.csv$",
  recursive=TRUE,
  full.names=TRUE
)

all.png <- list.files(
  results.dir,
  pattern="\\.png$",
  recursive=TRUE,
  full.names=TRUE
)

###############################################################
# DEG FILES
###############################################################

deg.files <- all.csv[
  grepl("Differential_Expression",all.csv)
]

###############################################################
# ENRICHMENT FILES
###############################################################

enrichment.files <- all.csv[
  grepl("Functional_Enrichment",all.csv)
]

###############################################################
# FIGURES
###############################################################

figure.files <- all.png

###############################################################
# CELL TYPES
###############################################################

deg.names <- basename(deg.files)

cell.types <- unique(
  str_extract(
    deg.names,
    "^[^_]+_[^_]+"
  )
)

###############################################################
# COMPARISONS
###############################################################

comparisons <- unique(
  
  str_extract(
    deg.names,
    "Healthy_vs_AChR_MG|
     Healthy_vs_SNMG_Pre|
     SNMG_Pre_vs_SNMG_Post"
  )
  
)

###############################################################
# DATABASES
###############################################################

databases <- unique(
  
  basename(
    dirname(
      dirname(
        enrichment.files
      )
    )
  )
  
)

###############################################################
# PIPELINE STATISTICS
###############################################################

pipeline.stats <- data.frame(
  
  Metric=c(
    
    "CSV Files",
    "PNG Figures",
    "DEG Files",
    "Enrichment Tables",
    "Cell Types",
    "Comparisons",
    "Databases"
    
  ),
  
  Value=c(
    
    length(all.csv),
    length(all.png),
    length(deg.files),
    length(enrichment.files),
    length(cell.types),
    length(comparisons),
    length(databases)
    
  )
  
)

###############################################################
# SAVE
###############################################################

write.csv(
  
  pipeline.stats,
  
  file.path(
    summary.dir,
    "Pipeline_Statistics.csv"
  ),
  
  row.names=FALSE
  
)

###############################################################
# PRINT SUMMARY
###############################################################

cat("=====================================\n")
cat("PROJECT SUMMARY\n")
cat("=====================================\n\n")

cat("CSV files              :",length(all.csv),"\n")
cat("PNG figures            :",length(all.png),"\n")
cat("DEG files              :",length(deg.files),"\n")
cat("Enrichment tables      :",length(enrichment.files),"\n")
cat("Cell types             :",length(cell.types),"\n")
cat("Comparisons            :",length(comparisons),"\n")
cat("Databases              :",length(databases),"\n\n")

cat("-------------------------------------\n")
cat("Cell Types Detected\n")
cat("-------------------------------------\n")

print(sort(cell.types))

cat("\n-------------------------------------\n")
cat("Comparisons\n")
cat("-------------------------------------\n")

print(sort(comparisons))

cat("\n-------------------------------------\n")
cat("Databases\n")
cat("-------------------------------------\n")

print(sort(databases))

cat("\n=====================================\n")
cat("PART 1 COMPLETE\n")
cat("=====================================\n")

###############################################################
# PART 2
# DEG MASTER SUMMARY
###############################################################

library(dplyr)
library(stringr)

###############################################################
# DIRECTORIES
###############################################################

results.dir <- "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS"
summary.dir <- file.path(results.dir,"Master_Summary")

deg.dir <- file.path(results.dir,"DEG")

###############################################################
# FIND DEG FILES
###############################################################

deg.files <- list.files(
  deg.dir,
  pattern="\\.csv$",
  recursive=TRUE,
  full.names=TRUE
)

###############################################################
# INITIALIZE
###############################################################

deg.summary <- data.frame()

###############################################################
# LOOP THROUGH FILES
###############################################################

for(file in deg.files){
  
  df <- read.csv(file)
  
  filename <- tools::file_path_sans_ext(basename(file))
  
  comparison <- str_extract(
    filename,
    "Healthy_vs_AChR_MG|Healthy_vs_SNMG_Pre|SNMG_Pre_vs_SNMG_Post"
  )
  
  cell.type <- filename
  
  cell.type <- str_remove(
    cell.type,
    "_Healthy_vs_AChR_MG|_Healthy_vs_SNMG_Pre|_SNMG_Pre_vs_SNMG_Post"
  )
  
  up <- sum(df$avg_log2FC > 0)
  
  down <- sum(df$avg_log2FC < 0)
  
  total <- nrow(df)
  
  mean.fc <- mean(abs(df$avg_log2FC))
  
  median.fc <- median(abs(df$avg_log2FC))
  
  volcano.score <- total * mean.fc
  
  deg.summary <- rbind(
    deg.summary,
    data.frame(
      Cell_Type=cell.type,
      Comparison=comparison,
      Upregulated=up,
      Downregulated=down,
      Total_DEGs=total,
      Mean_Log2FC=round(mean.fc,3),
      Median_Log2FC=round(median.fc,3),
      Volcano_Score=round(volcano.score,2)
    )
  )
  
}

###############################################################
# SAVE LONG TABLE
###############################################################

write.csv(
  deg.summary,
  file.path(summary.dir,"DEG_Long.csv"),
  row.names=FALSE
)

###############################################################
# CREATE WIDE TABLE
###############################################################

deg.wide <- reshape(
  deg.summary[,c("Cell_Type","Comparison","Total_DEGs")],
  idvar="Cell_Type",
  timevar="Comparison",
  direction="wide"
)

colnames(deg.wide) <- gsub(
  "Total_DEGs\\.",
  "",
  colnames(deg.wide)
)

deg.wide$Total_DEGs <- rowSums(
  deg.wide[,-1],
  na.rm=TRUE
)

deg.wide <- deg.wide[
  order(-deg.wide$Total_DEGs),
]

write.csv(
  deg.wide,
  file.path(summary.dir,"DEG_Wide.csv"),
  row.names=FALSE
)

###############################################################
# CELL TYPE RANKING
###############################################################

cell.rank <- deg.summary %>%
  group_by(Cell_Type) %>%
  summarise(
    Total_DEGs=sum(Total_DEGs),
    Mean_Volcano_Score=mean(Volcano_Score)
  ) %>%
  arrange(desc(Total_DEGs))

write.csv(
  cell.rank,
  file.path(summary.dir,"CellType_Ranking.csv"),
  row.names=FALSE
)

###############################################################
# COMPARISON RANKING
###############################################################

comparison.rank <- deg.summary %>%
  group_by(Comparison) %>%
  summarise(
    Total_DEGs=sum(Total_DEGs),
    Mean_Volcano_Score=mean(Volcano_Score)
  ) %>%
  arrange(desc(Total_DEGs))

write.csv(
  comparison.rank,
  file.path(summary.dir,"Comparison_Ranking.csv"),
  row.names=FALSE
)

###############################################################
# REPORT
###############################################################

sink(file.path(summary.dir,"DEG_Report.txt"))

cat("=========================================\n")
cat("DEG MASTER SUMMARY\n")
cat("=========================================\n\n")

cat("Total DEG files analysed :",length(deg.files),"\n")
cat("Total DEGs :",sum(deg.summary$Total_DEGs),"\n\n")

cat("Most affected cell type:\n")
print(cell.rank[1,])

cat("\nLeast affected cell type:\n")
print(cell.rank[nrow(cell.rank),])

cat("\nComparison ranking:\n")
print(comparison.rank)

sink()

###############################################################
# CONSOLE
###############################################################

cat("\n=========================================\n")
cat("PART 2 COMPLETE\n")
cat("=========================================\n")

cat("\nTop Cell Types:\n")
print(head(cell.rank))

cat("\nTop Comparisons:\n")
print(comparison.rank)

###############################################################
# PART 3
# ENRICHMENT MASTER SUMMARY
###############################################################

library(dplyr)
library(stringr)

###############################################################
# DIRECTORIES
###############################################################

results.dir <- "C:/Bioinformatics/scRNAseq_myasthenia_gravis/RESULTS"

summary.dir <- file.path(results.dir, "Master_Summary")

enrich.dir <- file.path(results.dir, "Functional_Enrichment")

if(!dir.exists(summary.dir)){
  dir.create(summary.dir, recursive = TRUE)
}

###############################################################
# FIND ALL ENRICHMENT TABLES
###############################################################

enrichment.files <- list.files(
  enrich.dir,
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

enrichment.files <- enrichment.files[
  grepl("Tables", enrichment.files)
]

cat("Enrichment tables found :", length(enrichment.files), "\n")

###############################################################
# INITIALIZE
###############################################################

master.enrichment <- data.frame()

###############################################################
# LOOP THROUGH TABLES
###############################################################

for(file in enrichment.files){
  
  df <- tryCatch(read.csv(file),
                 error = function(e) NULL)
  
  if(is.null(df)) next
  
  if(nrow(df) == 0) next
  
  #############################################################
  # FILE INFORMATION
  #############################################################
  
  filename <- tools::file_path_sans_ext(basename(file))
  
  database <- basename(dirname(dirname(file)))
  # Gives GO_BP, GO_CC, GO_MF, KEGG, Reactome, WikiPathways
  
  comparison <- str_extract(
    filename,
    "Healthy_vs_AChR_MG|Healthy_vs_SNMG_Pre|SNMG_Pre_vs_SNMG_Post"
  )
  
  cell.type <- str_remove(
    filename,
    "_Healthy_vs_AChR_MG|_Healthy_vs_SNMG_Pre|_SNMG_Pre_vs_SNMG_Post"
  )
  
  #############################################################
  # DIRECTION
  #############################################################
  
  direction <- if("Direction" %in% names(df))
    unique(df$Direction)[1]
  else
    NA
  
  #############################################################
  # SUMMARY STATISTICS
  #############################################################
  
  lowest.p <- min(df$p_value, na.rm = TRUE)
  
  top.term <- df$term_name[which.min(df$p_value)]
  
  mean.precision <- mean(df$precision, na.rm = TRUE)
  
  max.precision <- max(df$precision, na.rm = TRUE)
  
  mean.recall <- mean(df$recall, na.rm = TRUE)
  
  max.recall <- max(df$recall, na.rm = TRUE)
  
  largest.intersection <- max(df$intersection_size, na.rm = TRUE)
  
  #############################################################
  # ADD ROW
  #############################################################
  
  master.enrichment <- rbind(
    master.enrichment,
    data.frame(
      
      Cell_Type = cell.type,
      
      Comparison = comparison,
      
      Direction = direction,
      
      Database = database,
      
      Enriched_Terms = nrow(df),
      
      Top_Term = top.term,
      
      Lowest_P = lowest.p,
      
      Largest_Intersection = largest.intersection,
      
      Mean_Precision = round(mean.precision,3),
      
      Max_Precision = round(max.precision,3),
      
      Mean_Recall = round(mean.recall,3),
      
      Max_Recall = round(max.recall,3),
      
      stringsAsFactors = FALSE
      
    )
  )
  
}

###############################################################
# SAVE MASTER TABLE
###############################################################

write.csv(
  master.enrichment,
  file.path(summary.dir,"Enrichment_Master.csv"),
  row.names = FALSE
)

###############################################################
# DATABASE RANKING
###############################################################

database.rank <- master.enrichment %>%
  group_by(Database) %>%
  summarise(
    Tables = n(),
    Total_Terms = sum(Enriched_Terms),
    Mean_Precision = mean(Mean_Precision, na.rm = TRUE),
    Mean_Recall = mean(Mean_Recall, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Terms))

write.csv(
  database.rank,
  file.path(summary.dir,"Database_Ranking.csv"),
  row.names = FALSE
)

###############################################################
# CELL TYPE RANKING
###############################################################

cell.rank <- master.enrichment %>%
  group_by(Cell_Type) %>%
  summarise(
    Tables = n(),
    Total_Enriched_Terms = sum(Enriched_Terms),
    Mean_Precision = mean(Mean_Precision, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Enriched_Terms))

write.csv(
  cell.rank,
  file.path(summary.dir,"CellType_Enrichment_Ranking.csv"),
  row.names = FALSE
)

###############################################################
# COMPARISON RANKING
###############################################################

comparison.rank <- master.enrichment %>%
  group_by(Comparison) %>%
  summarise(
    Tables = n(),
    Total_Enriched_Terms = sum(Enriched_Terms),
    Mean_Precision = mean(Mean_Precision, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Enriched_Terms))

write.csv(
  comparison.rank,
  file.path(summary.dir,"Comparison_Enrichment_Ranking.csv"),
  row.names = FALSE
)

###############################################################
# DIRECTION SUMMARY
###############################################################

direction.rank <- master.enrichment %>%
  group_by(Direction) %>%
  summarise(
    Tables = n(),
    Total_Enriched_Terms = sum(Enriched_Terms),
    .groups = "drop"
  )

write.csv(
  direction.rank,
  file.path(summary.dir,"Direction_Ranking.csv"),
  row.names = FALSE
)

###############################################################
# REPORT
###############################################################

sink(file.path(summary.dir,"Enrichment_Report.txt"))

cat("=========================================\n")
cat("ENRICHMENT MASTER SUMMARY\n")
cat("=========================================\n\n")

cat("Tables analysed :", nrow(master.enrichment), "\n")
cat("Total enriched terms :", sum(master.enrichment$Enriched_Terms), "\n\n")

cat("Database Ranking\n")
print(database.rank)

cat("\nCell Type Ranking\n")
print(cell.rank)

cat("\nComparison Ranking\n")
print(comparison.rank)

cat("\nDirection Summary\n")
print(direction.rank)

sink()

###############################################################
# CONSOLE
###############################################################

cat("\n=========================================\n")
cat("PART 3 COMPLETE\n")
cat("=========================================\n")

cat("\nMaster table rows :", nrow(master.enrichment), "\n")
cat("Total enriched terms :", sum(master.enrichment$Enriched_Terms), "\n")