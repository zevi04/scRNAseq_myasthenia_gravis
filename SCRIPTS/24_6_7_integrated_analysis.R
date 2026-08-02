
###############################################################################
# SECTION 6 : BIOLOGICAL INTERPRETATION & RESULTS WRITING
###############################################################################

section_timer("SECTION 6 : Biological Interpretation", {
  
  log_message("Generating automated biological interpretation.")
  
  overall <- read.csv(file.path(TABLE_DIR,
                                "5.2_Overall_DEG_Summary.csv"))
  celltype <- read.csv(file.path(TABLE_DIR,
                                 "5.2_CellType_DEG_Summary.csv"))
  comparison <- read.csv(file.path(TABLE_DIR,
                                   "5.2_Comparison_DEG_Summary.csv"))
  signatures <- read.csv(file.path(TABLE_DIR,
                                   "5.5_Conserved_Disease_Signatures.csv"))
  
  report <- c(
    "###############################################################################",
    "# RESULTS",
    "###############################################################################",
    "",
    "6.1 Global Differential Expression",
    "",
    sprintf("A total of %s differential expression entries representing %s unique genes were identified across %s cell types and %s disease comparisons.",
            overall$Value[overall$Metric=="Total DEG Entries"],
            overall$Value[overall$Metric=="Unique Genes"],
            overall$Value[overall$Metric=="Unique Cell Types"],
            overall$Value[overall$Metric=="Unique Comparisons"]),
    "",
    "6.2 Cell-type Specific Changes",
    "",
    paste(
      "Among the analysed cell populations,",
      celltype$CellType[1],
      "displayed the highest number of differentially expressed genes (",
      celltype$Total_DEGs[1],
      "), followed by",
      celltype$CellType[2],
      "(",
      celltype$Total_DEGs[2],
      ")."
    ),
    "",
    "6.3 Disease Comparison",
    "",
    paste(
      "The comparison exhibiting the largest transcriptional response was",
      comparison$Comparison[1],
      "with",
      comparison$Total_DEGs[1],
      "differentially expressed genes."
    ),
    "",
    "6.4 Conserved Disease Signatures",
    "",
    paste(
      nrow(signatures),
      "conserved disease signatures satisfied the predefined filtering criteria and were ranked according to the composite signature score."
    ),
    "",
    "6.5 Summary",
    "",
    "The combined analyses demonstrate both shared and cell-type-specific transcriptional alterations associated with disease progression. Conserved signatures represent candidate molecular mechanisms for downstream validation."
  )
  
  writeLines(
    report,
    file.path(REPORT_DIR,
              "Section6_Biological_Interpretation.txt")
  )
  
  if(requireNamespace("officer", quietly=TRUE)){
    doc <- officer::read_docx()
    for(line in report){
      doc <- officer::body_add_par(doc,line)
    }
    print(doc,
          target=file.path(REPORT_DIR,
                           "Section6_Biological_Interpretation.docx"))
  }
  
  save_csv(
    data.frame(
      Output=c(
        "Section6_Biological_Interpretation.txt",
        "Section6_Biological_Interpretation.docx"
      ),
      Status=c(
        "Generated",
        ifelse(requireNamespace("officer",quietly=TRUE),
               "Generated","Skipped")
      )
    ),
    "Section6_Output_Index.csv"
  )
  
  log_message("Section 6 completed successfully.")
  
})


###############################################################################
# SECTION 7 : MANUSCRIPT-READY SUPPLEMENTARY TABLES
###############################################################################

section_timer("SECTION 7 : Supplementary Tables", {
  
  log_message("Generating supplementary tables.")
  
  dir.create(file.path(REPORT_DIR, "Supplementary"), showWarnings = FALSE)
  
  if (!exists("master_deg")) {
    stop("master_deg object not found. Please run Sections 5.1–5.5 first.")
  }
  
  sig_deg <- subset(
    master_deg,
    p_val_adj < ALPHA & abs(avg_log2FC) >= LOGFC_THRESHOLD
  )
  
  up_deg <- subset(sig_deg, Direction == "Upregulated")
  down_deg <- subset(sig_deg, Direction == "Downregulated")
  
  supp_dir <- file.path(REPORT_DIR, "Supplementary")
  
  write.csv(master_deg,
            file.path(supp_dir, "Supplementary_Table_S1_All_DEGs.csv"),
            row.names = FALSE)
  
  write.csv(sig_deg,
            file.path(supp_dir, "Supplementary_Table_S2_Significant_DEGs.csv"),
            row.names = FALSE)
  
  write.csv(up_deg,
            file.path(supp_dir, "Supplementary_Table_S3_Upregulated_DEGs.csv"),
            row.names = FALSE)
  
  write.csv(down_deg,
            file.path(supp_dir, "Supplementary_Table_S4_Downregulated_DEGs.csv"),
            row.names = FALSE)
  
  if(file.exists(file.path(TABLE_DIR,
                           "5.5_Conserved_Disease_Signatures.csv"))){
    
    file.copy(
      file.path(TABLE_DIR,
                "5.5_Conserved_Disease_Signatures.csv"),
      file.path(supp_dir,
                "Supplementary_Table_S5_Conserved_Signatures.csv"),
      overwrite = TRUE
    )
    
  }
  
  manifest <- data.frame(
    Table = c(
      "S1","S2","S3","S4","S5"
    ),
    Description = c(
      "Complete DEG table",
      "Significant DEGs",
      "Upregulated DEGs",
      "Downregulated DEGs",
      "Conserved disease signatures"
    ),
    stringsAsFactors = FALSE
  )
  
  save_csv(
    manifest,
    "Section7_Supplementary_Table_Index.csv"
  )
  
  log_message("Supplementary tables generated.")
  
})


###############################################################################
# SECTION 8 : METHODOLOGY DIAGRAM
###############################################################################

section_timer("SECTION 8 : Methodology Diagram", {
  
  log_message("Generating methodology workflow diagram.")
  
  suppressPackageStartupMessages({
    library(DiagrammeR)
    library(DiagrammeRsvg)
    library(rsvg)
  })
  
  graph <- grViz("
digraph workflow {

graph [layout = dot, rankdir = TB]

node [shape = box, style = rounded, fontname = Helvetica]

A [label='Raw scRNA-seq Dataset\n(GEO: GSE227835)']
B [label='Quality Control']
C [label='Normalization']
D [label='Dimensionality Reduction']
E [label='Cell Clustering']
F [label='Cell Annotation']
G [label='Differential Expression']
H [label='Integrated DEG Analysis']
I [label='Shared & Conserved\nDisease Signatures']
J [label='Biomarker Discovery']
K [label='Publication Figures']
L [label='Biological Interpretation']

A->B->C->D->E->F->G->H->I->J->K->L

}
")
  
  svg_txt <- export_svg(graph)
  
  svg_file <- file.path(FIGURE_DIR,
                        "Section8_Methodology_Workflow.svg")
  writeLines(svg_txt, svg_file)
  
  rsvg_pdf(svg_file,
           file.path(FIGURE_DIR,
                     "Section8_Methodology_Workflow.pdf"))
  
  rsvg_png(svg_file,
           file.path(FIGURE_DIR,
                     "Section8_Methodology_Workflow.png"),
           width = 2400,
           height = 3200)
  
  report <- data.frame(
    Output=c(
      "Section8_Methodology_Workflow.svg",
      "Section8_Methodology_Workflow.pdf",
      "Section8_Methodology_Workflow.png"
    ),
    Status="Generated",
    stringsAsFactors=FALSE
  )
  
  save_csv(
    report,
    "Section8_Methodology_Output_Index.csv"
  )
  
  log_message("Section 8 completed successfully.")
  
})
