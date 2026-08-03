<div align="center">
🧬 Single-Cell RNA-seq Analysis of Myasthenia Gravis
Cell-type-specific transcriptional profiling of GSE227835
![R](https://img.shields.io/badge/R-4.x-276DC3?style=flat-square&logo=r&logoColor=white)
![Seurat](https://img.shields.io/badge/Seurat-v5-6A5ACD?style=flat-square)
![scRNA-seq](https://img.shields.io/badge/Analysis-scRNA--seq-009688?style=flat-square)
![MG](https://img.shields.io/badge/Disease-Myasthenia%20Gravis-B23A48?style=flat-square)
![Status](https://img.shields.io/badge/Status-Completed-2E8B57?style=flat-square)
Independent reanalysis of publicly available single-cell RNA-sequencing data to characterize immune-cell-specific transcriptional alterations in Myasthenia Gravis.
</div>
---
📌 Overview
Myasthenia Gravis (MG) is an autoimmune disorder of the neuromuscular junction. This project reanalyzes GSE227835 to examine how MG-associated transcriptional alterations differ across immune-cell populations.
The workflow progresses from data inspection and quality control through clustering, immune-cell annotation, cell-type-specific differential expression, DEG recurrence, conserved-signature analysis, functional-enrichment reporting, and targeted AChR/MuSK investigation.
> **Central question:** Are MG-associated transcriptional alterations broadly shared across immune cells, or predominantly cell-type specific?
---
🎯 Biological Questions
Which immune-cell populations are represented?
How does AChR-positive MG differ transcriptionally from healthy controls?
How does pre-treatment seronegative MG differ from healthy controls?
What changes occur between seronegative MG pre- and post-treatment states?
Which cell types show the greatest transcriptional perturbation?
Which populations show predominantly upward or downward DEG patterns?
Which DEGs are cell-type specific, shared, recurrent, or conserved?
Can recurrent transcriptional signatures be prioritized?
Is there significant transcript-level evidence involving AChR receptor-subunit genes or `MUSK`?
---
🗂️ Dataset and Study Design
NCBI GEO accession: `GSE227835`
The analysis contains four biological groups:
Healthy Controls
AChR-positive MG
Seronegative MG — Pre-treatment
Seronegative MG — Post-treatment
Dataset scale
Metric	Value
Samples	40
Cells before QC	444,357
Genes	36,601
Cells after QC	420,538
Cells removed	23,819 (5.36%)
Initial transcriptional clusters	34
Downstream cell-type labels	10
Main comparisons
Comparison	Purpose
`Healthy_vs_AChR_MG`	AChR-positive MG versus healthy controls
`Healthy_vs_SNMG_Pre`	Pre-treatment seronegative MG versus healthy controls
`SNMG_Pre_vs_SNMG_Post`	Seronegative MG pre- versus post-treatment
Comparisons were evaluated within individual immune-cell populations.
---
🔬 Computational Workflow
```text
GSE227835
   ↓
Data inspection + metadata
   ↓
40 individual Seurat objects
   ↓
Merged dataset — 444,357 cells
   ↓
Quality control
   ↓
420,538 retained cells
   ↓
LogNormalize → Variable Features → Scaling
   ↓
PCA → Neighbors → Leiden Clustering → UMAP
   ↓
34 transcriptional clusters
   ↓
Cell-type annotation
   ↓
10 downstream cell populations
   ↓
Cell-type-specific differential expression
   ↓
Global DEG summaries
   ├── Shared / unique DEGs
   ├── Gene recurrence
   ├── Conserved signatures
   ├── Functional-enrichment reporting
   └── Directional transcriptional analysis
   ↓
AChR / MuSK targeted analysis
   ↓
Biological interpretation
```
---
🧹 Quality Control and Preprocessing
QC used `nFeature_RNA`, `nCount_RNA`, and `percent.mt`.
```r
nFeature_RNA > 200
nFeature_RNA < 6000
percent.mt < 15
```
After filtering, 420,538 cells were retained. Data were normalized using Seurat `LogNormalize` with a scale factor of 10,000, and 2,000 highly variable genes were selected using VST.
PCA used 50 principal components. The reduced representation was used for nearest-neighbor graph construction, Leiden clustering, and UMAP. The initial unsupervised analysis identified 34 clusters.
---
🧫 Cell-Type Annotation
The downstream analysis used 10 labels:
Cell types	Cell types	Cell types
B cells	Basophils	CD4 T cells
CD8 T cells	Dendritic cells	Monocytes
Neutrophils	NK cells	Progenitors
T cells		
Annotation allowed differential-expression testing within biologically comparable populations.
---
📊 Differential Expression
The integrated DEG framework retained gene, cell type, comparison, average log2 fold-change, adjusted P value, and direction.
Final reporting thresholds:
```text
Adjusted P < 0.05
|average log2 fold-change| ≥ 0.5
```
> **Important:** DEG direction is not a direct assay of physiological cellular activation. The project therefore uses **transcriptional perturbation** for broader cell-type comparisons.
---
🔎 Key Results
AChR-positive MG
In `Healthy_vs_AChR_MG`, neutrophils showed a prominent upward profile:
Cell type	Upregulated	Downregulated
Neutrophils	27	12
B cells	9	50
Dendritic cells	2	37
CD4 T cells	1	27
NK cells	3	24
Neutrophils showed a 69.2% upward directional bias among significant DEG entries. Several other populations were predominantly shifted downward.
---
Seronegative MG — Pre-treatment
In `Healthy_vs_SNMG_Pre`:
Cell type	Upregulated	Downregulated
Neutrophils	39	1
B cells	30	107
Dendritic cells	20	80
NK cells	9	40
Monocytes	59	57
Neutrophils showed a 97.5% upward directional bias. Monocytes showed substantial DEG burden but an almost balanced up/down distribution, consistent with extensive transcriptional remodeling rather than a simple directional response.
---
Seronegative MG — Pre vs Post
Cell type	Upregulated	Downregulated	Pattern
B cells	7	0	Upward
NK cells	17	5	Upward
Dendritic cells	13	6	Upward
Neutrophils	8	52	Downward
CD8 T cells	1	20	Downward
T cells	3	20	Downward
Neutrophils showed an 86.7% downward directional bias.
These pre/post differences are not automatically interpreted as treatment-induced normalization. Demonstrating normalization would require showing that the same disease-associated genes move toward healthy-control expression.
---
🔁 Shared, Recurrent and Conserved DEGs
A gene × cell-type occurrence framework calculated:
CellType_Frequency — number of cell types in which a gene occurred.
Comparison_Frequency — number of comparisons in which a gene occurred.
Genes were classified as cell-type specific, shared across two cell types, highly shared, or conserved across represented cell types.
The recurrence analysis indicated that many significant DEGs were restricted to relatively few cell types, while progressively fewer genes recurred across many populations. This supports a substantial cell-type-dependent component of the MG-associated transcriptional response.
Recurrent genes were additionally prioritized using effect size, adjusted significance, cell-type recurrence, and comparison recurrence. This conserved-signature analysis is a prioritization framework, not proof that a prioritized gene is a causal MG biomarker.
---
🎯 AChR / MuSK Targeted Analysis
The targeted panel included genes such as `CHRNA1`, `CHRNB1`, `CHRND`, `CHRNE`, `CHRNG`, `MUSK`, `LRP4`, `AGRN`, `DOK7`, and `RAPSN`.
Among examined AChR receptor-subunit genes, only `CHRNB1` appeared in the master DEG table:
Feature	Observation
Cell type	Progenitors
Comparison	Healthy vs AChR-MG
avg_log2FC	≈ −0.535
Adjusted P	1.0
Significant by final thresholds	No
`MUSK` was not represented as a significant DEG.
Target	Significant DEG records
AChR receptor-subunit genes	0
`MUSK`	0
> **AChR-positive MG refers to autoantibody status, not AChR-gene overexpression.** The same distinction applies to MuSK-associated antibody status.
---
🧠 Overall Interpretation
The principal observation is that MG-associated transcriptional alterations are strongly dependent on immune-cell identity.
Neutrophils showed prominent upward signatures in both disease-versus-healthy comparisons.
B cells showed substantial downward transcriptional changes in disease-versus-healthy comparisons.
Dendritic cells, NK cells, CD4 T cells, and CD8 T cells showed comparison-dependent patterns.
Monocytes could show substantial DEG burden without a strong net directional bias.
The seronegative pre/post comparison produced a different landscape, including pronounced downward changes in neutrophils.
DEG recurrence suggested that many alterations were restricted to relatively few cell types.
Significant AChR receptor-subunit or `MUSK` transcript overexpression was not supported.
---
🖼️ Figures
The repository contains:
QC distributions
PCA and UMAP
clustering visualizations
canonical-marker and top-marker dot plots
marker heatmaps
DEG summaries
cell-composition outputs
recurrence and conserved-signature analyses
directional transcriptional analyses
targeted AChR/MuSK figures
See `FIGURES/README.md` and `RESULTS/README.md` for organization details.
---
📁 Repository Structure
```text
scRNAseq_myasthenia_gravis/
│
├── DATA/
│   ├── PROCESSED/
│   └── RAW/
│
├── DOCS/
│   ├── Dataset_notes.md
│   └── README.md
│
├── FIGURES/
│   ├── CLUSTERING/
│   ├── PCA/
│   ├── qc/
│   ├── UMAP/
│   ├── Dotplot_CanonicalMarkers.png
│   ├── Dotplot_TopMarkers.png
│   ├── Heatmap_TopMarkers.png
│   └── README.md
│
├── RESULTS/
│   ├── Cell_Composition/
│   ├── DEG/
│   ├── figures/
│   ├── Functional_Enrichment/
│   ├── logs/
│   ├── Master_Summary/
│   ├── metadata/
│   ├── objects/
│   ├── qc/
│   ├── Script24_Integrated_Analysis/
│   └── README.md
│
├── SCRIPTS/
│   ├── Scripts 00–26
│   └── README.md
│
├── .gitattributes
├── .gitignore
├── JOURNAL.md
├── PROJECT.md
└── README.md
```
---
⚙️ Reproducibility
The workflow was developed in R / RStudio, primarily using Seurat.
Major packages include:
```r
Seurat
dplyr
tidyr
ggplot2
data.table
ComplexHeatmap
circlize
reshape2
patchwork
```
Large expression files and serialized Seurat objects are intentionally excluded from GitHub. Obtain GSE227835 data and place them under the local `DATA/` structure.
For the exact execution order, see `SCRIPTS/README.md`.
Analysis stages
Stage	Scripts	Purpose
Data preparation	`00–06`	Inspection, metadata, Seurat-object creation and merging
Quality control	`07–09`	QC calculation, visualization and filtering
Preprocessing	`10–12`	Normalization, variable features and scaling
Reduction/clustering	`13–16`	PCA, UMAP, neighbors and clustering
Annotation	`17–18`	Marker analysis and cell-type annotation
Primary analysis	`19–23`	DEGs, composition, enrichment and interpretation
Integrated analysis	`24.x`	Global integration, recurrence, signatures and reporting
Reporting repair	`25`	Refined downstream figures, legends and reports
MG directional analysis	`26`	Cell-type direction and AChR/MuSK target analysis
---
⚠️ Limitations
The >400,000-cell dataset imposed substantial RAM and storage requirements.
DEG burden and direction are not direct measurements of cellular activation or suppression.
Pre/post transcriptional differences alone do not establish treatment-induced normalization.
Autoantibody-defined disease categories should not be equated with target-gene overexpression.
Differential expression identifies association, not causation.
Functional conclusions require pathway-level, protein-level, experimental, or independent-cohort validation.
---
🔭 Future Directions
Potential extensions include pathway enrichment of directional cell-type DEGs, gene-set enrichment, pathway-activity scoring, deeper neutrophil/B-cell/T-cell investigation, explicit healthy-state normalization testing, independent MG-dataset validation, and multi-omic integration.
These are future extensions, not unfinished requirements of the current project.
---
<div align="center">
✅ Project Status
COMPLETED — 3 August 2026
Cell-type-specific differential-expression characterization and downstream descriptive analysis of Myasthenia Gravis using GSE227835
</div>
Completed analytical stages include data inspection, metadata construction, Seurat object creation and merging, QC, preprocessing, PCA/UMAP, clustering, annotation, differential expression, cell composition, functional-enrichment analysis, integrated DEG analysis, recurrence analysis, conserved signatures, downstream reporting, directional transcriptional analysis, targeted AChR/MuSK analysis, and biological interpretation.
---
💾 Data Availability
This project uses publicly available single-cell RNA-sequencing data from NCBI GEO accession GSE227835.
The original expression data are not redistributed through this repository. Large raw/processed datasets and serialized Seurat objects are excluded through `.gitignore`.
---
📚 Documentation
Document	Purpose
`PROJECT.md`	Scientific scope, objectives and analytical design
`JOURNAL.md`	Chronological research record and troubleshooting
`DOCS/Dataset_notes.md`	Dataset-specific notes
`DOCS/README.md`	Documentation guide
`FIGURES/README.md`	Figure organization
`RESULTS/README.md`	Generated-output guide
`SCRIPTS/README.md`	Script order, purpose and execution notes
---
👤 Author
Zeno Vimalan A.
M.Sc. Zoology  
Life Sciences · Single-Cell Transcriptomics · Bioinformatics
GitHub: @zevi04
Independent bioinformatics project focused on practical single-cell transcriptomics, reproducible computational analysis, and disease-oriented biological interpretation.
---
🙏 Acknowledgements
The original investigators who generated and publicly deposited GSE227835 are acknowledged for making the dataset available for secondary analysis.
The project uses open-source tools from the R and single-cell genomics communities, particularly Seurat.
> **Disclaimer:** This repository contains an independent computational reanalysis of publicly available data. The biological interpretations are exploratory and should not be considered clinical conclusions.
---
<div align="center">
GSE227835 · Myasthenia Gravis · scRNA-seq · Seurat · R
</div>