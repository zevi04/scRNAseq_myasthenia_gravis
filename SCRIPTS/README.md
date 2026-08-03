# Scripts

This directory contains the R scripts used for the **Single-Cell RNA-seq Analysis of Myasthenia Gravis (GSE227835)**.

The workflow progresses from inspection of the processed GEO files through Seurat object construction, quality control, preprocessing, dimensionality reduction, clustering, cell-type annotation, differential-expression analysis, functional enrichment, integrated DEG analysis, downstream reporting, and MG-focused directional analysis.

---

## Workflow Overview

```text
GSE227835 processed data
        ↓
Data inspection and import
        ↓
Metadata construction and validation
        ↓
Sample-level Seurat objects
        ↓
Merged Seurat object
        ↓
Quality control and filtering
        ↓
Normalization → Variable features → Scaling
        ↓
PCA → Nearest-neighbor graph → Clustering → UMAP
        ↓
Marker identification
        ↓
Cell-type annotation
        ↓
Differential expression
        ↓
Immune-cell composition
        ↓
Functional enrichment
        ↓
Biological interpretation and master summaries
        ↓
Integrated DEG analysis
        ↓
Downstream reporting / repair
        ↓
MG cell-type direction + AChR/MuSK-focused analysis
```

---

# Script Inventory

> **Note:** File extensions may be hidden in Windows Explorer. The entries below correspond to the R scripts visible in this directory.

## 00–06 · Data Inspection and Object Construction

| Script | Purpose |
| --- | --- |
| `SCRIPTS_00_inspect_data.R` | Inspects the processed GEO files and their tabular structure before import. |
| `01_import_data.R` | Imports the processed expression data into R. |
| `02_build_metadata.R` | Constructs sample-level metadata linking imported data to biological groups and sample identities. |
| `03_check_dataset_consistency.R` | Checks consistency between imported expression data, sample identities, and metadata. |
| `04_create_seurat_object.R` | Tests and establishes construction of a Seurat object from the imported data. |
| `05_create_all_seurat_objects.R` | Creates the full set of sample-level Seurat objects. |
| `06_merge_seurat_objects.R` | Merges the individual sample objects into the combined analysis object while retaining metadata. |

### Dataset state after merging

The completed import/object-construction stage produced:

- **40 individual Seurat objects**
- **444,357 cells**
- **36,601 genes**

---

## 07–09 · Quality Control

| Script | Purpose |
| --- | --- |
| `07_quality_control.R` | Calculates and evaluates cell-level quality-control metrics. |
| `08_qc_visualization.R` | Generates visual summaries of QC metrics and their distributions. |
| `09_qc_filtering.R` | Applies the final cell-filtering criteria and creates the QC-filtered dataset. |

The final QC criteria were:

```r
nFeature_RNA > 200
nFeature_RNA < 6000
percent.mt < 15
```

### QC outcome

| Metric | Value |
| --- | ---: |
| Cells before QC filtering | 444,357 |
| Cells after QC filtering | 420,538 |
| Cells removed | 23,819 |
| Percentage removed | 5.36% |

The **420,538 retained cells** formed the basis of the downstream workflow.

---

## 10–12 · Preprocessing

| Script | Purpose |
| --- | --- |
| `10_normalization.R` | Normalizes expression data for downstream analysis. |
| `11_variable_features.R` | Identifies highly variable features. |
| `12_scaling.R` | Scales the expression matrix before dimensionality reduction. |

These scripts prepare the QC-filtered expression data for PCA and graph-based analysis.

---

## 13–16 · Dimensionality Reduction and Clustering

| Script | Purpose |
| --- | --- |
| `13_PCA.R` | Performs Principal Component Analysis. |
| `14_UMAP_analysis.R` | Generates UMAP dimensionality-reduction coordinates and visualizations. |
| `15_Find_Neighbors.R` | Constructs the nearest-neighbor graph used for clustering. |
| `16_find_clusters.R` | Performs graph-based clustering of the cells. |

The project calculated **50 principal components** and ultimately identified **34 transcriptional clusters** during the unsupervised clustering stage.

> The numbering reflects the project development history. For exact computational dependencies, follow the objects loaded and saved within each script rather than assuming that filename numbering alone defines every Seurat dependency.

---

## 17–18 · Marker Identification and Cell-Type Annotation

| Script | Purpose |
| --- | --- |
| `17_find_markers.R` | Identifies cluster-associated marker genes. |
| `18_cell_type_marker_annotation.R` | Uses marker-expression evidence to assign biological cell identities to transcriptional clusters. |

The 34 computational clusters were consolidated into **10 downstream cell-type labels**:

- B cells
- Basophils
- CD4 T cells
- CD8 T cells
- Dendritic cells
- Monocytes
- Neutrophils
- NK cells
- Progenitors
- T cells

Annotation was supported by canonical-marker and top-marker visualizations stored under [`../FIGURES/`](../FIGURES/).

---

## 19–20 · Differential Expression and Cell Composition

| Script | Purpose |
| --- | --- |
| `19_differential_gene_expression_analysis.R` | Performs cell-type-specific differential-expression analysis across the biological comparisons. |
| `20_immune_cell_composition_analysis.R` | Examines the distribution of annotated immune-cell populations across the dataset and biological groups. |

The principal comparisons are:

```text
Healthy_vs_AChR_MG
Healthy_vs_SNMG_Pre
SNMG_Pre_vs_SNMG_Post
```

These represent:

1. healthy controls versus AChR-positive MG,
2. healthy controls versus pre-treatment seronegative MG,
3. seronegative MG pre-treatment versus post-treatment.

---

## 21–22 · Functional Enrichment

| Script | Purpose |
| --- | --- |
| `21_1_functional_enrichment_analysis_1.R` | Performs functional-enrichment analysis on selected DEG-derived gene sets. |
| `22_functional_enrichment_visualization.R` | Generates visual summaries of functional-enrichment results. |

Functional enrichment provides pathway- and biological-process context for transcriptional changes.

Enrichment should be interpreted as **over-representation within a gene set**, not as direct proof of pathway activation, inhibition, or causality.

---

## 23 · Biological Interpretation and Master Summary

| Script | Purpose |
| --- | --- |
| `23_biological_interpretation_and_master_summary.R` | Consolidates downstream biological information and project-level summary outputs from the preceding analyses. |

This stage acts as a bridge between the primary analyses and the more extensive integrated DEG reporting performed in Script 24.

---

# 24 · Integrated Analysis

The integrated analysis was split into multiple R files:

```text
24_1_2_integrated_analysis.R
24_3_4_integrated_analysis.R
24_5_integrated_analysis.R
24_6_7_integrated_analysis.R
```

Splitting the workflow allowed the downstream analysis to be executed and debugged in manageable sections rather than as one very large script.

---

## `24_1_2_integrated_analysis.R`

Contains the first stages of the integrated analysis.

Broadly, these stages establish the integrated downstream workspace and consolidate the differential-expression information needed by later sections.

The resulting integrated structures provide the basis for comparison of genes across:

```text
Gene
× Cell type
× Biological comparison
× Fold-change
× Adjusted significance
× Direction
```

---

## `24_3_4_integrated_analysis.R`

Contains the global DEG visualization and shared/unique DEG stages.

### Section 5.3 — Global DEG Visualization

Produces high-level summaries such as:

- DEG burden by cell type,
- DEG burden across comparisons,
- upregulated versus downregulated distributions,
- significant-DEG heatmap summaries.

### Section 5.4 — Shared and Unique DEGs

Constructs gene-occurrence structures across cell types and comparisons.

Important derived variables include:

```text
CellType_Frequency
Comparison_Frequency
```

where:

- `CellType_Frequency` = number of cell types in which a gene occurs,
- `Comparison_Frequency` = number of biological comparisons in which a gene occurs.

This stage also supports classification of genes as cell-type-specific, shared, highly shared, or conserved across represented populations.

---

## `24_5_integrated_analysis.R`

Contains the conserved disease-signature stage.

### Section 5.5 — Conserved Disease Signatures

The script summarizes gene-level statistics such as:

- mean log2 fold-change,
- mean absolute log2 fold-change,
- median log2 fold-change,
- adjusted significance,
- cell-type frequency,
- comparison frequency.

These metrics are used to prioritize recurrent transcriptional signals.

The resulting rankings represent **computationally prioritized transcriptional signatures** and should not be interpreted as experimentally validated MG biomarkers.

---

## `24_6_7_integrated_analysis.R`

Contains the later integrated-analysis and reporting stages following conserved-signature analysis.

These sections complete the Script 24 downstream workflow and generate additional summary/reporting outputs required for interpretation of the integrated analysis.

---

# 25 · Downstream Reporting Repair

| Script | Purpose |
| --- | --- |
| `25_Downstream_Reporting_Repair.R` | Repairs and regenerates selected downstream reporting outputs after problems were identified in the original integrated visualizations/tables. |

Script 25 was added after reviewing the Script 24 outputs.

The repair stage addressed downstream reporting problems such as:

- unclear or poorly scaled visualizations,
- recurrence visualization problems,
- conserved-signature reporting,
- figure readability,
- downstream table/figure generation issues.

This script should therefore be understood as a **reporting and downstream-output repair stage**, rather than a replacement for the biological preprocessing pipeline.

---

# 26 · MG Directional and AChR/MuSK Analysis

The final script visible in the repository is:

```text
26_MG_CellType_Directional_AChR_MuSK_A...
```

The filename is truncated in the Windows Explorer screenshot, but it represents the final MG-focused downstream analysis script.

Its purpose is to examine:

- significant transcriptional changes by immune-cell type,
- upregulated versus downregulated DEG burden,
- relative transcriptional perturbation across comparisons,
- cell-type-specific directional bias,
- AChR-associated target genes,
- `MUSK`,
- selected neuromuscular-junction-related genes.

The targeted panel included genes such as:

```text
CHRNA1
CHRNB1
CHRND
CHRNE
CHRNG
MUSK
LRP4
AGRN
DOK7
RAPSN
```

The final targeted analysis did **not** identify significant AChR receptor-subunit or `MUSK` transcript overexpression under the final downstream thresholds.

This distinction is important:

> **AChR-positive MG is an autoantibody-defined clinical category; it does not mean that AChR genes are transcriptionally overexpressed.**

---

# Core Pipeline vs Downstream Analysis

For repository navigation, the scripts can be divided into three broad layers.

| Layer | Scripts | Role |
| --- | --- | --- |
| **Core scRNA-seq pipeline** | 00–18 | Import, QC, preprocessing, dimensionality reduction, clustering, markers and annotation |
| **Primary biological analysis** | 19–23 | Differential expression, composition, enrichment and biological summaries |
| **Integrated / reporting analysis** | 24–26 | Cross-analysis integration, recurrence, conserved signatures, reporting repair and targeted MG interpretation |

This distinction is useful when reproducing only part of the project.

---

# Recommended Execution Order

For a complete reconstruction from the processed GEO input files, use:

```text
SCRIPTS_00_inspect_data.R
        ↓
01_import_data.R
        ↓
02_build_metadata.R
        ↓
03_check_dataset_consistency.R
        ↓
04_create_seurat_object.R
        ↓
05_create_all_seurat_objects.R
        ↓
06_merge_seurat_objects.R
        ↓
07_quality_control.R
        ↓
08_qc_visualization.R
        ↓
09_qc_filtering.R
        ↓
10_normalization.R
        ↓
11_variable_features.R
        ↓
12_scaling.R
        ↓
13_PCA.R
        ↓
14_UMAP_analysis.R
        ↓
15_Find_Neighbors.R
        ↓
16_find_clusters.R
        ↓
17_find_markers.R
        ↓
18_cell_type_marker_annotation.R
        ↓
19_differential_gene_expression_analysis.R
        ↓
20_immune_cell_composition_analysis.R
        ↓
21_1_functional_enrichment_analysis_1.R
        ↓
22_functional_enrichment_visualization.R
        ↓
23_biological_interpretation_and_master_summary.R
        ↓
24_1_2_integrated_analysis.R
        ↓
24_3_4_integrated_analysis.R
        ↓
24_5_integrated_analysis.R
        ↓
24_6_7_integrated_analysis.R
        ↓
25_Downstream_Reporting_Repair.R
        ↓
26_MG_CellType_Directional_AChR_MuSK_A...
```

---

# Running the Workflow

## 1. Clone the repository

Clone the repository locally and open the project from its root directory.

## 2. Obtain the source dataset

The original expression data are not intended to be duplicated in the GitHub repository.

Obtain the required GSE227835 processed files from the original GEO record and organize them according to the structure documented in:

[`../DOCS/Dataset_notes.md`](../DOCS/Dataset_notes.md)

The local structure is:

```text
DATA/
├── RAW/
└── PROCESSED/
```

## 3. Check paths before execution

Before running a script, confirm that:

- the project root is correct,
- expected input files exist,
- expected upstream `.rds` objects exist when resuming,
- output directories are available,
- required R packages are installed.

## 4. Run scripts sequentially

The safest full-reproduction strategy is to execute the scripts in order.

Later scripts depend on tables, metadata, serialized objects, or other outputs generated upstream.

---

# Restarting R Between Major Stages

Because this project contains more than **420,000 QC-retained cells**, memory use can become substantial.

For long downstream stages, a clean R session can be useful:

```text
Save required upstream object
        ↓
Restart R
        ↓
Clear environment
        ↓
Load only required packages
        ↓
Load required checkpoint object
        ↓
Run next stage
```

Do not clear the environment in the middle of a script that depends on objects created only in memory unless those objects have already been saved.

---

# Intermediate Objects

Large Seurat objects and `.rds` checkpoints may be excluded from GitHub.

This is intentional.

The repository is designed primarily to preserve:

- analysis code,
- documentation,
- selected tables,
- selected figures,
- reproducible workflow logic.

Large computational objects can be regenerated from upstream scripts when required.

---

# Temporary R Files

Files such as the following are not analysis scripts and should not be version controlled:

```text
.RData
.Rhistory
.RDataTmp
.RDataTmp1
.RDataTmp2
.Rproj.user/
```

These should be covered by the repository `.gitignore`.

The `.Rproj` project file may be retained in the repository because it provides a convenient reproducible RStudio project entry point.

---

# Output Locations

Scripts write outputs primarily under:

```text
RESULTS/
FIGURES/
```

Major result categories include:

```text
RESULTS/
├── Cell_Composition/
├── DEG/
├── figures/
├── Functional_Enrichment/
├── logs/
├── Master_Summary/
├── metadata/
├── objects/
├── qc/
└── Script24_Integrated_Analysis/
```

See [`../RESULTS/README.md`](../RESULTS/README.md) for result interpretation.

---

# Important Interpretation Rules

## DEG direction is not cellular activation

An increased number of upregulated DEGs indicates an upward transcriptional bias among the genes meeting the analysis thresholds. It does not independently prove physiological activation.

Similarly, a predominance of downregulated DEGs does not independently prove cellular suppression.

## AChR-positive status is not gene overexpression

The AChR-positive group is defined by disease-associated autoantibody status, not by increased expression of AChR receptor-subunit transcripts.

## Pre vs post is not automatically treatment recovery

`SNMG_Pre_vs_SNMG_Post` identifies differences between sampled pre- and post-treatment states.

A post-treatment difference should not automatically be described as normalization or recovery without demonstrating movement toward the healthy-control state.

## scRNA-seq measures transcripts

The analysis evaluates RNA-level changes. It does not directly measure:

- autoantibody concentration,
- receptor abundance,
- protein activity,
- neuromuscular-junction function,
- clinical severity.

---

# Troubleshooting Notes

Several downstream scripts were developed iteratively as issues became visible during integrated reporting.

Examples included:

- duplicated column names following joins,
- `CellType_Frequency.x` / `.y` suffixes,
- missing expected result files,
- empty or uninformative recurrence heatmaps,
- conserved-signature filtering issues,
- `slice_head()` errors caused by non-constant `n`,
- clustering failures when matrices contained `NA`, `NaN`, or `Inf`,
- visualizations requiring redesign for clearer interpretation.

These issues and the corresponding project decisions are documented in greater detail in:

[`../JOURNAL.md`](../JOURNAL.md)

The later Script 25 and Script 26 stages should therefore be understood in the context of this iterative validation and reporting process.

---

# Repository Navigation

| Location | Contents |
| --- | --- |
| [`../README.md`](../README.md) | Main repository overview and key findings |
| [`../PROJECT.md`](../PROJECT.md) | Scientific rationale, objectives and analytical scope |
| [`../JOURNAL.md`](../JOURNAL.md) | Chronological development and troubleshooting record |
| [`../DOCS/`](../DOCS/) | Dataset and technical documentation |
| [`../FIGURES/`](../FIGURES/) | Primary QC, clustering, UMAP and annotation figures |
| [`../RESULTS/`](../RESULTS/) | Analytical outputs |
| `SCRIPTS/` | Complete R workflow described in this document |

---

# Reproducibility Notes

This repository documents an independent secondary analysis of a publicly available single-cell RNA-seq dataset.

Exact reproduction can depend on:

- R version,
- Seurat version,
- package versions,
- available memory,
- source-file organization,
- presence of intermediate objects,
- random seeds where applicable.

For a fresh reconstruction, begin with the source processed files and execute the numbered scripts sequentially.

For downstream reanalysis, begin from the appropriate saved checkpoint only after confirming that its metadata and object structure match the expectations of the downstream script.

---

# Project Status

**Completed — 3 August 2026**

The script collection records the complete computational progression of the project from initial inspection of GSE227835 through the final MG-focused cell-type directional and AChR/MuSK analysis.
