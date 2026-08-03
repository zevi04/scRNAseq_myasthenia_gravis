# Project Charter

## Project Title

**Single-Cell Transcriptomic Analysis of Immune Cell Populations in Myasthenia Gravis**

**Dataset:** GSE227835  
**Analysis type:** Single-cell RNA sequencing (scRNA-seq) reanalysis  
**Primary framework:** R / Seurat  
**Project status:** Completed  
**Completion date:** 3 August 2026

---

## Project Overview

This project presents a single-cell RNA-sequencing reanalysis of the publicly available **GSE227835** dataset to investigate immune-cell-specific transcriptional alterations associated with Myasthenia Gravis (MG). The analysis was designed to move beyond whole-dataset differential expression by examining disease-associated transcriptional changes within individual immune-cell populations.

The workflow included dataset inspection, metadata construction, Seurat object generation, quality control, normalization, dimensionality reduction, unsupervised clustering, cell-type annotation, differential gene-expression analysis, immune-cell composition analysis, functional enrichment, DEG recurrence analysis, conserved-signature identification, and targeted evaluation of AChR/MuSK-associated genes.

The study evaluates three principal biological comparisons: **healthy controls versus AChR-positive MG**, **healthy controls versus pre-treatment seronegative MG**, and **seronegative MG pre-treatment versus post-treatment**. This design enables both disease-versus-healthy characterization and assessment of transcriptional differences between the two seronegative MG time points.

The final analysis focuses not only on the number of differentially expressed genes, but also on their **direction, cell-type distribution, recurrence across immune populations, and consistency across comparisons**. Particular attention is given to identifying immune-cell populations showing substantial transcriptional perturbation and distinguishing broadly recurrent signals from cell-type-specific responses.

---

## Scientific Rationale

Myasthenia Gravis is an autoimmune disorder in which disruption of neuromuscular transmission is associated with pathogenic autoantibodies, most commonly against the acetylcholine receptor (AChR), while other disease subgroups include seronegative MG and antibody-defined forms involving other neuromuscular-junction targets.

Although MG is clinically defined by impaired neuromuscular transmission, the autoimmune response involves multiple immune-cell populations. An analysis performed only at the aggregate transcriptomic level can obscure this heterogeneity because transcriptional changes occurring within one population may be diluted or combined with unrelated changes occurring in another.

Single-cell RNA sequencing provides a framework for resolving this heterogeneity by examining transcriptional states at the level of individual cells and subsequently comparing defined immune-cell populations between biological conditions. This makes it possible to ask not simply whether gene expression differs in MG, but **which immune-cell populations carry those differences and whether the observed transcriptional programs are shared or cell-type restricted**.

The inclusion of both **AChR-positive MG and seronegative MG** allows the analysis to examine transcriptional patterns across clinically distinct MG groups rather than assuming a single molecular profile for the disease. The seronegative pre- versus post-treatment comparison provides an additional view of how the immune transcriptome differs between the two sampled time points.

A further objective is to distinguish **autoantibody status from transcript-level expression of the corresponding antigen genes**. AChR-positive disease does not imply transcriptional overexpression of AChR receptor-subunit genes, and the same principle applies to MuSK-associated disease. Targeted analysis of AChR/MuSK and related neuromuscular-junction genes was therefore incorporated to test whether such transcript-level evidence was present in the DEG dataset rather than inferring it from disease classification.

Finally, recurrence and conserved-signature analyses were incorporated to distinguish genes detected in isolated cell populations from genes repeatedly observed across multiple immune-cell types or biological comparisons. This provides a structured approach for prioritizing potentially broader MG-associated transcriptional signatures while retaining the cell-type-specific context in which those signals occur.

---

## Primary Research Question

**How do cell-type-specific transcriptional profiles differ between healthy controls, AChR-positive Myasthenia Gravis, and seronegative Myasthenia Gravis, and which differentially expressed genes and transcriptional patterns are shared, recurrent, or condition-specific across immune-cell populations?**

### Secondary Research Questions

1. Which immune-cell populations are represented in the analyzed dataset?
2. Which cell types show the greatest transcriptional perturbation in MG-associated comparisons?
3. Which cell types show predominantly upward or downward differential-expression patterns?
4. How does pre-treatment seronegative MG differ from healthy controls?
5. How does the immune-cell transcriptional landscape differ between seronegative MG pre- and post-treatment states?
6. Which differentially expressed genes are cell-type specific, shared, recurrent, or conserved?
7. Which recurrent genes can be prioritized as broader transcriptional signatures?
8. Do AChR receptor-subunit genes, `MUSK`, or related neuromuscular-junction genes show significant differential expression in the analyzed immune-cell DEG dataset?

---

## Project Objectives

### Objective 1 — Construct a reproducible single-cell analysis workflow

Develop a structured R/Seurat workflow beginning with processed public data and progressing through metadata construction, Seurat object creation, merging, quality control, preprocessing, dimensionality reduction, clustering, annotation, and downstream analysis.

### Objective 2 — Characterize the cellular landscape

Identify transcriptionally distinct clusters and assign biologically interpretable immune-cell identities using marker-based annotation.

### Objective 3 — Perform cell-type-specific differential-expression analysis

Evaluate gene-expression differences within individual cell populations rather than relying solely on aggregate comparisons across all cells.

### Objective 4 — Compare major MG-associated biological states

Evaluate:

- healthy controls versus AChR-positive MG,
- healthy controls versus pre-treatment seronegative MG,
- seronegative MG pre-treatment versus post-treatment.

### Objective 5 — Characterize transcriptional direction and perturbation

Quantify significant upward and downward DEG patterns across cell types and comparisons to determine which populations exhibit the greatest transcriptional remodeling.

### Objective 6 — Identify shared and recurrent transcriptional signals

Construct gene-occurrence summaries across cell types and comparisons to distinguish cell-type-specific genes from more widely recurrent genes.

### Objective 7 — Prioritize conserved disease-associated signatures

Integrate recurrence, effect size, and statistical information to prioritize genes showing broader transcriptional consistency.

### Objective 8 — Examine AChR/MuSK-related targets

Investigate whether genes encoding AChR receptor subunits, `MUSK`, and selected neuromuscular-junction-associated components are represented among significant immune-cell DEGs.

### Objective 9 — Produce reproducible outputs and documentation

Maintain scripts, figures, tables, logs, project notes, and a chronological research journal so that the analytical process and its limitations remain transparent.

---

## Dataset

### Source

**NCBI Gene Expression Omnibus (GEO)**  
**Accession:** `GSE227835`

The project uses publicly available processed single-cell RNA-sequencing data. Large source files and serialized Seurat objects are not redistributed through the GitHub repository.

### Biological Groups

The analysis includes four study groups:

1. **Healthy Controls**
2. **AChR-positive Myasthenia Gravis**
3. **Seronegative MG — Pre-treatment**
4. **Seronegative MG — Post-treatment**

### Dataset Scale

| Metric | Value |
| --- | ---: |
| Individual sample objects | 40 |
| Cells before QC | 444,357 |
| Genes | 36,601 |
| Cells after QC | 420,538 |
| Cells removed during QC | 23,819 |
| Percentage removed | 5.36% |
| Initial transcriptional clusters | 34 |
| Downstream cell-type labels | 10 |

---

## Study Design

### Primary Comparisons

| Comparison | Analytical purpose |
| --- | --- |
| `Healthy_vs_AChR_MG` | Identify transcriptional differences associated with AChR-positive MG relative to healthy controls |
| `Healthy_vs_SNMG_Pre` | Identify transcriptional differences associated with pre-treatment seronegative MG relative to healthy controls |
| `SNMG_Pre_vs_SNMG_Post` | Characterize transcriptional differences between seronegative MG pre- and post-treatment states |

The analysis was structured around **within-cell-type comparisons** so that transcriptional changes could be interpreted in the context of specific immune populations.

---

## Quality-Control Strategy

Quality control was based on standard Seurat metrics:

- `nFeature_RNA`
- `nCount_RNA`
- `percent.mt`

Cells were retained using the following thresholds:

```r
nFeature_RNA > 200
nFeature_RNA < 6000
percent.mt < 15
```

The merged dataset contained **444,357 cells** before filtering. Following QC, **420,538 cells** were retained for downstream analysis.

---

## Preprocessing and Dimensionality Reduction

Following quality control:

1. Expression data were normalized using Seurat's `LogNormalize` procedure.
2. A scale factor of 10,000 was used during normalization.
3. **2,000 highly variable genes** were selected using the VST method.
4. Data were scaled for downstream dimensionality reduction.
5. Principal Component Analysis was performed using **50 principal components**.
6. A nearest-neighbor graph was constructed.
7. Leiden clustering was performed.
8. UMAP was used to visualize the reduced transcriptional space.

The unsupervised analysis initially identified **34 transcriptional clusters**.

---

## Cell-Type Annotation

The 34 computational clusters were interpreted and consolidated into **10 downstream cell-type labels**:

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

Marker-based annotation was supported by canonical-marker and top-marker visualization, including dot plots and heatmaps.

---

## Differential-Expression Framework

Differential-expression results were organized into an integrated framework containing information such as:

- gene,
- cell type,
- comparison,
- average log2 fold-change,
- adjusted P value,
- transcriptional direction.

For final downstream reporting, significant DEG records were defined using:

```text
Adjusted P value < 0.05
|average log2 fold-change| ≥ 0.5
```

Direction was assigned according to the sign of the fold-change:

```text
avg_log2FC > 0  → Upregulated
avg_log2FC < 0  → Downregulated
```

These classifications describe **transcriptional direction** and are not treated as direct measurements of physiological cellular activation or inhibition.

---

## Analytical Modules

### 1. Data Inspection and Object Construction

The initial scripts inspect the downloaded files, establish dataset structure, construct metadata, verify consistency, create individual Seurat objects, and merge the samples into a unified dataset.

### 2. Quality Control

QC scripts calculate and visualize cell-level metrics and apply filtering thresholds before downstream analysis.

### 3. Preprocessing

Normalization, variable-feature identification, and scaling prepare the expression matrix for dimensionality reduction.

### 4. Dimensionality Reduction and Clustering

PCA, nearest-neighbor graph construction, Leiden clustering, and UMAP are used to identify and visualize transcriptional structure.

### 5. Marker Identification and Cell-Type Annotation

Cluster markers and canonical immune markers are used to assign downstream cell-type identities.

### 6. Differential Gene Expression

Cell-type-specific differential expression is evaluated across the three principal biological comparisons.

### 7. Immune-Cell Composition

Cell-population distributions are summarized to provide complementary information about the cellular composition of the dataset.

### 8. Functional Enrichment

Functional-enrichment analyses provide a pathway-oriented framework for interpreting selected transcriptional changes.

### 9. Integrated DEG Analysis

Downstream DEG results are consolidated to generate global summaries, directional comparisons, shared/unique gene analyses, recurrence matrices, and reporting outputs.

### 10. Gene Recurrence Analysis

Genes are evaluated according to their occurrence across cell types and biological comparisons.

Two principal recurrence measures are used:

- **CellType_Frequency** — number of cell types in which a gene occurs.
- **Comparison_Frequency** — number of biological comparisons in which a gene occurs.

This allows genes to be separated into more cell-type-restricted and more broadly recurrent categories.

### 11. Conserved-Signature Prioritization

Recurrent genes are further evaluated using information including effect size, statistical significance, and recurrence.

This module is intended for **prioritization**, not causal biomarker identification.

### 12. Directional Cell-Type Analysis

Upregulated and downregulated DEG counts are compared across cell types and biological conditions to characterize the direction and extent of transcriptional perturbation.

### 13. AChR/MuSK Targeted Analysis

A targeted panel includes AChR receptor-subunit and neuromuscular-junction-associated genes such as:

`CHRNA1`, `CHRNB1`, `CHRND`, `CHRNE`, `CHRNG`, `MUSK`, `LRP4`, `AGRN`, `DOK7`, and `RAPSN`.

This analysis tests whether these genes are represented among significant immune-cell DEGs rather than assuming transcript-level overexpression from antibody-defined disease status.

---

## Interpretation Framework

The project follows several interpretive rules intended to prevent overstatement of computational findings.

### DEG burden is not equivalent to cellular activation

A cell type with many upregulated DEGs is not automatically described as physiologically "activated," and a cell type with many downregulated DEGs is not automatically described as "suppressed."

The preferred description is **transcriptional perturbation**, **upward transcriptional bias**, or **downward transcriptional bias**, depending on the result.

### Pre/post change is not automatically treatment response

Differences between seronegative MG pre- and post-treatment samples demonstrate transcriptional differences between the sampled states.

They are not automatically interpreted as treatment-induced normalization. Demonstrating normalization would require showing that disease-associated expression changes move toward the healthy-control state.

### Autoantibody status is distinct from gene expression

AChR-positive MG describes antibody status and does not imply that AChR receptor-subunit genes are transcriptionally overexpressed.

Similarly, MuSK-associated autoimmunity should not be equated with `MUSK` overexpression.

### Recurrence does not establish causality

A recurrent gene may represent a potentially broader transcriptional signature, but recurrence alone does not establish that the gene is causal, disease-specific, or clinically useful as a biomarker.

### Computational association requires validation

Differential-expression and enrichment analyses generate biological hypotheses. Functional conclusions require independent transcriptomic cohorts, protein-level evidence, functional assays, or other experimental validation.

---

## Major Observational Outcomes

The completed downstream analysis identified several broad patterns.

### AChR-positive MG

Neutrophils showed a prominent upward DEG pattern in the AChR-positive MG versus healthy comparison, while B cells, dendritic cells, CD4 T cells, and NK cells showed predominantly downward DEG patterns.

### Pre-treatment Seronegative MG

Neutrophils again showed a strong upward directional pattern. B cells, dendritic cells, and NK cells showed substantial downward changes, while monocytes displayed a high but approximately balanced DEG burden.

### Seronegative MG Pre vs Post

The pre/post comparison showed a different transcriptional configuration. Neutrophils displayed a pronounced downward pattern, while B cells, NK cells, and dendritic cells showed predominantly upward changes.

### Gene Recurrence

Many significant DEGs were restricted to relatively few cell types, while a smaller set recurred more broadly. This supports a substantial cell-type-specific component to the observed MG-associated transcriptional landscape.

### AChR/MuSK Targets

The targeted analysis did not identify significant AChR receptor-subunit or `MUSK` DEG records under the final reporting thresholds. This reinforces the distinction between antibody-defined MG categories and transcript-level expression of the corresponding antigen genes.

Detailed numerical results and biological interpretation are maintained in the root `README.md` and downstream result outputs rather than duplicated extensively in this charter.

---

## Expected and Generated Outputs

The completed project generates or maintains:

### Data and Metadata

- processed local data structure,
- sample metadata,
- cell-level metadata,
- intermediate Seurat objects where required locally.

### Quality-Control Outputs

- QC summaries,
- QC distributions,
- filtering reports.

### Dimensionality-Reduction Outputs

- PCA visualizations,
- UMAP visualizations,
- clustering figures.

### Annotation Outputs

- cluster-marker tables,
- canonical-marker dot plots,
- top-marker dot plots,
- marker heatmaps,
- annotated cell identities.

### Differential-Expression Outputs

- cell-type-specific DEG tables,
- comparison summaries,
- upregulated/downregulated summaries,
- integrated master DEG outputs.

### Downstream Outputs

- cell-composition summaries,
- functional-enrichment outputs,
- shared and unique DEG tables,
- gene-occurrence matrices,
- recurrence summaries,
- conserved-signature tables,
- directional cell-type analyses,
- AChR/MuSK targeted-analysis outputs.

### Documentation

- root project README,
- project charter,
- chronological research journal,
- dataset notes,
- script documentation,
- figure documentation,
- results documentation,
- analysis logs.

---

## Repository Organization

The repository is organized into the following major components:

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
├── JOURNAL.md
├── PROJECT.md
└── README.md
```

---

## Reproducibility Strategy

The project was structured as a numbered R-script workflow rather than a single monolithic analysis file.

Broad script stages are:

| Stage | Scripts | Purpose |
| --- | --- | --- |
| Data preparation | `00–06` | Data inspection, metadata, Seurat object construction and merging |
| Quality control | `07–09` | QC calculation, visualization and filtering |
| Preprocessing | `10–12` | Normalization, variable-feature selection and scaling |
| Reduction and clustering | `13–16` | PCA, UMAP, neighbors and clustering |
| Annotation | `17–18` | Marker identification and cell-type annotation |
| Primary analysis | `19–23` | Differential expression, composition, enrichment and interpretation |
| Integrated analysis | `24.x` | Global DEG analysis, recurrence, signatures and reporting |
| Reporting refinement | `25` | Downstream figure and reporting refinement |
| MG directional analysis | `26` | Cell-type directional comparison and AChR/MuSK analysis |

Large raw/processed data and serialized R objects are excluded from GitHub where appropriate. The repository therefore contains the **analysis logic and documentation**, while large reproducibility inputs must be obtained separately.

---

## Project Scope

### In Scope

- public scRNA-seq data reanalysis,
- Seurat-based preprocessing,
- QC and filtering,
- dimensionality reduction,
- clustering,
- immune-cell annotation,
- cell-type-specific differential expression,
- cell-composition analysis,
- functional-enrichment analysis,
- DEG recurrence analysis,
- conserved-signature prioritization,
- directional transcriptional comparison,
- AChR/MuSK targeted analysis,
- reproducible documentation and reporting.

### Outside the Current Scope

The following are not required for completion of the present project:

- experimental validation,
- protein-level validation,
- clinical diagnostic model development,
- causal inference,
- prospective patient recruitment,
- treatment-efficacy claims,
- definitive biomarker validation,
- multi-omic integration,
- independent-cohort replication.

These may be considered future extensions.

---

## Limitations

### Computational Scale

The dataset contains more than 400,000 QC-filtered cells. Several Seurat operations therefore impose substantial RAM, storage, and processing requirements on desktop hardware.

### Dependence on Available Metadata

Interpretation is limited to the biological and clinical information represented in the public dataset and project metadata.

### Cell-Type Annotation

Cell-type labels are computational interpretations based on transcriptional markers. Closely related immune populations may require more detailed subclustering or orthogonal validation for finer resolution.

### Differential Expression

DEG results are dependent on the analytical design, statistical method, thresholds, cell representation, and available samples.

### Functional Enrichment

Enrichment results indicate over-representation of annotated biological functions or pathways; they do not demonstrate pathway activity experimentally.

### Pre/Post Interpretation

Differences between pre- and post-treatment states cannot by themselves establish therapeutic response, causality, or normalization.

### Targeted AChR/MuSK Analysis

Absence of significant transcript-level differential expression does not imply absence of AChR- or MuSK-related disease mechanisms, because MG pathogenesis involves autoantibodies and protein-level interactions that are not directly measured by this transcriptomic analysis.

---

## Future Directions

Potential future extensions include:

1. Separate pathway analysis of cell-type-specific upregulated and downregulated genes.
2. Gene-set enrichment analysis.
3. Pathway-activity scoring at the single-cell or cell-type level.
4. Deeper investigation of neutrophil-associated transcriptional programs.
5. Focused B-cell, T-cell, NK-cell, dendritic-cell, and monocyte analyses.
6. Explicit testing of whether disease-associated genes move toward healthy expression after the post-treatment time point.
7. Validation using an independent MG transcriptomic dataset.
8. Integration with proteomic, repertoire, or other omics data where compatible datasets are available.
9. More detailed subclustering of selected immune populations.
10. Experimental validation of prioritized transcriptional signatures.

These are extensions beyond the completion criteria of the current project.

---

## Relationship to Repository Documentation

This file serves as the **scientific and analytical charter** for the project.

For other levels of documentation:

- [`README.md`](README.md) — public-facing project overview, major results, repository navigation, and reproducibility summary.
- [`JOURNAL.md`](JOURNAL.md) — chronological record of project development, troubleshooting, decisions, and observations.
- [`DOCS/Dataset_notes.md`](DOCS/Dataset_notes.md) — dataset-specific notes.
- [`SCRIPTS/README.md`](SCRIPTS/README.md) — script order and execution guidance.
- [`FIGURES/README.md`](FIGURES/README.md) — figure organization.
- [`RESULTS/README.md`](RESULTS/README.md) — downstream-output organization.

---

## Completion Criteria

The project was considered complete when the workflow had progressed through:

- data inspection and import,
- metadata construction,
- Seurat object creation and merging,
- quality control,
- normalization and preprocessing,
- PCA and UMAP,
- clustering,
- marker identification,
- cell-type annotation,
- cell-type-specific differential expression,
- cell-composition analysis,
- functional-enrichment analysis,
- integrated downstream DEG analysis,
- shared/unique and recurrent DEG analysis,
- conserved-signature prioritization,
- downstream reporting refinement,
- directional MG cell-type analysis,
- targeted AChR/MuSK analysis,
- final documentation and interpretation.

---

## Project Status

**Status: COMPLETED**  
**Completion date: 3 August 2026**

The project is complete at the level of:

> **Cell-type-specific differential-expression characterization and downstream descriptive analysis of Myasthenia Gravis using GSE227835.**

Further pathway-level, experimental, or cross-dataset analyses are considered independent future extensions rather than unfinished components of the current project.

---

## Original Charter

The project began with the following concise charter:

> **Title:** Single-cell transcriptomic analysis of immune cell populations in Myasthenia Gravis  
> **Dataset:** GSE227835  
> **Original research question:** How does the immune cell landscape differ between patients with Myasthenia Gravis and healthy individuals at single-cell resolution?  
> **Initial status:** Dataset selected. Project locked. No change of topic until completion.

This original charter is retained here as a record of the project's starting scope. The completed project subsequently developed into the broader analytical framework documented above.
