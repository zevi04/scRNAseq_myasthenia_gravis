# Results

This directory contains the primary analytical outputs generated during the **Single-Cell RNA-seq Analysis of Myasthenia Gravis (GSE227835)** project.

The `RESULTS/` directory stores tables, figures, metadata, logs, intermediate analysis objects, quality-control outputs, differential-expression results, functional-enrichment outputs, and integrated downstream analyses produced across the R workflow.

The directory should be interpreted together with the project scripts and documentation rather than as a collection of independent results.

---

## Directory Structure

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
├── Script24_Integrated_Analysis/
└── README.md
```

---

## `Cell_Composition/`

Contains outputs describing the distribution of annotated cell populations across samples or biological groups.

Cell-composition analysis provides a complementary view of the dataset by asking whether the relative representation of annotated populations differs across conditions.

The downstream annotation framework contains **10 cell-type labels**:

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

### Interpretation

Cell-composition differences describe the representation of annotated populations in the analyzed dataset.

They should not automatically be interpreted as evidence of biological expansion, depletion, migration, or clinical severity without appropriate statistical and experimental support.

---

## `DEG/`

Contains outputs from differential gene-expression analysis.

The project evaluates three principal comparisons:

| Comparison | Purpose |
| --- | --- |
| `Healthy_vs_AChR_MG` | AChR-positive MG versus healthy controls |
| `Healthy_vs_SNMG_Pre` | Pre-treatment seronegative MG versus healthy controls |
| `SNMG_Pre_vs_SNMG_Post` | Seronegative MG pre- versus post-treatment |

Differential expression was analyzed in a **cell-type-specific framework**, allowing transcriptional changes to be evaluated within individual annotated populations.

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

### Important interpretation note

The number of upregulated or downregulated genes is not a direct measurement of physiological cellular activation or suppression.

The preferred interpretation is **transcriptional perturbation**, with upward or downward directional bias described where appropriate.

---

## `figures/`

Contains downstream figures generated from result-level analyses.

These figures complement the primary preprocessing and annotation figures stored under the root [`../FIGURES/`](../FIGURES/) directory.

Depending on the analysis stage, this directory may contain visual summaries of:

- differential-expression patterns,
- cell-type comparisons,
- transcriptional direction,
- integrated downstream results,
- reporting outputs.

For the primary QC, PCA, UMAP, clustering, and marker-annotation figures, see [`../FIGURES/README.md`](../FIGURES/README.md).

---

## `Functional_Enrichment/`

Contains outputs generated during functional-enrichment analysis.

Functional enrichment was included to provide pathway- and process-level context for selected transcriptional changes.

Outputs may include enrichment tables and associated visualizations generated from DEG-derived gene sets.

### Interpretation

Enrichment results identify biological terms or pathways that are statistically over-represented within the analyzed gene sets.

They do **not** independently demonstrate:

- pathway activation,
- pathway inhibition,
- causal involvement in MG,
- clinical relevance,
- protein-level activity.

Functional-enrichment results should therefore be interpreted as **hypothesis-generating biological context**.

---

## `logs/`

Contains analysis logs generated during execution of the computational workflow.

Logs provide a record of script execution, progress messages, warnings, completion states, and other diagnostic information.

They are useful for:

- troubleshooting,
- checking whether an analysis stage completed,
- tracing the origin of an output,
- documenting computational execution.

Logs should not be interpreted as biological results.

---

## `Master_Summary/`

Contains consolidated summary outputs produced across the project.

This directory is intended for high-level tables or reports that bring together information from multiple analytical stages.

Depending on the completed workflow, these summaries may include:

- project-level metrics,
- analysis summaries,
- cell-type summaries,
- comparison summaries,
- consolidated biological interpretation tables.

These files are useful entry points when reviewing the overall analysis without opening every individual output directory.

---

## `metadata/`

Contains metadata generated or used during the analysis.

Metadata provide the link between expression data and biological interpretation.

Relevant information may include:

- sample identity,
- GSM/sample identifiers,
- biological group,
- condition,
- cell-level identifiers,
- cluster assignments,
- annotated cell type.

Preserving metadata is essential because the merged expression object alone does not provide all information required to reconstruct the biological comparisons.

For dataset-specific metadata context, see [`../DOCS/Dataset_notes.md`](../DOCS/Dataset_notes.md).

---

## `objects/`

Contains serialized R/Seurat objects or object-related outputs used as computational checkpoints.

Large `.rds` files are excluded from GitHub where appropriate through `.gitignore`.

These objects can be useful locally because the complete dataset is computationally large:

| Stage | Cells |
| --- | ---: |
| Merged dataset before QC | **444,357** |
| Dataset after QC | **420,538** |

Saving intermediate objects can avoid repeating expensive upstream calculations when performing downstream analyses.

### Reproducibility note

Because large serialized objects may not be available in the public repository, full reproduction may require executing the relevant upstream scripts to regenerate them locally.

---

## `qc/`

Contains result-level quality-control outputs.

QC was based on:

- `nFeature_RNA`,
- `nCount_RNA`,
- `percent.mt`.

The final filtering criteria were:

```r
nFeature_RNA > 200
nFeature_RNA < 6000
percent.mt < 15
```

### QC outcome

| Metric | Value |
| --- | ---: |
| Cells before filtering | **444,357** |
| Cells after filtering | **420,538** |
| Cells removed | **23,819** |
| Percentage removed | **5.36%** |

Primary QC visualizations are also organized under [`../FIGURES/qc/`](../FIGURES/qc/).

---

# `Script24_Integrated_Analysis/`

This directory contains the major integrated downstream outputs produced by the Script 24 analysis series.

The integrated analysis was divided across multiple scripts rather than executed as one monolithic file.

Broadly, this stage consolidated cell-type-specific DEG results and generated higher-level summaries for comparison across the complete project.

Outputs from this analysis include or support:

- global DEG summaries,
- cell-type DEG summaries,
- comparison-level DEG summaries,
- upregulated versus downregulated summaries,
- significant DEG visualization,
- shared and unique DEG analysis,
- gene × cell-type occurrence analysis,
- gene × comparison occurrence analysis,
- gene-recurrence summaries,
- conserved disease-signature prioritization,
- publication/reporting-oriented downstream outputs.

---

## Global DEG Analysis

The integrated DEG workflow combines differential-expression information across cell types and comparisons.

This makes it possible to compare:

```text
Cell type
    ×
Biological comparison
    ×
Gene
    ×
Direction
    ×
Effect size
    ×
Adjusted significance
```

rather than interpreting each DEG table independently.

---

## Shared and Unique DEG Analysis

Genes were evaluated according to how many annotated cell populations contained a corresponding DEG record.

A gene-occurrence framework was constructed using:

**CellType_Frequency**

> Number of cell types in which a gene occurs.

and:

**Comparison_Frequency**

> Number of biological comparisons in which a gene occurs.

This enabled classification of genes into categories such as:

- cell-type specific,
- shared across two cell types,
- highly shared,
- conserved across represented cell types.

### Interpretation

A recurrent gene may represent a more broadly distributed transcriptional alteration, while a low-frequency gene may represent a more cell-type-restricted signal.

Recurrence does not by itself establish biological importance, disease specificity, or causality.

---

## Conserved Disease Signatures

The downstream workflow further evaluated recurrent genes using gene-level statistics incorporating information such as:

- mean log2 fold-change,
- mean absolute log2 fold-change,
- adjusted significance,
- cell-type recurrence,
- comparison recurrence.

The purpose of this analysis was to prioritize genes showing stronger and more recurrent transcriptional patterns.

The resulting signature rankings should be interpreted as **computational prioritization**, not validated MG biomarkers.

---

## Directional Cell-Type Results

Later downstream analysis compared upregulated and downregulated significant DEG records across cell types and conditions.

Several major patterns emerged from the completed analysis.

### AChR-positive MG vs Healthy

| Cell type | Upregulated | Downregulated |
| --- | ---: | ---: |
| Neutrophils | **27** | 12 |
| B cells | 9 | **50** |
| Dendritic cells | 2 | **37** |
| CD4 T cells | 1 | **27** |
| NK cells | 3 | **24** |

Neutrophils showed a **69.2% upward directional bias** among significant DEG entries in this comparison.

---

### Seronegative MG Pre-treatment vs Healthy

| Cell type | Upregulated | Downregulated |
| --- | ---: | ---: |
| Neutrophils | **39** | 1 |
| B cells | 30 | **107** |
| Dendritic cells | 20 | **80** |
| NK cells | 9 | **40** |
| Monocytes | **59** | **57** |

Neutrophils showed a **97.5% upward directional bias**.

Monocytes displayed a substantial DEG burden but an approximately balanced directional distribution.

---

### Seronegative MG Pre vs Post

| Cell type | Upregulated | Downregulated |
| --- | ---: | ---: |
| B cells | **7** | 0 |
| NK cells | **17** | 5 |
| Dendritic cells | **13** | 6 |
| Neutrophils | 8 | **52** |
| CD8 T cells | 1 | **20** |
| T cells | 3 | **20** |

Neutrophils showed an **86.7% downward directional bias** in the pre/post comparison.

### Overall interpretation

The direction and magnitude of transcriptional perturbation varied substantially by cell type and biological comparison.

The data therefore do not support a single uniform immune transcriptional response across MG-associated conditions.

---

## AChR / MuSK Targeted Results

A targeted downstream analysis examined genes associated with the acetylcholine receptor, MuSK, and selected neuromuscular-junction components.

The target panel included genes such as:

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

Among the examined AChR receptor-subunit genes, only `CHRNB1` appeared in the master DEG table.

The observed record was:

| Feature | Observation |
| --- | --- |
| Gene | `CHRNB1` |
| Cell type | Progenitors |
| Comparison | Healthy vs AChR-MG |
| avg_log2FC | approximately −0.535 |
| Adjusted P value | 1.0 |
| Significant by final thresholds | **No** |

`MUSK` was not identified as a significant DEG.

### Targeted summary

| Target | Significant DEG records |
| --- | ---: |
| AChR receptor-subunit genes | **0** |
| `MUSK` | **0** |

This analysis therefore does not support significant AChR receptor-subunit or `MUSK` transcript overexpression in the analyzed immune-cell DEG framework.

AChR-positive MG is an **autoantibody-defined disease category** and should not be interpreted as evidence of AChR-gene overexpression.

---

# How to Interpret the Results Directory

A useful reading order is:

```text
1. qc/
      ↓
2. metadata/
      ↓
3. Cell_Composition/
      ↓
4. DEG/
      ↓
5. Functional_Enrichment/
      ↓
6. Script24_Integrated_Analysis/
      ↓
7. Master_Summary/
      ↓
8. figures/
```

The `logs/` and `objects/` directories primarily support reproducibility and troubleshooting rather than biological interpretation.

---

## Result Interpretation Principles

### Statistical significance is not biological importance

A statistically significant DEG is not automatically a biologically important or disease-driving gene.

### Fold change should be interpreted with significance and context

Effect size, adjusted significance, recurrence, cell identity, and biological comparison should be considered together.

### DEG count is not cell activity

A cell type with many upregulated genes is not automatically physiologically activated. Likewise, many downregulated genes do not prove functional suppression.

### Recurrence is not causality

Genes occurring across multiple cell types may be useful candidates for further investigation, but recurrence does not establish a causal MG mechanism.

### Pre/post difference is not automatically recovery

The `SNMG_Pre_vs_SNMG_Post` comparison identifies transcriptional differences between sampled states. Demonstrating normalization would require showing movement of disease-associated genes toward the healthy-control state.

### Transcript abundance is not protein abundance

scRNA-seq measures RNA. It does not directly quantify autoantibodies, receptor abundance, protein activity, signaling activity, or neuromuscular-junction function.

---

## File Formats

Depending on the analysis stage, result directories may contain:

| Format | Typical use |
| --- | --- |
| `.csv` | Tables, summaries, metadata and DEG outputs |
| `.png` | Directly viewable figures |
| `.pdf` | High-resolution/vector-oriented figures |
| `.rds` | Serialized R/Seurat objects |
| `.txt` / log files | Execution and diagnostic records |

Large `.rds` files are intentionally excluded from GitHub where appropriate.

---

## Reproducibility

The results in this directory were generated programmatically from the numbered R workflow under:

[`../SCRIPTS/`](../SCRIPTS/)

The broad workflow is:

```text
Data preparation
      ↓
Quality control
      ↓
Preprocessing
      ↓
PCA / UMAP / clustering
      ↓
Marker identification
      ↓
Cell-type annotation
      ↓
Differential expression
      ↓
Cell composition + functional enrichment
      ↓
Integrated DEG analysis
      ↓
Recurrence + conserved signatures
      ↓
Directional analysis
      ↓
AChR / MuSK targeted analysis
```

For exact script order and execution notes, see [`../SCRIPTS/README.md`](../SCRIPTS/README.md).

---

## Related Documentation

| Document | Purpose |
| --- | --- |
| [`../README.md`](../README.md) | Main project overview and key findings |
| [`../PROJECT.md`](../PROJECT.md) | Scientific rationale, objectives, design and scope |
| [`../JOURNAL.md`](../JOURNAL.md) | Chronological research and troubleshooting record |
| [`../DOCS/Dataset_notes.md`](../DOCS/Dataset_notes.md) | Dataset-specific reference |
| [`../FIGURES/README.md`](../FIGURES/README.md) | Primary figure guide |
| [`../SCRIPTS/README.md`](../SCRIPTS/README.md) | Script execution and workflow guide |

---

## Project Status

**Completed — 3 August 2026**

The `RESULTS/` directory contains the analytical outputs supporting the completed cell-type-specific transcriptional characterization and downstream descriptive analysis of Myasthenia Gravis using GSE227835.
