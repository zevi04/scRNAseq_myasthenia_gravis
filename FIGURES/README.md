# Figures

This directory contains the principal visual outputs generated during the **Single-Cell RNA-seq Analysis of Myasthenia Gravis (GSE227835)** project.

The figures document the progression of the analysis from quality control and dimensionality reduction through clustering and cell-type annotation. Additional downstream figures generated during differential-expression, recurrence, functional-enrichment, and integrated analyses are stored under the corresponding `RESULTS/` directories.

---

## Directory Structure

```text
FIGURES/
├── CLUSTERING/
├── PCA/
├── qc/
├── UMAP/
├── Dotplot_CanonicalMarkers.png
├── Dotplot_TopMarkers.png
├── Heatmap_TopMarkers.png
└── README.md
```

---

## Figure Categories

### `qc/` — Quality-Control Figures

Contains visualizations used to assess the quality of the merged single-cell dataset before and after filtering.

The QC workflow evaluated standard Seurat metrics including:

- `nFeature_RNA` — number of detected genes per cell,
- `nCount_RNA` — total RNA counts per cell,
- `percent.mt` — percentage of mitochondrial transcripts.

The final filtering criteria were:

```r
nFeature_RNA > 200
nFeature_RNA < 6000
percent.mt < 15
```

The merged dataset contained **444,357 cells before QC** and **420,538 cells after filtering**.

QC figures should be interpreted as diagnostic visualizations used to assess cell quality and the consequences of the filtering thresholds.

---

### `PCA/` — Principal Component Analysis

Contains figures generated during Principal Component Analysis (PCA).

PCA was used to reduce the dimensionality of the normalized and scaled expression matrix before graph construction and clustering.

The project calculated **50 principal components** for downstream dimensionality-reduction and neighborhood analyses.

PCA figures provide a global view of major transcriptional variation in the dataset and support selection and evaluation of dimensions used in subsequent analyses.

---

### `UMAP/` — UMAP Visualizations

Contains Uniform Manifold Approximation and Projection (UMAP) figures generated during the analysis.

UMAP provides a two-dimensional representation of transcriptional similarity among cells.

Depending on the analytical stage, UMAP figures may display cells according to:

- computational cluster,
- sample,
- biological condition,
- annotated cell type.

UMAP coordinates are a visualization of high-dimensional transcriptional structure and should not be interpreted as direct physical or developmental distances between cells.

---

### `CLUSTERING/` — Clustering Figures

Contains figures related to graph-based clustering and visualization of transcriptionally distinct cell populations.

The unsupervised analysis initially identified **34 transcriptional clusters**.

These clusters were subsequently evaluated using marker-gene expression and consolidated into **10 downstream cell-type labels**:

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

Clustering figures represent computationally identified transcriptional populations and should be interpreted together with marker-based annotation.

---

## Marker-Based Annotation Figures

Three principal annotation figures are stored directly in the `FIGURES/` directory.

### `Dotplot_CanonicalMarkers.png`

Displays expression patterns of canonical marker genes across clusters or annotated populations.

This figure was used to support biological interpretation of the computational clusters by comparing observed expression with expected immune-cell marker patterns.

The dot-plot representation generally combines information about:

- the proportion of cells expressing a marker,
- the relative expression level of that marker.

Canonical-marker expression should be considered collectively rather than assigning a cell identity from a single gene.

---

### `Dotplot_TopMarkers.png`

Displays selected top marker genes associated with the transcriptional clusters.

This figure complements the canonical-marker analysis by showing genes identified from the dataset itself as prominent cluster-associated markers.

It provides an additional layer of evidence for distinguishing transcriptionally related populations and evaluating proposed annotations.

---

### `Heatmap_TopMarkers.png`

Displays expression patterns of selected top marker genes across clusters or cell populations in heatmap form.

The heatmap provides a broader view of marker-expression structure and helps identify:

- cluster-specific expression patterns,
- shared marker programs,
- transcriptional similarity between clusters,
- marker sets supporting downstream annotation.

The heatmap should be interpreted alongside the dot plots and marker tables rather than as an independent annotation method.

---

## Relationship Between Clusters and Cell Types

The project distinguishes between **computational clusters** and **final downstream cell-type labels**.

```text
High-dimensional expression data
        ↓
Graph-based clustering
        ↓
34 transcriptional clusters
        ↓
Marker-gene evaluation
        ↓
Biological annotation
        ↓
10 downstream cell-type labels
```

The 34 clusters therefore do not represent 34 final biological cell types.

Multiple transcriptionally related clusters can contribute to the same broader annotated cell population.

---

## Downstream Figures

Not every figure generated by the project is stored in the root `FIGURES/` directory.

Later analysis stages generated additional figures under `RESULTS/`, particularly within:

```text
RESULTS/
├── Cell_Composition/
├── DEG/
├── figures/
├── Functional_Enrichment/
├── qc/
└── Script24_Integrated_Analysis/
```

These downstream outputs include visualizations associated with:

- immune-cell composition,
- differential-expression summaries,
- upregulated versus downregulated DEG patterns,
- comparison-level DEG distributions,
- functional enrichment,
- gene recurrence,
- shared and unique DEGs,
- conserved transcriptional signatures,
- cell-type directional analysis,
- AChR/MuSK targeted analysis,
- integrated reporting.

See [`../RESULTS/README.md`](../RESULTS/README.md) for the organization and interpretation of downstream outputs.

---

## Figure Interpretation Principles

### UMAP proximity is not physical distance

Cells appearing close together on UMAP have similar representations in the reduced transcriptional space. UMAP geometry should not be interpreted as anatomical distance, lineage proof, or direct cell-cell interaction.

### Clusters are computational populations

Cluster numbers are generated from graph-based analysis. Biological identities require marker-based interpretation.

### Marker annotation should use multiple genes

A single marker gene is rarely sufficient to establish a robust cell identity. Annotation was therefore supported using canonical-marker patterns, top markers, dot plots, heatmaps, and broader transcriptional context.

### DEG counts are not direct measures of cellular activity

Downstream figures showing upregulated and downregulated DEG counts describe transcriptional direction and perturbation. They should not automatically be interpreted as physiological activation or suppression.

### Heatmap intensity is context dependent

Heatmap values may represent scaled expression, fold change, occurrence, or another analysis-specific quantity. The corresponding legend and source analysis should always be checked before biological interpretation.

### Pre/post differences are not automatically normalization

Figures comparing seronegative MG pre- and post-treatment states demonstrate transcriptional differences between the sampled states. They do not independently establish treatment-induced recovery or normalization toward healthy controls.

---

## Recommended Figure Reading Order

For readers following the analytical workflow, the figures are best examined in approximately this order:

```text
1. QC figures
       ↓
2. PCA
       ↓
3. Clustering
       ↓
4. UMAP
       ↓
5. Canonical-marker dot plot
       ↓
6. Top-marker dot plot
       ↓
7. Top-marker heatmap
       ↓
8. Cell-composition figures
       ↓
9. Differential-expression figures
       ↓
10. Functional-enrichment figures
       ↓
11. Gene-recurrence / conserved-signature figures
       ↓
12. Directional and AChR/MuSK downstream figures
```

This sequence follows the logic of the computational workflow from data quality to biological interpretation.

---

## Figure Formats

The repository primarily uses image formats suitable for GitHub viewing and downstream reporting.

Depending on the analysis stage, figures may be available as:

- `.png` — raster images suitable for direct viewing,
- `.pdf` — vector or publication-oriented output where generated.

Where both formats exist, the PNG version is convenient for rapid inspection while the PDF version is generally preferable for high-quality figure assembly or publication-oriented use.

---

## Reproducibility

Figures were generated programmatically from the R analysis workflow rather than manually reconstructed.

The relevant code is maintained under:

[`../SCRIPTS/`](../SCRIPTS/)

For script sequence and purpose, see:

[`../SCRIPTS/README.md`](../SCRIPTS/README.md)

Because figures are derived from intermediate and final analysis objects, reproducing every figure may require the corresponding local Seurat `.rds` objects or upstream scripts. Large serialized objects are excluded from GitHub where appropriate.

---

## Related Documentation

| Document | Purpose |
| --- | --- |
| [`../README.md`](../README.md) | Main project overview and key results |
| [`../PROJECT.md`](../PROJECT.md) | Scientific rationale and analytical design |
| [`../JOURNAL.md`](../JOURNAL.md) | Chronological analysis and troubleshooting record |
| [`../DOCS/Dataset_notes.md`](../DOCS/Dataset_notes.md) | Dataset-specific reference |
| [`../SCRIPTS/README.md`](../SCRIPTS/README.md) | Script execution guide |
| [`../RESULTS/README.md`](../RESULTS/README.md) | Downstream output guide |

---

## Project Status

**Completed — 3 August 2026**

The figures in this directory document the primary preprocessing, dimensionality-reduction, clustering, and cell-annotation stages of the completed GSE227835 Myasthenia Gravis single-cell transcriptomic analysis.
