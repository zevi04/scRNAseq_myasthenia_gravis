# RESEARCH JOURNAL

## Single-Cell RNA-seq Analysis of Myasthenia Gravis --- GSE227835

**Project period:** 22 July 2026 -- 3 August 2026\
**Repository:** `scRNAseq_myasthenia_gravis`\
**Primary framework:** R / Seurat

> This journal documents the computational progression, technical
> challenges, analytical decisions, observations, and final
> interpretation of the project.

------------------------------------------------------------------------

# 22 July

Today's objective - Repository setup

Completed - Created GitHub repository - Locked dataset - Defined
research question

Questions - What is the sample composition of GSE227835?

Tomorrow - Download processed files - Read associated publication

# 23 July

## Objective

Set up the computational workflow for the single-cell RNA-seq analysis
of Myasthenia Gravis (GSE227835) and inspect the processed data.

## Completed

-   Created and published the GitHub repository.
-   Defined the research question and locked the project scope.
-   Created the project structure:
    -   DATA/
    -   DOCS/
    -   FIGURES/
    -   NOTEBOOKS/
    -   RESULTS/
    -   SCRIPTS/
-   Added `PROJECT.md`, `README.md`, and `JOURNAL.md`.
-   Created a `.gitignore` to exclude downloaded datasets from version
    control.
-   Downloaded the processed GEO dataset (GSE227835).
-   Extracted the supplementary archive into `DATA/PROCESSED/`.
-   Examined the study design and identified the sample groups:
    -   Healthy Controls
    -   AChR-positive Myasthenia Gravis
    -   Seronegative MG (Pre-treatment)
    -   Seronegative MG (Post-treatment)

## Observations

-   The processed dataset consists of multiple compressed `.txt.gz`
    files.
-   The files are not immediately readable using `fread()` as initially
    expected.
-   A preview of one file showed barcode-like entries rather than a
    standard gene × cell expression matrix.
-   The processed data format needs to be understood before importing it
    into R.

## Challenges

-   Assumed the processed files were simple expression matrices.
-   Initial attempts to read the files in R were unsuccessful.
-   Realized that understanding the file structure is the necessary next
    step before writing analysis code.

## Lessons Learned

-   Never assume the format of a public dataset.
-   Inspect the dataset structure before choosing import functions.
-   Good computational biology begins with understanding the data, not
    writing code.

## Next Steps

-   Determine the exact structure of the processed files.
-   Identify what each `.txt.gz` file represents.
-   Build the sample metadata table.
-   Import the dataset into R using the appropriate workflow.

# 26 JULY

## Progress

Today marked the completion of the data preparation phase of the
scRNA-seq pipeline.

### Completed

-   Successfully created **40 individual Seurat objects** from the
    processed GEO expression matrices.
-   Merged all samples into a single Seurat object.

### Merged Dataset

-   **Samples:** 40
-   **Cells:** 444,357
-   **Genes:** 36,601

The merged object was saved to: RESULTS/objects/merged_seurat_object.rds

### Quality Control

Calculated standard QC metrics for every cell:

-   `nFeature_RNA`
-   `nCount_RNA`
-   `percent.mt`

QC summary:

  Metric            Min   Median      Max
  -------------- ------ -------- --------
  nFeature_RNA       12     1474     8380
  nCount_RNA        500     5168   107729
  percent.mt       0.00     7.37    98.26

### QC Visualization

Generated and saved:

-   QC Violin Plot
-   Counts vs Features
-   Counts vs Mitochondrial %
-   Features vs Mitochondrial %

Location: FIGURES/qc/

Initial violin plots were cluttered because \~444k cells were plotted
with individual points. Updated the visualization by removing point
overlays (`pt.size = 0`), resulting in cleaner distributions.

### Observations

-   Strong positive correlation between UMI counts and detected genes
    (\~0.81).
-   Most cells have mitochondrial percentages below \~10%.
-   A subset of cells shows very high mitochondrial content, suggesting
    low-quality or dying cells that will be addressed during filtering.

## Next Steps

-   Determine QC filtering thresholds.
-   Remove low-quality cells.
-   Normalize the filtered dataset.
-   Begin downstream dimensionality reduction and clustering.

# 27 JULY

Work Completed 1. Quality Control Filtering

Applied filtering criteria to remove low-quality cells.

Filtering thresholds

nFeature_RNA \> 200 nFeature_RNA \< 6000 percent.mt \< 15%

Results

Cells before filtering: 444,357 Cells after filtering: 420,538 Cells
removed: 23,819 Percentage removed: 5.36%

The removal percentage indicates that the dataset was already of good
quality while eliminating cells likely to introduce technical noise.

2.  Data Normalization

Performed library-size normalization using Seurat's NormalizeData()
function.

Method

LogNormalize Scale factor = 10,000

Normalization adjusts gene expression values so that differences in
sequencing depth between cells do not bias downstream analyses.

3.  Highly Variable Gene Selection

Identified genes with the greatest biological variability using
FindVariableFeatures().

Method

Selection: VST Number of genes selected: 2,000

These genes will be used for dimensionality reduction because they
contribute most to distinguishing different cell populations.

4.  Data Scaling

Scaled the highly variable genes using ScaleData().

Scaling standardizes each gene to:

Mean = 0 Standard deviation = 1

This prevents highly expressed genes from dominating Principal Component
Analysis.

Conceptual Understanding Gained

Today's session focused heavily on understanding the reasoning behind
each preprocessing step instead of simply executing Seurat functions.

Key concepts learned Quality control removes poor-quality cells, not
genes. Normalization corrects differences in sequencing depth between
cells. Highly variable genes capture the biological differences between
cells. Scaling standardizes expression values across genes before
dimensionality reduction. PCA does not classify or cluster cells; it
compresses thousands of correlated gene expression measurements into a
smaller number of informative principal components. UMAP is a
visualization method that projects cells into two dimensions while
preserving neighborhood relationships. Cell clustering occurs after PCA
and UMAP and is based on similarity between cells. Cell annotation is a
biological interpretation step performed using known marker genes.
Project Status Completed

✓ Dataset Exploration ✓ Data Import ✓ Metadata Construction ✓ Dataset
Consistency Check ✓ Seurat Object Creation ✓ Merge Seurat Objects ✓
Quality Control Metrics ✓ QC Visualization ✓ QC Filtering ✓
Normalization ✓ Highly Variable Gene Selection ✓ Scaling Repository
Improvements

During today's work, the project structure was reorganized for improved
readability.

Current organization:

RESULTS/ ├── metadata/ ├── objects/ ├── qc/ ├── tables/ ├── figures/ └──
logs/

This organization will make future analyses and figure management
easier.

Challenges Encountered

Several technical issues were encountered and resolved.

Corrupted intermediate Seurat object caused by an incomplete saveRDS()
operation. Initial confusion regarding Seurat object storage and large
file sizes. Realization that intermediate Seurat objects occupy
substantial disk space (approximately 12 GB total). Discussed future
optimization by retaining only essential checkpoint objects after the
pipeline is finalized. Reflections

An important realization today was that simply running scripts does not
necessarily lead to understanding the analysis. Moving forward, the
workflow will emphasize understanding the biological and computational
rationale behind each step before implementing it in code.

This approach is intended to improve confidence in explaining the
pipeline during interviews, research discussions, and future projects.

Next Steps

The next stage begins dimensionality reduction.

Planned analyses:

Principal Component Analysis (PCA) Elbow Plot generation PCA loading
visualization Selection of informative principal components UMAP
Clustering Cell type annotation Personal Reflection

Today marked a shift from simply executing a bioinformatics pipeline to
understanding how each analytical step contributes to interpreting
single-cell RNA-seq data. The focus moved away from memorizing Seurat
functions toward understanding the underlying concepts, with the goal of
being able to explain the complete workflow confidently in a research
interview or laboratory discussion.

# 28 JULY

**Session:** Differential Expression & Marker Gene Identification
Troubleshooting

------------------------------------------------------------------------

# Objective

To identify cluster-specific marker genes following successful
clustering of the integrated scRNA-seq dataset.

------------------------------------------------------------------------

# Work Completed

### 1. PCA Completed

-   Successfully performed PCA using 50 principal components.

-   Generated:

    -   Elbow Plot
    -   PCA Heatmap

-   Saved PCA object for downstream analysis.

------------------------------------------------------------------------

### 2. UMAP Completed

Successfully generated the UMAP embedding using the selected principal
components.

Outputs:

-   UMAP projection
-   Publication-quality TIFF figure

------------------------------------------------------------------------

### 3. Nearest Neighbor Graph

Constructed the Shared Nearest Neighbor (SNN) graph using:

-   FindNeighbors()
-   Leiden clustering algorithm

No issues encountered after installing the required Leiden dependencies.

------------------------------------------------------------------------

### 4. Cell Clustering

Successfully clustered the complete dataset.

Dataset statistics:

-   **Genes:** 36,601
-   **Cells:** 420,538
-   **Clusters identified:** 34

Largest cluster:

-   Cluster 1
-   65,031 cells

Smallest cluster:

-   Cluster 34
-   307 cells

Cluster UMAP generated successfully.

------------------------------------------------------------------------

# Major Issue 1

## Corrupted PCA Object

The PCA object became corrupted after an interrupted save.

Error:

    error reading from connection

### Resolution

The PCA step was rerun from the previous checkpoint.

Lesson learned:

> Never interrupt `saveRDS()` while it is writing large Seurat objects.

------------------------------------------------------------------------

# Major Issue 2

## Seurat v5 Layer Architecture

Initial execution of:

``` r
FindAllMarkers()
```

returned:

    No DE genes identified

    data layers are not joined.
    Please run JoinLayers

Investigation showed that the merged Assay5 object contained:

    counts.A1
    counts.A2
    ...
    counts.N9b

    data.A1
    data.A2
    ...
    data.N9b

Attempting:

``` r
JoinLayers(scaled_seurat)
```

failed because Windows could not allocate approximately 3.5 GB of
additional memory.

------------------------------------------------------------------------

## Resolution

Instead of joining the entire Seurat object, only the normalized RNA
data layer was joined:

``` r
scaled_seurat[["RNA"]] <- JoinLayers(
    object = scaled_seurat[["RNA"]],
    layers = "data",
    new = "data"
)
```

Verification:

Before:

    data.A1
    ...
    data.N9b

After:

    data
    counts.A1
    ...
    counts.N9b

This successfully resolved the Seurat v5 layer issue.

------------------------------------------------------------------------

# Major Issue 3

## Differential Expression Memory Limitation

Despite resolving the layer issue,

    FindAllMarkers()

failed with:

    cannot allocate vector of size 2.3–2.6 GB

The error occurred during:

    rowSums()

inside the differential expression workflow.

Testing the smallest cluster (307 cells) with:

``` r
FindMarkers()
```

produced the same error, indicating that Seurat constructs a large
temporary matrix involving the entire dataset regardless of cluster
size.

Conclusion:

The limitation is due to available system memory rather than incorrect
preprocessing.

------------------------------------------------------------------------

# Additional Investigation

Attempted to install the **Presto** package to enable a more
memory-efficient Wilcoxon implementation.

Installation failed because build tools were unavailable.

Diagnostic checks showed:

    where make
    → not found

    where gcc
    → not found

indicating that Rtools was not correctly configured for compilation.

------------------------------------------------------------------------

# Storage Management

Due to limited storage on the C: drive, intermediate Seurat checkpoints
were removed.

Retained:

-   `merged_seurat_object.rds`
-   `clustered_seurat_object.rds`

Deleted:

-   Filtered object
-   Normalized object
-   Variable features object
-   Scaled object
-   PCA object
-   Other intermediate checkpoints

This reduced disk usage while preserving the ability to restart from
both the merged dataset and the clustered dataset.

------------------------------------------------------------------------

# Key Lessons Learned

-   Seurat v5 stores merged samples as multiple assay layers.
-   Differential expression requires normalized data layers to be
    joined.
-   Joining only the RNA `data` layer is substantially more
    memory-efficient than joining the entire Seurat object.
-   Large Seurat objects (\>400,000 cells) can exceed the memory
    available on standard desktop hardware during differential
    expression analysis.
-   Saving checkpoints after major computational steps prevents loss of
    work following unexpected failures.

------------------------------------------------------------------------

# Next Steps

-   Evaluate a memory-efficient strategy for marker identification
    (e.g., downsampling, alternative differential expression methods, or
    reference-based annotation).
-   Complete cell type annotation.
-   Validate annotations using canonical marker genes.
-   Document the Seurat v5 `JoinLayers()` solution and memory
    considerations in the project README.

------------------------------------------------------------------------

## Personal Reflection

Today's session was dominated by debugging rather than generating
biological results, but it resolved one of the most significant
technical challenges encountered so far. The preprocessing
pipeline---from quality control through clustering---is complete and
reproducible. The remaining challenge is selecting a differential
expression strategy that is compatible with the available hardware while
maintaining biological validity. This marks a transition from building
the computational pipeline to extracting and interpreting biological
insights.

------------------------------------------------------------------------

# 29 JULY

## Objective

Move from unsupervised clustering toward biologically interpretable
immune-cell identities and establish the framework for
cell-type-specific disease analysis.

## Work Completed

The clustered dataset was carried forward into cell-type annotation and
downstream interpretation. The 34 computational clusters were
consolidated into biologically interpretable immune-cell populations
using marker-expression information and canonical immune-cell
identities.

The final downstream analysis used 10 cell-type labels:

-   B cells
-   Basophils
-   CD4 T cells
-   CD8 T cells
-   Dendritic cells
-   Monocytes
-   Neutrophils
-   NK cells
-   Progenitors
-   T cells

This was an important transition in the project. Cluster numbers
describe transcriptionally similar groups, but they are not themselves
biological identities. Annotation made it possible to ask which
immune-cell populations were associated with the transcriptional
alterations observed in MG.

## Analytical Direction

The project was reorganized around cell-type-specific comparisons rather
than treating the complete dataset as one homogeneous population. This
was necessary because different immune-cell populations have
fundamentally different baseline transcriptomes.

## Next Steps

-   Finalize condition metadata.
-   Define disease comparisons.
-   Perform differential-expression analysis within cell types.
-   Integrate DEG outputs into a common downstream framework.

------------------------------------------------------------------------

# 30 JULY

## Objective

Finalize the biological comparison framework and begin systematic
cell-type-specific differential-expression analysis.

## Study Groups

Four biological groups were retained:

-   Healthy Controls
-   AChR-positive Myasthenia Gravis
-   Seronegative MG --- Pre-treatment
-   Seronegative MG --- Post-treatment

## Main Comparisons

Three comparisons were used for downstream analysis:

``` text
Healthy_vs_AChR_MG
Healthy_vs_SNMG_Pre
SNMG_Pre_vs_SNMG_Post
```

These comparisons answer different biological questions.

**Healthy vs AChR-MG** examines transcriptional differences associated
with AChR-positive MG relative to healthy controls.

**Healthy vs SNMG Pre** examines the transcriptional state of
seronegative MG before the post-treatment time point relative to healthy
controls.

**SNMG Pre vs SNMG Post** examines transcriptional changes between the
two seronegative MG time points.

## Differential-Expression Framework

The downstream workflow retained gene identity, cell type, comparison,
average log2 fold-change, adjusted P value, and direction of change.

For final reporting, significant transcriptional changes were defined
using:

``` text
Adjusted P < 0.05
|average log2 fold-change| >= 0.5
```

Upregulated and downregulated DEGs were considered separately.

## Important Interpretation

A cell type with more upregulated genes was not automatically described
as "more active." DEG direction describes transcriptional change, not
direct physiological activation. The term **transcriptional
perturbation** was therefore adopted for broader comparisons.

------------------------------------------------------------------------

# 31 JULY

## Objective

Integrate differential-expression results across cell types and disease
comparisons and build a reproducible downstream analysis.

## Work Completed

A master DEG framework was assembled so that results from different cell
types and comparisons could be analyzed together.

Summary outputs were generated for:

-   total DEG burden by cell type
-   total DEG burden by comparison
-   upregulated and downregulated genes
-   significant DEG distributions
-   cell-type × comparison summaries

This enabled the project to move beyond individual DEG tables and ask
broader questions about the structure of the MG-associated
transcriptional response.

## Shared and Unique DEG Analysis

A gene-by-cell-type occurrence matrix was constructed.

For every gene, two recurrence measures were calculated:

-   **CellType_Frequency** --- number of cell types in which the gene
    occurred
-   **Comparison_Frequency** --- number of comparisons in which the gene
    occurred

Genes were classified as:

-   Cell-type specific
-   Shared across two cell types
-   Highly shared
-   Conserved across all represented cell types

The recurrence analysis showed that many significant DEGs occurred in
only a limited number of immune-cell populations, while progressively
fewer genes were shared across many cell types.

This suggested that a substantial part of the MG-associated
transcriptional response in this dataset was **cell-type dependent
rather than universally conserved**.

## Conserved Disease Signatures

A conserved-signature analysis was developed to prioritize recurrent
genes using information from:

-   average absolute log2 fold-change
-   adjusted statistical significance
-   cell-type recurrence

## Technical Issue

Joining gene-level statistics with occurrence statistics produced
duplicate frequency columns:

``` text
CellType_Frequency.x
CellType_Frequency.y
Comparison_Frequency.x
Comparison_Frequency.y
```

An earlier version of the script attempted to use a non-existent
unsuffixed `CellType_Frequency` column.

The section was rewritten to explicitly handle the duplicated columns.

## Visualization Review

Several early heatmaps were technically generated but were not
scientifically useful because of:

-   excessive numbers of genes
-   unreadable labels
-   sparse matrices
-   large blank regions
-   clustering instability

The visualization strategy was revised to prioritize interpretable
summaries rather than retaining every automatically generated plot.

## Lesson Learned

A figure being successfully generated does not mean that it communicates
the biology effectively.

------------------------------------------------------------------------

# 1 AUGUST

## Objective

Improve downstream visualization, reporting, figure legends, methodology
documentation, and publication-style organization.

## Work Completed

The integrated-analysis outputs were reviewed for clarity and
reproducibility.

The reporting stage included:

-   global DEG visualization
-   comparison-level DEG summaries
-   upregulated-versus-downregulated plots
-   recurrence analysis
-   conserved-signature visualization
-   multi-panel figure assembly
-   figure legends
-   methodology documentation
-   manuscript-oriented supplementary tables
-   biological interpretation outputs

## Reporting Problems Identified

Some figures were difficult to interpret even though the underlying
analysis had completed.

Examples included:

-   a gene-recurrence heatmap that appeared effectively empty
-   a global DEG heatmap that was visually overcrowded
-   conserved-signature plots with excessive unused space
-   figure panels that lacked sufficient explanation
-   legends that did not clearly state what colors or values represented

The recurrence visualization was therefore shifted toward a clearer
frequency-based representation showing how many genes were detected
across increasing numbers of cell types.

## Methodology Structure

The computational workflow was organized conceptually as:

``` text
GEO dataset acquisition
        ↓
Sample-level Seurat objects
        ↓
Merge and metadata integration
        ↓
Quality control
        ↓
Normalization
        ↓
Highly variable gene selection
        ↓
Scaling
        ↓
PCA
        ↓
Neighbor graph + Leiden clustering
        ↓
UMAP
        ↓
Cell-type annotation
        ↓
Cell-type-specific differential expression
        ↓
Integrated DEG analysis
        ↓
Shared / unique / recurrent DEG analysis
        ↓
Conserved-signature analysis
        ↓
Directional biological interpretation
```

## Reflection

This stage emphasized that scientific reporting is part of the analysis
itself. Clear figures, explicit legends, and defensible terminology are
necessary for another reader to understand what was actually measured.

------------------------------------------------------------------------

# 2 AUGUST

## Objective

Repair remaining downstream-reporting problems and make the analysis
robust enough for final interpretation.

## Downstream Reporting Repair

A dedicated repair script was used to correct problems in the downstream
reporting workflow.

The work focused on:

-   improving DEG summary figures
-   repairing recurrence and conserved-signature outputs
-   improving figure assembly
-   generating clearer legends
-   improving metadata/reporting outputs
-   organizing final result files

## Error --- `slice_head()`

A script used `dplyr::n()` inside the `n` argument of `slice_head()`.

This failed because `slice_head()` requires a constant value for `n`,
whereas `dplyr::n()` is intended for data-masking operations.

The affected sections were corrected by determining the required number
of rows before calling `slice_head()`.

## Error --- Heatmap Clustering

A later heatmap failed during hierarchical clustering with:

``` text
NA/NaN/Inf in foreign function call
```

This highlighted the need to validate matrices for missing and
non-finite values before clustering.

## Final Figure Philosophy

The figures were reviewed according to whether an independent reader
could understand:

-   what was being compared
-   what the axes represented
-   what the colors meant
-   whether values represented DEG counts, direction, recurrence, or
    effect size
-   what biological conclusion could and could not be drawn

This reduced reliance on visually complex but poorly interpretable
heatmaps.

------------------------------------------------------------------------

# 3 AUGUST

## Objective

Answer the main biological questions that emerged from the project and
define a clear analytical endpoint.

## Final Downstream Analysis --- Script 26

The final analysis focused on four questions:

1.  Which immune-cell populations show the greatest significant
    transcriptional perturbation?
2.  Which populations are predominantly upregulated or downregulated?
3.  How do these patterns differ across the three MG comparisons?
4.  Is there transcript-level evidence for significant AChR- or
    MuSK-related target-gene changes?

The analysis used:

``` text
Adjusted P < 0.05
|average log2 fold-change| >= 0.5
10 cell types
3 comparisons
```

## Healthy vs AChR-positive MG

Neutrophils showed the strongest upward DEG profile:

-   27 upregulated
-   12 downregulated
-   39 significant DEG entries
-   69.2% upward directional bias

Several other populations showed predominantly downward changes.

**B cells** - 9 upregulated - 50 downregulated

**Dendritic cells** - 2 upregulated - 37 downregulated

**CD4 T cells** - 1 upregulated - 27 downregulated

**NK cells** - 3 upregulated - 24 downregulated

### Observation

The AChR-MG-associated transcriptional landscape was heterogeneous.
Neutrophils showed a prominent upward directional response, while B
cells, dendritic cells, CD4 T cells, and NK cells were predominantly
shifted in the opposite direction.

These results describe transcriptional direction and should not be
interpreted as direct functional activation or suppression.

------------------------------------------------------------------------

## Healthy vs Seronegative MG Pre-treatment

Neutrophils again showed the strongest upward profile:

-   39 upregulated
-   1 downregulated
-   40 significant DEG entries
-   97.5% upward directional bias

**B cells** - 30 upregulated - 107 downregulated

**Dendritic cells** - 20 upregulated - 80 downregulated

**NK cells** - 9 upregulated - 40 downregulated

**Monocytes** - 59 upregulated - 57 downregulated

### Observation

Pre-treatment seronegative MG showed a particularly strong
neutrophil-associated upward signature. B cells, dendritic cells, NK
cells, and several T-cell populations showed substantial downward
transcriptional changes.

Monocytes had a high DEG burden but an almost balanced up/down
distribution. This was interpreted as extensive transcriptional
remodeling rather than simple activation or suppression.

------------------------------------------------------------------------

## Seronegative MG Pre vs Post

The pre/post comparison showed a distinctly different directional
pattern.

**B cells** - 7 upregulated - 0 downregulated

**NK cells** - 17 upregulated - 5 downregulated

**Dendritic cells** - 13 upregulated - 6 downregulated

**Neutrophils** - 8 upregulated - 52 downregulated - 86.7% downward
directional bias

**CD8 T cells** - 1 upregulated - 20 downregulated

**T cells** - 3 upregulated - 20 downregulated

### Observation

The seronegative pre/post comparison showed a redistribution of
transcriptional direction. Neutrophils displayed a strong downward
profile, while B cells, NK cells, and dendritic cells showed
predominantly upward changes.

This cannot yet be described as "normalization." Demonstrating
normalization would require showing that the same genes altered before
treatment moved back toward the healthy-control expression state after
treatment.

------------------------------------------------------------------------

# Targeted AChR / MuSK Analysis

## Objective

Determine whether the immune-cell DEG dataset contained significant
transcript-level evidence involving AChR receptor-subunit genes or
`MUSK`.

## Results

Among the examined AChR receptor-subunit genes, only **CHRNB1** appeared
in the master DEG table.

The CHRNB1 observation:

-   occurred in progenitor cells
-   belonged to the Healthy vs AChR-MG comparison
-   average log2 fold-change ≈ −0.535
-   adjusted P value = 1.0
-   did not meet the final significance criteria

`CHRNA1`, `CHRND`, `CHRNE`, and `CHRNG` were not represented as
significant DEG evidence.

`MUSK` was not represented as a significant DEG.

### Final Target-Gene Summary

-   Significant AChR receptor-subunit DEG records: **0**
-   Significant `MUSK` DEG records: **0**

## Interpretation

The analyzed immune-cell DEG data do not support AChR- or MuSK-gene
overexpression.

An important conceptual distinction was established during this
analysis:

**AChR-positive MG refers to autoantibody status, not overexpression of
AChR genes.**

Likewise, antibody status involving MuSK should not be interpreted as
evidence that `MUSK` itself is transcriptionally overexpressed.

------------------------------------------------------------------------

# FINAL BIOLOGICAL INTERPRETATION

The central result of the project is that MG-associated transcriptional
alterations are **strongly cell-type specific**.

The immune system did not show a single uniform transcriptional
response.

Instead:

-   Neutrophils showed a prominent upward transcriptional signature in
    both AChR-associated MG and pre-treatment seronegative MG
    comparisons.
-   B cells showed substantial downward transcriptional changes in
    disease-versus-healthy comparisons.
-   Dendritic cells, NK cells, and T-cell populations displayed
    comparison-dependent directional patterns.
-   Monocytes could exhibit a large DEG burden without a strong net
    directional bias.
-   The seronegative pre/post comparison showed a different landscape,
    including strong downward changes in neutrophils and predominantly
    upward changes in B cells, NK cells, and dendritic cells.
-   DEG recurrence analysis suggested that many changes were restricted
    to relatively few cell types rather than being universally shared.
-   Direct AChR/MuSK transcript overexpression was not supported by the
    targeted DEG analysis.

These observations identify candidate immune-cell populations and
transcriptional programs for future pathway-level or experimental
investigation.

------------------------------------------------------------------------

# PROJECT LIMITATIONS

## Computational Limitations

The dataset contained more than 400,000 cells after QC and imposed
substantial memory and storage demands on desktop hardware.

These constraints affected:

-   Seurat layer joining
-   marker identification
-   large matrix operations
-   heatmap clustering
-   storage of intermediate Seurat objects

The workflow was adapted to remain computationally feasible.

## Biological Interpretation

DEG counts and DEG direction are not direct assays of cellular activity.

A predominance of upregulated genes does not automatically mean that a
cell is physiologically activated, and a predominance of downregulated
genes does not automatically mean that it is functionally suppressed.

Further interpretation would require pathway-level analysis, gene-set
scoring, protein-level evidence, functional assays, or independent
validation.

## Pre/Post Interpretation

The seronegative pre/post comparison demonstrates transcriptional change
between the two time points. The present analysis alone does not
establish that these changes represent treatment-induced normalization.

------------------------------------------------------------------------

# SKILLS AND CONCEPTS DEVELOPED

This project provided practical experience in:

-   GEO dataset exploration
-   public scRNA-seq data handling
-   R and Seurat
-   metadata construction
-   large-object management
-   quality control
-   normalization
-   variable-feature selection
-   scaling
-   PCA
-   UMAP
-   nearest-neighbor graph construction
-   Leiden clustering
-   Seurat v5 assay-layer handling
-   cell-type annotation
-   differential-expression analysis
-   DEG filtering and interpretation
-   cell-type-specific transcriptional analysis
-   recurrence analysis
-   conserved-signature analysis
-   scientific visualization
-   reproducible project organization
-   Git/GitHub workflow
-   computational troubleshooting
-   biological interpretation of transcriptomic data

------------------------------------------------------------------------

# MAJOR LESSONS FROM THE PROJECT

1.  Public datasets must be inspected before assumptions are made about
    their structure.
2.  Correct metadata is fundamental to every downstream biological
    comparison.
3.  Large single-cell datasets require resource-aware computational
    strategies.
4.  Seurat v5 assay layers must be handled deliberately.
5.  Strategic checkpoints prevent loss of computational work, but
    retaining every intermediate object is unnecessary.
6.  A technically successful figure is not automatically a useful
    scientific figure.
7.  DEG burden is not equivalent to cellular activation.
8.  Statistical significance, effect size, recurrence, cell type, and
    biological context should be interpreted together.
9.  Negative results can prevent incorrect biological claims.
10. Reproducibility includes code, metadata, directory organization,
    figures, tables, and documentation.

------------------------------------------------------------------------

# FINAL PROJECT STATUS --- 3 AUGUST 2026

## Analytical Work Completed

-   Dataset selection and acquisition
-   Processed-data inspection
-   Metadata construction
-   Individual Seurat object creation
-   Dataset merging
-   QC calculation and visualization
-   QC filtering
-   Normalization
-   Highly variable gene selection
-   Scaling
-   PCA
-   UMAP
-   nearest-neighbor graph construction
-   Leiden clustering
-   cell-type annotation
-   disease-comparison framework
-   cell-type-specific differential-expression analysis
-   integrated DEG analysis
-   global DEG visualization
-   shared and unique DEG analysis
-   gene-recurrence analysis
-   conserved-signature analysis
-   downstream reporting repair
-   figure refinement
-   directional transcriptional analysis
-   transcriptional perturbation ranking
-   targeted AChR/MuSK analysis
-   final biological interpretation

## Remaining Closure Tasks

The remaining work is project documentation and repository maintenance
rather than expansion of the biological analysis:

-   clean the GitHub repository
-   remove obsolete or duplicate files
-   standardize final script organization
-   finalize the root `README.md`
-   document the `DATA`, `SCRIPTS`, and `RESULTS` directories
-   ensure the repository contains only interpretable and reproducible
    outputs
-   make the final project commit

## Future Extensions

Possible future analyses include:

-   pathway enrichment of cell-type-specific upregulated and
    downregulated DEGs
-   gene-set or pathway activity scoring
-   explicit testing of whether disease-associated genes move toward
    healthy expression after the post-treatment time point
-   independent-dataset validation
-   deeper analysis of neutrophil, B-cell, dendritic-cell, NK-cell, and
    T-cell programs

These are **future extensions**, not unfinished requirements of the
present project.

------------------------------------------------------------------------

# CLOSING REFLECTION

This project began as an attempt to learn single-cell RNA-seq analysis
using a publicly available Myasthenia Gravis dataset. It developed into
a complete computational workflow involving more than 400,000 cells
after quality control, 34 initial clusters, multiple disease groups,
cell-type-specific differential-expression analysis, gene-recurrence
analysis, conserved-signature analysis, and targeted biological
interpretation.

A substantial part of the project involved troubleshooting rather than
simply executing a standard pipeline. Challenges included unfamiliar
public-data formats, corrupted checkpoint files, Seurat v5 assay-layer
behavior, memory limitations, storage constraints,
differential-expression failures, ambiguous joined columns, unsuitable
heatmaps, clustering errors, and downstream-reporting problems.

These difficulties became part of the learning outcome.

The most important development was conceptual. The project moved from
simply generating outputs toward understanding what each analytical step
measured, what conclusions were justified by the data, and where
interpretation had to stop.

The final analysis suggests that MG-associated transcriptional
alterations in this dataset are heterogeneous across immune-cell
populations. Neutrophils, B cells, dendritic cells, NK cells, and T-cell
populations show distinct directional patterns across disease and
pre/post comparisons. At the same time, the targeted analysis does not
support the simplistic interpretation that AChR-positive disease
corresponds to AChR-gene overexpression.

The present project is therefore considered complete at the level of
**cell-type-specific differential-expression characterization and
downstream descriptive interpretation**.

Further pathway-level or validation analyses may be undertaken as
independent future extensions, but they are not required to consider
this computational project closed.

**Analytical endpoint reached: 3 August 2026.**
