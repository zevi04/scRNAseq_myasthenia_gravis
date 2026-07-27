# 22 July

Today's objective
- Repository setup

Completed
- Created GitHub repository
- Locked dataset
- Defined research question

Questions
- What is the sample composition of GSE227835?

Tomorrow
- Download processed files
- Read associated publication


# 23 July 2026

## Objective

Set up the computational workflow for the single-cell RNA-seq analysis of Myasthenia Gravis (GSE227835) and inspect the processed data.

## Completed

- Created and published the GitHub repository.
- Defined the research question and locked the project scope.
- Created the project structure:
  - DATA/
  - DOCS/
  - FIGURES/
  - NOTEBOOKS/
  - RESULTS/
  - SCRIPTS/
- Added `PROJECT.md`, `README.md`, and `JOURNAL.md`.
- Created a `.gitignore` to exclude downloaded datasets from version control.
- Downloaded the processed GEO dataset (GSE227835).
- Extracted the supplementary archive into `DATA/PROCESSED/`.
- Examined the study design and identified the sample groups:
  - Healthy Controls
  - AChR-positive Myasthenia Gravis
  - Seronegative MG (Pre-treatment)
  - Seronegative MG (Post-treatment)

## Observations

- The processed dataset consists of multiple compressed `.txt.gz` files.
- The files are not immediately readable using `fread()` as initially expected.
- A preview of one file showed barcode-like entries rather than a standard gene × cell expression matrix.
- The processed data format needs to be understood before importing it into R.

## Challenges

- Assumed the processed files were simple expression matrices.
- Initial attempts to read the files in R were unsuccessful.
- Realized that understanding the file structure is the necessary next step before writing analysis code.

## Lessons Learned

- Never assume the format of a public dataset.
- Inspect the dataset structure before choosing import functions.
- Good computational biology begins with understanding the data, not writing code.

## Next Steps

- Determine the exact structure of the processed files.
- Identify what each `.txt.gz` file represents.
- Build the sample metadata table.
- Import the dataset into R using the appropriate workflow.


# 2026-07-26

## Progress

Today marked the completion of the data preparation phase of the scRNA-seq pipeline.

### Completed

- Successfully created **40 individual Seurat objects** from the processed GEO expression matrices.
- Merged all samples into a single Seurat object.

### Merged Dataset

- **Samples:** 40
- **Cells:** 444,357
- **Genes:** 36,601

The merged object was saved to: RESULTS/objects/merged_seurat_object.rds


### Quality Control

Calculated standard QC metrics for every cell:

- `nFeature_RNA`
- `nCount_RNA`
- `percent.mt`

QC summary:

| Metric | Min | Median | Max |
|---------|----:|-------:|----:|
| nFeature_RNA | 12 | 1474 | 8380 |
| nCount_RNA | 500 | 5168 | 107729 |
| percent.mt | 0.00 | 7.37 | 98.26 |

### QC Visualization

Generated and saved:

- QC Violin Plot
- Counts vs Features
- Counts vs Mitochondrial %
- Features vs Mitochondrial %

Location: FIGURES/qc/


Initial violin plots were cluttered because ~444k cells were plotted with individual points. Updated the visualization by removing point overlays (`pt.size = 0`), resulting in cleaner distributions.

### Observations

- Strong positive correlation between UMI counts and detected genes (~0.81).
- Most cells have mitochondrial percentages below ~10%.
- A subset of cells shows very high mitochondrial content, suggesting low-quality or dying cells that will be addressed during filtering.

## Next Steps

- Determine QC filtering thresholds.
- Remove low-quality cells.
- Normalize the filtered dataset.
- Begin downstream dimensionality reduction and clustering.



Work Completed
1. Quality Control Filtering

Applied filtering criteria to remove low-quality cells.

Filtering thresholds

nFeature_RNA > 200
nFeature_RNA < 6000
percent.mt < 15%

Results

Cells before filtering: 444,357
Cells after filtering: 420,538
Cells removed: 23,819
Percentage removed: 5.36%

The removal percentage indicates that the dataset was already of good quality while eliminating cells likely to introduce technical noise.

2. Data Normalization

Performed library-size normalization using Seurat's NormalizeData() function.

Method

LogNormalize
Scale factor = 10,000

Normalization adjusts gene expression values so that differences in sequencing depth between cells do not bias downstream analyses.

3. Highly Variable Gene Selection

Identified genes with the greatest biological variability using FindVariableFeatures().

Method

Selection: VST
Number of genes selected: 2,000

These genes will be used for dimensionality reduction because they contribute most to distinguishing different cell populations.

4. Data Scaling

Scaled the highly variable genes using ScaleData().

Scaling standardizes each gene to:

Mean = 0
Standard deviation = 1

This prevents highly expressed genes from dominating Principal Component Analysis.

Conceptual Understanding Gained

Today's session focused heavily on understanding the reasoning behind each preprocessing step instead of simply executing Seurat functions.

Key concepts learned
Quality control removes poor-quality cells, not genes.
Normalization corrects differences in sequencing depth between cells.
Highly variable genes capture the biological differences between cells.
Scaling standardizes expression values across genes before dimensionality reduction.
PCA does not classify or cluster cells; it compresses thousands of correlated gene expression measurements into a smaller number of informative principal components.
UMAP is a visualization method that projects cells into two dimensions while preserving neighborhood relationships.
Cell clustering occurs after PCA and UMAP and is based on similarity between cells.
Cell annotation is a biological interpretation step performed using known marker genes.
Project Status
Completed

✓ Dataset Exploration
✓ Data Import
✓ Metadata Construction
✓ Dataset Consistency Check
✓ Seurat Object Creation
✓ Merge Seurat Objects
✓ Quality Control Metrics
✓ QC Visualization
✓ QC Filtering
✓ Normalization
✓ Highly Variable Gene Selection
✓ Scaling
Repository Improvements

During today's work, the project structure was reorganized for improved readability.

Current organization:

RESULTS/
├── metadata/
├── objects/
├── qc/
├── tables/
├── figures/
└── logs/

This organization will make future analyses and figure management easier.

Challenges Encountered

Several technical issues were encountered and resolved.

Corrupted intermediate Seurat object caused by an incomplete saveRDS() operation.
Initial confusion regarding Seurat object storage and large file sizes.
Realization that intermediate Seurat objects occupy substantial disk space (approximately 12 GB total).
Discussed future optimization by retaining only essential checkpoint objects after the pipeline is finalized.
Reflections

An important realization today was that simply running scripts does not necessarily lead to understanding the analysis. Moving forward, the workflow will emphasize understanding the biological and computational rationale behind each step before implementing it in code.

This approach is intended to improve confidence in explaining the pipeline during interviews, research discussions, and future projects.

Next Steps

The next stage begins dimensionality reduction.

Planned analyses:

Principal Component Analysis (PCA)
Elbow Plot generation
PCA loading visualization
Selection of informative principal components
UMAP
Clustering
Cell type annotation
Personal Reflection

Today marked a shift from simply executing a bioinformatics pipeline to understanding how each analytical step contributes to interpreting single-cell RNA-seq data. The focus moved away from memorizing Seurat functions toward understanding the underlying concepts, with the goal of being able to explain the complete workflow confidently in a research interview or laboratory discussion.

