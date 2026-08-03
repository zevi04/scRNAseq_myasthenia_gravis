# GSE227835 — Dataset Notes

## Dataset Overview

**GEO accession:** `GSE227835`  
**Disease:** Myasthenia Gravis (MG)  
**Technology:** Single-cell RNA sequencing (scRNA-seq)  
**Biological material:** Peripheral blood mononuclear cells (PBMCs)  
**Organism:** Human

This dataset was selected for an independent single-cell transcriptomic reanalysis of immune-cell populations in Myasthenia Gravis.

The project uses the processed expression files associated with **GSE227835** to construct individual Seurat objects, combine samples into a unified dataset, perform quality control and preprocessing, identify transcriptional clusters, annotate cell populations, and carry out cell-type-specific downstream analyses.

---

## Biological Context

Myasthenia Gravis is an autoimmune disorder affecting neuromuscular transmission. The dataset contains samples representing healthy controls and different MG groups, allowing immune-cell transcriptional profiles to be compared across disease states.

Because the analyzed material is peripheral blood, the project focuses on the **circulating immune-cell transcriptional landscape** represented in the dataset.

The analysis does not treat AChR-positive status as evidence of AChR-gene overexpression. AChR-positive MG refers to the relevant autoantibody-defined disease group.

---

## Study Groups

Four biological groups were represented in the project metadata:

| Group | Project label / interpretation |
| --- | --- |
| Healthy Controls | Reference group |
| AChR-positive MG | AChR autoantibody-positive MG |
| Seronegative MG — Pre-treatment | Seronegative MG at the pre-treatment time point |
| Seronegative MG — Post-treatment | Seronegative MG at the post-treatment time point |

These groups formed the basis for the downstream biological comparisons.

---

## Biological Comparisons Used

Three principal comparisons were evaluated:

| Comparison | Purpose |
| --- | --- |
| `Healthy_vs_AChR_MG` | Compare AChR-positive MG with healthy controls |
| `Healthy_vs_SNMG_Pre` | Compare pre-treatment seronegative MG with healthy controls |
| `SNMG_Pre_vs_SNMG_Post` | Compare seronegative MG pre- and post-treatment states |

Differential-expression analysis was subsequently interpreted within individual annotated cell populations.

---

## GEO File Organization

The processed GEO download consisted of compressed text files with names containing **GSM accession identifiers**.

A GSM identifier is a GEO **sample accession**. It identifies a deposited sample record within the larger GSE study.

The project therefore treated the GSM-associated processed files as sample-level inputs rather than assuming that each filename represented a separate cell type or biological condition.

Biological group information had to be represented separately in the project metadata so that each imported sample could later be assigned to the appropriate experimental condition.

---

## What Do the `.txt.gz` Files Contain?

The `.txt.gz` files are **gzip-compressed text files** supplied as processed data associated with individual GEO sample records.

During initial inspection of the dataset, these files were found to contain expression information that could be read into R and converted into matrices suitable for construction of Seurat objects.

An important observation from the early inspection was that the processed files did **not necessarily resemble a standard 10x Genomics directory** containing separate:

```text
barcodes.tsv.gz
features.tsv.gz
matrix.mtx.gz
```

files.

For example, inspection of files such as:

```text
GSM7266236_A1.txt.gz
```

showed that the downloaded processed-data organization had to be examined directly before choosing an import strategy.

The project therefore did not assume a standard `Read10X()` input structure. Instead, the files were inspected and imported according to their actual tabular organization.

---

## How Were the Files Imported?

The import workflow was developed in stages rather than attempting to load the complete dataset immediately.

The general strategy was:

```text
Processed GEO .txt.gz file
        ↓
Inspect file structure
        ↓
Read tabular expression data into R
        ↓
Verify dimensions and identifiers
        ↓
Convert expression values into a matrix
        ↓
Construct an individual Seurat object
        ↓
Attach sample-level metadata
        ↓
Repeat for all samples
        ↓
Merge individual Seurat objects
```

The early scripts were specifically used to determine how the processed files were structured before object construction.

Broadly, the relevant script stages were:

| Script stage | Purpose |
| --- | --- |
| `00` | Inspect the processed data files |
| `01` | Import data |
| `02` | Build sample metadata |
| `03` | Check dataset consistency |
| `04` | Create a Seurat object |
| `05` | Create all sample-level Seurat objects |
| `06` | Merge the Seurat objects |

For the exact current filenames and execution order, refer to:

[`../SCRIPTS/README.md`](../SCRIPTS/README.md)

---

## GSM Accessions and Biological Samples

### What is a GSM?

In the NCBI Gene Expression Omnibus hierarchy:

- **GSE** identifies the overall study or series.
- **GSM** identifies an individual GEO sample record associated with that study.

Therefore:

```text
GSE227835
    │
    ├── GSM...
    ├── GSM...
    ├── GSM...
    └── ...
```

The GSM identifier is useful for tracking the origin of each processed input file.

### How was this handled in the project?

Each imported sample was linked to project metadata describing its biological group.

This separation was important because the expression matrix provides the molecular measurements, whereas the metadata provides the biological interpretation required for comparisons such as:

```text
Healthy
AChR-positive MG
Seronegative MG Pre
Seronegative MG Post
```

The GSM/sample identity was retained so that cells could remain traceable to their source sample after the individual Seurat objects were merged.

---

## Metadata Construction

Sample metadata was constructed separately from the expression matrices.

The metadata provided the information needed to associate each imported sample with its biological condition and to preserve sample identity during downstream analysis.

Conceptually:

```text
GSM / sample identifier
        +
Biological group
        +
Imported expression matrix
        ↓
Sample-level Seurat object
```

After all sample-level objects were created, they were merged into the complete project object while retaining the relevant sample and condition information.

This metadata structure later enabled condition-specific and cell-type-specific comparisons.

---

## Initial Dataset Dimensions

Following sample-level object construction and merging, the combined dataset contained:

| Metric | Value |
| --- | ---: |
| Individual Seurat objects | **40** |
| Cells | **444,357** |
| Genes | **36,601** |

These values represent the merged dataset before the final QC filtering step.

---

## Quality-Control Dataset

Quality control used:

- `nFeature_RNA`
- `nCount_RNA`
- `percent.mt`

The filtering criteria were:

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

The **420,538 retained cells** formed the basis for subsequent preprocessing and downstream analysis.

---

## Downstream Cell-Type Representation

Unsupervised analysis initially identified **34 transcriptional clusters**.

These were subsequently interpreted and consolidated into **10 downstream cell-type labels**:

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

These labels were used for cell-type-specific downstream comparisons.

---

## Important Dataset Terminology

### AChR-positive MG

In this project, **AChR-positive MG** refers to the disease group characterized by acetylcholine-receptor autoantibody positivity.

It does **not** mean that genes encoding AChR receptor subunits are necessarily overexpressed in the scRNA-seq data.

### Seronegative MG

**Seronegative MG** refers to the seronegative disease group represented in the dataset/project metadata.

The label should be interpreted according to the study grouping rather than as a statement about transcriptional expression of any single target gene.

### Pre-treatment

`Pre` refers to the pre-treatment time point represented in the seronegative MG group.

### Post-treatment

`Post` refers to the corresponding post-treatment time point represented in the project metadata.

The `SNMG_Pre_vs_SNMG_Post` comparison therefore measures transcriptional differences between these two sampled states.

A difference between pre and post should not automatically be described as treatment-induced normalization without demonstrating movement toward the healthy-control expression state.

---

## Local Data Organization

The repository uses:

```text
DATA/
├── RAW/
└── PROCESSED/
```

### `DATA/RAW/`

Intended for the original downloaded dataset files retained locally.

### `DATA/PROCESSED/`

Intended for processed or analysis-ready data generated or organized during the workflow.

Large data files are excluded from GitHub through `.gitignore`.

The public repository therefore documents the data structure and analysis workflow without redistributing the complete expression dataset.

---

## Why the Data Are Not Stored in GitHub

Single-cell expression datasets and serialized Seurat objects can be very large.

To keep the repository manageable:

- raw data are maintained locally,
- large processed data are maintained locally,
- large `.rds` objects are excluded,
- analysis scripts and documentation are version controlled.

The public repository is therefore intended to preserve the **analysis logic, documentation, and selected outputs**, rather than function as a duplicate host for the GEO dataset.

---

## Resolved Initial Questions

The project originally began with three dataset-level questions.

### 1. What does each `.txt.gz` file contain?

The files are gzip-compressed processed text data associated with GEO GSM sample records. Their internal tabular structure was inspected directly before import rather than assuming a standard 10x directory layout.

### 2. How should these files be imported?

They were imported through a staged R workflow: inspect the file structure, read the tabular expression data, verify the resulting matrix, construct sample-level Seurat objects, attach metadata, and merge the individual objects.

### 3. What is the relationship between GSM files and biological samples?

`GSM` is the GEO sample accession. GSM-associated files provide sample-level processed data, while project metadata links those sample identities to biological groups such as healthy control, AChR-positive MG, seronegative MG pre-treatment, and seronegative MG post-treatment.

These questions are retained here because resolving them was necessary before the downstream single-cell workflow could be constructed correctly.

---

## Dataset-Specific Caveats

### Processed rather than assumed raw 10x structure

The downloaded GEO processed files required inspection before import. The workflow should not be reproduced by assuming that the dataset is organized as a conventional `matrix.mtx.gz` / `features.tsv.gz` / `barcodes.tsv.gz` directory unless the files being used actually follow that format.

### Metadata are essential

Expression matrices alone do not encode all biological comparisons. Correct sample-to-condition mapping is required before disease-group comparisons can be interpreted.

### Sample identity should be preserved

Cells from different samples are merged for downstream analysis, but their original sample identities should remain available in metadata.

### Large object size

With more than 400,000 QC-retained cells, the complete Seurat object can require substantial memory and storage. Intermediate `.rds` checkpoints may therefore be useful locally even though they are excluded from GitHub.

### Interpretation is limited by the source dataset

The analysis can only evaluate biological and clinical variables represented in the available dataset and associated metadata.

---

## Related Documentation

- [`../README.md`](../README.md) — complete project overview and major results
- [`../PROJECT.md`](../PROJECT.md) — scientific rationale, objectives, scope, and analytical design
- [`../JOURNAL.md`](../JOURNAL.md) — chronological analysis and troubleshooting record
- [`README.md`](README.md) — documentation-directory guide
- [`../SCRIPTS/README.md`](../SCRIPTS/README.md) — script execution order
- [`../RESULTS/README.md`](../RESULTS/README.md) — generated-output guide

---

## Data Source

**Dataset:** GSE227835  
**Repository:** NCBI Gene Expression Omnibus (GEO)  
**Data type:** Single-cell RNA sequencing  
**Biological material:** Human peripheral blood mononuclear cells  
**Disease context:** Myasthenia Gravis

The original dataset remains available through its public GEO record. This repository contains an independent secondary computational analysis and does not redistribute the original complete dataset.
