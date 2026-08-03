# Documentation

This directory contains supporting documentation for the **Single-Cell RNA-seq Analysis of Myasthenia Gravis (GSE227835)** project.

The files in `DOCS/` provide dataset-specific reference information that complements the main project documentation without duplicating the full analysis, results, or chronological research record.

---

## Directory Contents

```text
DOCS/
├── Dataset_notes.md
└── README.md
```

### `Dataset_notes.md`

Detailed reference notes for the **GSE227835** dataset.

This document covers:

- dataset overview and biological context,
- study groups,
- GEO and GSM organization,
- interpretation of the processed `.txt.gz` files,
- relationship between GSM accessions and biological samples,
- data-import strategy,
- metadata construction,
- initial and post-QC dataset dimensions,
- downstream cell-type representation,
- terminology used for AChR-positive and seronegative MG groups,
- local `RAW/` and `PROCESSED/` data organization,
- dataset-specific reproducibility considerations and caveats.

It also records the answers to the three dataset questions that guided the initial data-inspection stage:

1. What does each `.txt.gz` file contain?
2. How should these files be imported?
3. What is the relationship between GSM files and biological samples?

See [`Dataset_notes.md`](Dataset_notes.md).

---

## Purpose of the `DOCS/` Directory

The `DOCS/` directory is intended for **supporting technical and dataset documentation**.

It is separate from the other project documents because each serves a different purpose:

| Document | Purpose |
| --- | --- |
| [`../README.md`](../README.md) | Public-facing overview of the complete project, major results, workflow, and repository |
| [`../PROJECT.md`](../PROJECT.md) | Scientific rationale, research questions, objectives, analytical design, scope, and limitations |
| [`../JOURNAL.md`](../JOURNAL.md) | Chronological record of project development, troubleshooting, decisions, and observations |
| [`Dataset_notes.md`](Dataset_notes.md) | Dataset-specific technical reference |
| `DOCS/README.md` | Guide to the documentation directory |

This separation keeps the root `README.md` concise enough to function as the repository landing page while preserving more detailed supporting information elsewhere.

---

## Dataset Reference

The project analyzes:

**GEO accession:** `GSE227835`  
**Data type:** Single-cell RNA sequencing  
**Biological material:** Human peripheral blood mononuclear cells (PBMCs)  
**Disease context:** Myasthenia Gravis

The project includes four biological groups:

- Healthy Controls
- AChR-positive MG
- Seronegative MG — Pre-treatment
- Seronegative MG — Post-treatment

Three principal comparisons are used downstream:

```text
Healthy_vs_AChR_MG
Healthy_vs_SNMG_Pre
SNMG_Pre_vs_SNMG_Post
```

Detailed dataset notes are maintained in [`Dataset_notes.md`](Dataset_notes.md).

---

## Related Repository Documentation

For analysis-specific information, use the corresponding repository documentation:

- [`../SCRIPTS/README.md`](../SCRIPTS/README.md) — script sequence and execution guidance
- [`../RESULTS/README.md`](../RESULTS/README.md) — organization of generated results
- [`../FIGURES/README.md`](../FIGURES/README.md) — organization of project figures
- [`../DATA/`](../DATA/) — local raw and processed data structure
- [`../PROJECT.md`](../PROJECT.md) — complete scientific project charter
- [`../JOURNAL.md`](../JOURNAL.md) — chronological research journal

---

## Documentation Principle

The repository documentation is organized so that information is recorded at the appropriate level:

```text
README.md
    │
    └── What is the project and what was found?

PROJECT.md
    │
    └── Why was the project performed and how was it designed?

JOURNAL.md
    │
    └── How did the analysis develop over time?

DOCS/
    │
    └── What technical and dataset-specific information is needed?

SCRIPTS/README.md
    │
    └── How is the computational workflow executed?

RESULTS/README.md
    │
    └── What outputs were generated and where are they stored?

FIGURES/README.md
    │
    └── What visual outputs are available and how are they organized?
```

This structure reduces duplication and makes the repository easier to navigate.

---

## Project Status

**Completed — 3 August 2026**

The documentation in this directory accompanies the completed cell-type-specific single-cell transcriptomic analysis of Myasthenia Gravis using GSE227835.
