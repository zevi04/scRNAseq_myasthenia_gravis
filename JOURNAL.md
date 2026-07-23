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