# CSE763 Assignment: SoupX Ambient RNA Replication Workflow

This repository contains the functional implementation and technical replication report for identifying, quantifying, and subtracting cell-free extracellular ambient mRNA technical contamination from droplet-based single-cell RNA sequencing (scRNA-seq) matrices. The architecture follows the statistical estimation models engineered by Matthew D. Young and Sam Behjati (2020).

## Academic Profile
* **Student Name:** Meher Afroz
* **Student ID:** 1000058897
* **Course Code:** CSE763 (Advanced Bioinformatics)
* **Institution:** BRAC University, Dhaka, Bangladesh

## Core Objective
To evaluate background transcript expression noise across microfluidic droplet structures by profiling cell-free empty partition vectors, constructing maximum likelihood diagnostic curves, and outputting an adjusted gene-expression matrix.

## Dataset & Technical Controls
* **Dataset Benchmarked:** 10X Genomics Peripheral Blood Mononuclear Cell (PBMC) validation control dataset.
* **Lineage Exclusion Marker Panels:** Hemoglobin and Immunoglobin sequence vectors (`HBB`, `HBA1`, `HBA2`, `IGKC`, `IGHG1`, `IGHG2`, `IGHG3`).

## Project Replication Metrics
The analysis pipeline evaluates trace mRNA background distribution and converges on an explicit calculated metric matching the validation curves:
* **Replicated Contamination Fraction ($\rho$):** `6.2%` (Posterior Estimation Interval Range: 1.3% to 14%)

## Repository Manifest
* `soupX_script.R`: The complete automated single-cell pre-processing execution script featuring Seurat graph-clustering layers and live diagnostic density plotting functions.

## Implementation Steps
To run this pipeline inside an R/RStudio environment:
1. Open the project script file `soupX_script.R`.
2. Ensure the required single-cell toolkit extensions are pre-installed via the R console:
   ```R
   install.packages(c("SoupX", "Seurat"))
   ```
3. Source the script to execute the integrated validation matrix array and output the diagnostic curve:
   ```R
   source("soupX_script.R")
   ```
