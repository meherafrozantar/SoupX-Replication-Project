# Replication Project: SoupX Ambient RNA Contamination Removal

This repository hosts the background correction replication code for processing droplet-based single-cell transcriptomic datasets using **SoupX**, evaluating the control metrics outlined by Matthew D. Young and Sam Behjati (2020).

## Methodology & Dataset
* **Replicated Control Target:** 10X Genomics PBMC 4k validation control sample dataset.
* **Lineage Exclusion Background Profiles:** Hemoglobin and Immunoglobin genes (`HBB`, `HBA1`, `HBA2`, `IGKC`, `IGHG1`, `IGHG2`, `IGHG3`).

## Project Replication Results
Running the automated pipeline isolates the background empty droplet profiles to calculate the final ambient mRNA technical noise value:
* **Replicated Contamination Fraction ($\rho$):** `4.21%`

## Script Execution Dependencies
* R (>= 4.0)
* SoupX 
* Seurat
