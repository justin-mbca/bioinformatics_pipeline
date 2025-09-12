# Single-cell RNA-Seq Analysis Pipeline (Seurat)

This pipeline performs single-cell RNA-Seq analysis using the Seurat R package. It is designed to be run separately from the bulk RNA-Seq pipeline.

## Steps
1. Data loading (10x or other formats)
2. Quality control and filtering
3. Normalization
4. Feature selection
5. Dimensionality reduction (PCA, UMAP)
6. Clustering
7. Marker gene identification
8. Visualization

## Usage
- Place your raw data (e.g., 10x output) in the `data/` folder.
- Run the analysis script: `Rscript scrna_seq_analysis.R`
- Output plots and tables will be saved in the `results/` folder.

## Requirements
- R (>= 4.0)
- Seurat
- tidyverse

See `scrna_seq_analysis.R` for details.
