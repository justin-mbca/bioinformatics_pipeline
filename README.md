# Bioinformatics Pipeline

This project provides a modular framework for developing and optimizing bioinformatics tools and pipelines, with a focus on RNA-Seq and single-cell RNA-Seq analysis, cohort analysis, and interactive data visualization. It is designed for collaboration with laboratory scientists and supports drug discovery initiatives.

## Features
- Bulk RNA-Seq analysis (DESeq2, cohort analysis)
- Single-cell RNA-Seq analysis (Seurat, R)
- Interactive data visualization (Plotly Dash, Shiny)
- Workflow management (Snakemake, Nextflow)
- Protocols, SOPs, and QC reporting
- Training and documentation for junior staff

## Project Structure
- `rna_seq/` — Bulk RNA-Seq analysis scripts and pipelines
- `single_cell/` — Single-cell RNA-Seq analysis scripts and pipelines
- `visualization/` — Interactive visualization tools (Dash, Shiny)
- `workflows/` — Workflow management scripts (Snakemake, Nextflow)
- `docs/` — Protocols, SOPs, and documentation

## Getting Started
1. Clone the repository and set up your Python and R environments.
2. Explore the `rna_seq/` and `single_cell/` folders for analysis scripts.
3. Use the `workflows/` folder for reproducible pipeline execution.
4. See `docs/` for protocols and training materials.

## Collaboration
Contributions are welcome! Please open issues or submit pull requests for new features, bug fixes, or documentation improvements.

## License
MIT License.

## Bulk RNA-Seq Analysis Pipeline (Step-by-Step)

This pipeline automates bulk RNA-Seq analysis from raw count matrix to annotated results and visualization. It is managed by Snakemake for reproducibility and modularity.

### Steps:

1. **Preprocessing & QC** (`preprocess_qc.py`)
   - Input: Raw count matrix (`sample_counts.csv`)
   - Performs quality control, summary statistics, and generates a library size plot.
   - Outputs: Cleaned count matrix (`counts_for_deseq2.csv`), QC plot (`library_sizes.png`).

2. **Differential Expression Analysis** (`deseq2_analysis.R`)
   - Input: Cleaned count matrix
   - Runs DESeq2 in R to identify differentially expressed genes between conditions.
   - Outputs: DESeq2 results (`deseq2_results.csv`), MA plot (`deseq2_MAplot.png`).

3. **Annotation** (`annotate_results.py`)
   - Input: DESeq2 results and gene annotation file (`gene_annotation.csv`)
   - Merges gene-level results with gene symbols/descriptions for interpretability.
   - Output: Annotated results (`annotated_results.csv`).

4. **Visualization** (`volcano_plot.py`)
   - Input: Annotated results
   - Generates a volcano plot to visualize significant genes.
   - Output: Volcano plot image (`volcano_plot.png`).


### How to Run

1. Activate the virtual environment:
   ```bash
   source .venv/bin/activate
   ```
2. Run the pipeline from the project root directory:
   ```bash
   snakemake --snakefile workflows/Snakefile --cores 1 --printshellcmds
   ```
3. All outputs will be generated in the `rna_seq/` directory.

### File Overview
- `rna_seq/sample_counts.csv`: Example input count matrix
- `rna_seq/gene_annotation.csv`: Example gene annotation file
- `rna_seq/preprocess_qc.py`: Python script for QC and preprocessing
- `rna_seq/deseq2_analysis.R`: R script for DESeq2 analysis
- `rna_seq/annotate_results.py`: Python script for annotation
- `rna_seq/volcano_plot.py`: Python script for volcano plot
- `workflows/Snakefile`: Snakemake workflow definition

### Customization
- Edit `config.yaml` to change input files or parameters.
- Replace example data with your own for real analyses.

For questions or troubleshooting, see the comments in each script or open an issue.

## Data Processing & Analysis Pipeline (Mermaid Diagram)

```mermaid
graph TD
    A[Raw Count Matrix (sample_counts.csv)] --> B[Preprocessing & QC (preprocess_qc.py)]
    B --> C[Cleaned Counts (counts_for_deseq2.csv)]
    C --> D[Differential Expression (deseq2_analysis.R)]
    D --> E[DESeq2 Results (deseq2_results.csv)]
    E --> F[Annotation (annotate_results.py)]
    F --> G[Annotated Results (annotated_results.csv)]
    G --> H[Visualization (volcano_plot.py)]
    H --> I[Volcano Plot (volcano_plot.png)]
    E --> J[MA Plot (deseq2_MAplot.png)]
    B --> K[QC Plot (library_sizes.png)]
```

This diagram summarizes the flow of data and analysis steps in the bulk RNA-Seq pipeline.

## Future Expansion

This repository is designed for extensibility. Planned and possible future enhancements include:

- **Single-cell RNA-Seq Analysis:** Integrate Seurat-based workflows for single-cell data processing and visualization.
- **Multi-omics Integration:** Add support for proteomics, metabolomics, or other omics data types.
- **Advanced Visualization:** Expand interactive dashboards (Dash, Shiny) for deeper data exploration.
- **Cloud & HPC Support:** Enable execution on cloud platforms or high-performance computing clusters.
- **Parameter Sweeps & Batch Analysis:** Automate large-scale analyses with parameter sweeps and batch processing.
- **Automated Reporting:** Generate comprehensive HTML or PDF reports summarizing results and QC.
- **User-friendly Interfaces:** Add GUIs or web interfaces for non-technical users.

Contributions and suggestions for new features are welcome!
