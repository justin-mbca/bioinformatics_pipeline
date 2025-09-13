
# Interactive Data Visualization

This folder contains interactive dashboards for exploring RNA-Seq data:


## Python: Dash App (Bulk RNA-Seq)
- **File:** `python/dash_app.py`
- **Purpose:** Interactive dashboard for PCA, heatmap, bar, and gene expression plots.
- **Requirements:** dash, plotly, scikit-learn, pandas
- **How to run:**
  1. Open a terminal and navigate to the Python visualization folder:
	  ```bash
	  cd visualization/python
	  ```
  2. Run the Dash app:
	  ```bash
	  python dash_app.py
	  ```
  3. Open your browser and go to: [http://127.0.0.1:8050/](http://127.0.0.1:8050/)
  4. Use the dropdowns to explore all plot types interactively.

## Python: Dash App (Single-cell RNA-Seq, SCVI)
- **File:** `python/scvi_dashboard.py`
- **Purpose:** Interactive dashboard for UMAP and gene expression plots using SCVI unsupervised embedding.
- **Requirements:** dash, plotly, scanpy, scvi-tools, anndata
- **How to run:**
  1. Generate the SCVI embedding:
	  ```bash
	  python scanvi_annotate.py
	  ```
	  This will create `scrna_data_scvi_annotated.h5ad` in `scrna_seq/results/`.
  2. Run the Dash app:
	  ```bash
	  python scvi_dashboard.py
	  ```
  3. Open your browser and go to: [http://127.0.0.1:8050/](http://127.0.0.1:8050/)
  4. Use the dropdowns to explore UMAP and gene expression interactively.

## R: Shiny App (Single-cell RNA-Seq)
- **File:** `r/shiny_app.R`
- **Purpose:** Interactive dashboard for UMAP, violin, dot, bar, and gene expression plots (Seurat object).
- **Requirements:** shiny, Seurat, ggplot2
- **How to run:**
  1. Open R or RStudio.
  2. Set your working directory to the project root or use the full path.
  3. Run the app with:
	  ```R
	  shiny::runApp('visualization/r/shiny_app.R')
	  ```
  4. The app will open in your browser. Use the dropdowns to explore all plot types interactively.

---
**Note:**
- These dashboards are designed to run locally. GitHub does not support hosting or running interactive Dash or Shiny apps directly in the browser. To share interactive dashboards online, consider deploying to [Dash Enterprise](https://plotly.com/dash/app-manager/) or [ShinyApps.io](https://www.shinyapps.io/).

Feel free to extend these apps for additional visualizations or data types.
