# Interactive Data Visualization

This folder contains interactive visualization tools for exploring RNA-Seq data:

## Python: Dash App
- **File:** `python/dash_app.py`
- **Purpose:** Interactive PCA plot for bulk RNA-Seq data.
- **Requirements:** dash, plotly, scikit-learn, pandas
- **How to run:**
	```bash
	cd visualization/python
	python dash_app.py
	```
	The app will launch in your browser at `http://127.0.0.1:8050/`.

## R: Shiny App
- **File:** `r/shiny_app.R`
- **Purpose:** Interactive UMAP and feature plots for single-cell RNA-Seq data (Seurat object).
- **Requirements:** shiny, Seurat, ggplot2
- **How to run:**
	Open R and run:
	```R
	shiny::runApp('r/shiny_app.R')
	```

---
Feel free to extend these apps for additional visualizations or data types.
