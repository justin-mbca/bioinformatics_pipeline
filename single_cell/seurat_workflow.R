# Seurat Single-Cell RNA-Seq Analysis Example
# This script demonstrates a basic Seurat workflow for clustering and visualization.

library(Seurat)
library(ggplot2)

# Load example data (replace with your own data)
# pbmc_small is a built-in Seurat object for demonstration
sc_data <- pbmc_small

# Standard preprocessing
sc_data <- NormalizeData(sc_data)
sc_data <- FindVariableFeatures(sc_data)
sc_data <- ScaleData(sc_data)
sc_data <- RunPCA(sc_data)
sc_data <- FindNeighbors(sc_data, dims = 1:10)
sc_data <- FindClusters(sc_data, resolution = 0.5)
sc_data <- RunUMAP(sc_data, dims = 1:10)

# Visualization: UMAP plot colored by cluster
png('seurat_umap.png')
DimPlot(sc_data, reduction = 'umap', label = TRUE) + ggtitle('Seurat UMAP Clustering')
dev.off()

cat('Seurat single-cell analysis complete. UMAP plot saved as seurat_umap.png\n')
