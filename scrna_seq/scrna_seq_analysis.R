# Single-cell RNA-Seq Analysis Pipeline using Seurat
# Author: Your Name
# Date: 2025-09-12

library(Seurat)
library(tidyverse)

# 1. Load data (10x Genomics format)
data_dir <- "data/demo_10x"  # Path to your demo 10x data
sc_data <- Read10X(data.dir = data_dir)
seurat_obj <- CreateSeuratObject(counts = sc_data, project = "scRNAseq")

# 2. Quality control
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 10 & nFeature_RNA < 500 & percent.mt < 100)

# 3. Normalization
seurat_obj <- NormalizeData(seurat_obj)

# 4. Feature selection
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)


# 5. Scaling and dimensionality reduction
pca_dims <- 1:10  # You can adjust this range for tuning
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs = max(pca_dims))
seurat_obj <- RunUMAP(seurat_obj, dims = pca_dims)

# 6. Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = pca_dims)
seurat_obj <- FindClusters(seurat_obj, resolution = 1.0)

# 7. Marker gene identification
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file = "results/markers.csv")

# 8. Visualization
png("results/umap_clusters.png")
DimPlot(seurat_obj, reduction = "umap", label = TRUE)
dev.off()

png("results/feature_plot.png")
FeaturePlot(seurat_obj, features = head(VariableFeatures(seurat_obj), 4))
dev.off()

# Save Seurat object
saveRDS(seurat_obj, file = "results/seurat_obj.rds")
