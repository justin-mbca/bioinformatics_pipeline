import scanpy as sc

adata = sc.read_h5ad('../../scrna_seq/results/scrna_data.h5ad')
print("Shape (cells x genes):", adata.shape)
print("X dtype:", adata.X.dtype)
print("X min/max:", adata.X.min(), adata.X.max())

# Handle sparse matrix for first 5 values
if hasattr(adata.X, 'data'):
    print("First 5 values in X:", adata.X.data[:5])
    # Check if all values are integers
    print("All integer values:", (adata.X.data % 1 == 0).all())
else:
    print("First 5 values in X:", adata.X.flatten()[:5])
    print("All integer values:", (adata.X % 1 == 0).all())

print("AnnData layers:", adata.layers.keys())
print("AnnData obs columns:", adata.obs.columns)
print("AnnData var columns:", adata.var.columns)