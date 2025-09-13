# Unsupervised SCVI cell embedding for single-cell RNA-Seq data
# Requires: scvi-tools, scanpy, anndata

import scanpy as sc
import scvi
import anndata
import pandas as pd

# Load AnnData object (replace with your file path)
adata = sc.read_h5ad('../../scrna_seq/results/scrna_data.h5ad')

# Prepare AnnData for scvi-tools (no labels needed)
scvi.model.SCVI.setup_anndata(adata)

# Train SCVI model
vae = scvi.model.SCVI(adata)
vae.train()

# Get latent representation
adata.obsm['X_scvi'] = vae.get_latent_representation()

# Save results
adata.write('../../scrna_seq/results/scrna_data_scvi_annotated.h5ad')
print('SCVI unsupervised embedding complete. Results saved to scrna_data_scvi_annotated.h5ad')
