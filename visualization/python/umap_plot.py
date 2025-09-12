# UMAP Plot in Python for single-cell RNA-Seq (Scanpy style)
import scanpy as sc
import matplotlib.pyplot as plt
adata = sc.read_10x_mtx('../../scrna_seq/data/demo_10x/', var_names='gene_symbols', cache=True)
sc.pp.filter_cells(adata, min_genes=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata, n_comps=10)
sc.pp.neighbors(adata, n_pcs=10)
sc.tl.umap(adata)
sc.pl.umap(adata, show=True)
