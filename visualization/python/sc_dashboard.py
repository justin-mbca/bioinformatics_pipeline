# Dash app for interactive single-cell RNA-Seq visualization with scANVI cell type annotation
import dash
from dash import dcc, html, Input, Output
import plotly.express as px
import scanpy as sc
import anndata
import pandas as pd
import numpy as np

# Load annotated AnnData (output from scanvi_annotate.py)
adata = sc.read_h5ad('../../scrna_seq/results/scrna_data_scanvi_annotated.h5ad')

# UMAP coordinates (assume already computed)
if 'X_umap' not in adata.obsm:
    sc.tl.umap(adata)
umap = adata.obsm['X_umap']

# Cell type labels from scANVI
cell_types = adata.obs['scanvi_labels'] if 'scanvi_labels' in adata.obs else None

# Genes and cells
gene_list = adata.var_names.tolist()
cell_list = adata.obs_names.tolist()

# Build Dash app
app = dash.Dash(__name__)

app.layout = html.Div([
    html.H1('Single-cell RNA-Seq Dashboard (scANVI Annotation)'),
    dcc.Dropdown(
        id='sc-plot-type',
        options=[
            {'label': 'UMAP', 'value': 'umap'},
            {'label': 'Gene Expression', 'value': 'gene'}
        ],
        value='umap',
        clearable=False
    ),
    dcc.Dropdown(
        id='sc-gene-dropdown',
        options=[{'label': g, 'value': g} for g in gene_list],
        value=gene_list[0],
        clearable=False,
        style={'marginTop': 10, 'display': 'none'}
    ),
    dcc.Graph(id='sc-main-plot')
])

@app.callback(
    Output('sc-main-plot', 'figure'),
    Output('sc-gene-dropdown', 'style'),
    Input('sc-plot-type', 'value'),
    Input('sc-gene-dropdown', 'value')
)
def update_sc_plot(plot_type, gene):
    if plot_type == 'umap':
        if cell_types is not None:
            fig = px.scatter(x=umap[:,0], y=umap[:,1], color=cell_types,
                             labels={'x': 'UMAP1', 'y': 'UMAP2', 'color': 'Cell Type'},
                             title='UMAP (scANVI Cell Type Annotation)')
        else:
            fig = px.scatter(x=umap[:,0], y=umap[:,1],
                             labels={'x': 'UMAP1', 'y': 'UMAP2'},
                             title='UMAP')
        return fig, {'marginTop': 10, 'display': 'none'}
    elif plot_type == 'gene':
        expr = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, 'toarray') else adata[:, gene].X.flatten()
        fig = px.scatter(x=umap[:,0], y=umap[:,1], color=expr,
                         labels={'x': 'UMAP1', 'y': 'UMAP2', 'color': f'{gene} Expression'},
                         title=f'UMAP colored by {gene} expression')
        return fig, {'marginTop': 10, 'display': 'block'}
    else:
        return {}, {'marginTop': 10, 'display': 'none'}

if __name__ == '__main__':
    app.run(debug=True)
