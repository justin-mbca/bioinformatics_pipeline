# Dash app for interactive single-cell RNA-Seq visualization using SCVI embedding
import dash
from dash import dcc, html, Input, Output
import plotly.express as px
import scanpy as sc
import numpy as np

# Load SCVI-annotated AnnData
adata = sc.read_h5ad('../../scrna_seq/results/scrna_data_scvi_annotated.h5ad')

# Get SCVI latent embedding
scvi_latent = adata.obsm['X_scvi']

# Compute UMAP on SCVI embedding if not present
if 'X_umap' not in adata.obsm:
    import scanpy as sc
    sc.pp.neighbors(adata, use_rep='X_scvi')
    sc.tl.umap(adata)
umap = adata.obsm['X_umap']

# Genes and cells
gene_list = adata.var_names.tolist()
cell_list = adata.obs_names.tolist()

app = dash.Dash(__name__)

app.layout = html.Div([
    html.H1('Single-cell RNA-Seq Dashboard (SCVI Embedding)'),
    dcc.Dropdown(
        id='scvi-plot-type',
        options=[
            {'label': 'UMAP (SCVI)', 'value': 'umap'},
            {'label': 'Gene Expression', 'value': 'gene'}
        ],
        value='umap',
        clearable=False
    ),
    dcc.Dropdown(
        id='scvi-gene-dropdown',
        options=[{'label': g, 'value': g} for g in gene_list],
        value=gene_list[0],
        clearable=False,
        style={'marginTop': 10, 'display': 'none'}
    ),
    dcc.Graph(id='scvi-main-plot')
])

@app.callback(
    Output('scvi-main-plot', 'figure'),
    Output('scvi-gene-dropdown', 'style'),
    Input('scvi-plot-type', 'value'),
    Input('scvi-gene-dropdown', 'value')
)
def update_scvi_plot(plot_type, gene):
    if plot_type == 'umap':
        fig = px.scatter(x=umap[:,0], y=umap[:,1],
                         labels={'x': 'UMAP1', 'y': 'UMAP2'},
                         title='UMAP (SCVI Embedding)')
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
