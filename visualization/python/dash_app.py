
# Dash app for interactive RNA-Seq visualization (PCA, Heatmap, Bar, Gene Expression)
import dash
from dash import dcc, html, Input, Output
import plotly.express as px
import plotly.figure_factory as ff
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

# Demo data loading
counts = pd.read_csv('../../rna_seq/sample_counts.csv', index_col=0)
gene_list = counts.index.tolist()
sample_list = counts.columns.tolist()

# Precompute PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(counts.T)
df_pca = pd.DataFrame(X_pca, columns=['PC1', 'PC2'])
df_pca['Sample'] = counts.columns

# Demo gene expression (pick first gene)
gene_expr = counts.iloc[0]
df_gene = pd.DataFrame({'Sample': sample_list, 'Expression': gene_expr.values})

app = dash.Dash(__name__)

app.layout = html.Div([
    html.H1('Bulk RNA-Seq Interactive Dashboard'),
    dcc.Dropdown(
        id='plot-type',
        options=[
            {'label': 'PCA', 'value': 'pca'},
            {'label': 'Heatmap', 'value': 'heatmap'},
            {'label': 'Bar Plot', 'value': 'bar'},
            {'label': 'Gene Expression', 'value': 'gene'}
        ],
        value='pca',
        clearable=False
    ),
    dcc.Dropdown(
        id='gene-dropdown',
        options=[{'label': g, 'value': g} for g in gene_list],
        value=gene_list[0],
        clearable=False,
        style={'marginTop': 10, 'display': 'none'}
    ),
    dcc.Graph(id='main-plot')
])

@app.callback(
    Output('main-plot', 'figure'),
    Output('gene-dropdown', 'style'),
    Input('plot-type', 'value'),
    Input('gene-dropdown', 'value')
)
def update_plot(plot_type, gene):
    if plot_type == 'pca':
        fig = px.scatter(df_pca, x='PC1', y='PC2', text='Sample', title='PCA of Bulk RNA-Seq Samples')
        return fig, {'marginTop': 10, 'display': 'none'}
    elif plot_type == 'heatmap':
        fig = px.imshow(counts, labels=dict(x='Sample', y='Gene', color='Count'),
                        aspect='auto', title='Expression Heatmap')
        return fig, {'marginTop': 10, 'display': 'none'}
    elif plot_type == 'bar':
        # Bar plot of total counts per sample
        total_counts = counts.sum(axis=0)
        fig = px.bar(x=sample_list, y=total_counts, labels={'x': 'Sample', 'y': 'Total Counts'},
                     title='Total Counts per Sample')
        return fig, {'marginTop': 10, 'display': 'none'}
    elif plot_type == 'gene':
        # Bar plot of expression for selected gene
        expr = counts.loc[gene]
        fig = px.bar(x=sample_list, y=expr, labels={'x': 'Sample', 'y': 'Expression'},
                     title=f'Expression of {gene} across Samples')
        return fig, {'marginTop': 10, 'display': 'block'}
    else:
        return {}, {'marginTop': 10, 'display': 'none'}

if __name__ == '__main__':
    app.run_server(debug=True)
