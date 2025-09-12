# Basic Dash app for interactive PCA plot (bulk RNA-Seq)
import dash
from dash import dcc, html
import plotly.express as px
import pandas as pd
from sklearn.decomposition import PCA

counts = pd.read_csv('../../rna_seq/sample_counts.csv', index_col=0)
pca = PCA(n_components=2)
X_pca = pca.fit_transform(counts.T)
df = pd.DataFrame(X_pca, columns=['PC1', 'PC2'])
df['Sample'] = counts.columns

app = dash.Dash(__name__)

app.layout = html.Div([
    html.H1('Bulk RNA-Seq PCA Explorer'),
    dcc.Graph(
        id='pca-plot',
        figure=px.scatter(df, x='PC1', y='PC2', text='Sample', title='PCA of Bulk RNA-Seq Samples')
    )
])

if __name__ == '__main__':
    app.run_server(debug=True)
