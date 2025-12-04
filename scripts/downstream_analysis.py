# Downstream analysis and visualization for integrated multi-omics features
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import seaborn as sns

# Load integrated features
df = pd.read_csv('results/integration/integrated_features.csv', index_col=0)

# Separate RNA and mutation features (assume RNA first 50 columns, rest are mutation)
rna = df.iloc[:, :50]
mut = df.iloc[:, 50:]

# --- Summary statistics ---
print('RNA-seq feature summary:')
print(rna.describe())
print('\nMutation feature summary:')
print(mut.sum(axis=0).sort_values(ascending=False).head(10))

# --- PCA on RNA features ---
pca = PCA(n_components=2)
X_pca = pca.fit_transform(rna)
plt.figure(figsize=(6,5))
plt.scatter(X_pca[:,0], X_pca[:,1], s=40)
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)')
plt.title('PCA of RNA-seq Features')
plt.tight_layout()
plt.savefig('results/integration/pca_rna.png')
plt.close()

# --- Mutation burden distribution ---
mut_burden = mut.sum(axis=1)
sns.histplot(mut_burden, bins=30, kde=True)
plt.xlabel('Mutation Burden (total mutated genes)')
plt.title('Mutation Burden Distribution')
plt.tight_layout()
plt.savefig('results/integration/mutation_burden.png')
plt.close()

print('Analysis and plots saved to results/integration/')
