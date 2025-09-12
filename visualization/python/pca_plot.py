# PCA Plot in Python for bulk RNA-Seq
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

counts = pd.read_csv('../../rna_seq/sample_counts.csv', index_col=0)
pca = PCA(n_components=2)
X_pca = pca.fit_transform(counts.T)
plt.figure(figsize=(6,5))
plt.scatter(X_pca[:,0], X_pca[:,1], s=50)
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)')
plt.title('PCA of Bulk RNA-Seq Samples')
plt.tight_layout()
plt.show()
