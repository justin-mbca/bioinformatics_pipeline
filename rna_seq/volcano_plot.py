"""
volcano_plot.py
--------------
Generate a volcano plot from annotated DESeq2 results.
Input: Annotated results CSV (with log2FoldChange, pvalue, gene_symbol columns)
Output: volcano_plot.png
"""

import pandas as pd
import matplotlib.pyplot as plt
import sys

# Usage: python volcano_plot.py annotated_results.csv volcano_plot.png

if __name__ == "__main__":
    import numpy as np
    # Use the global 'snakemake' variable if present (set by Snakemake), else fallback to CLI args
    if 'snakemake' in globals():
        in_file = snakemake.input[0]
        out_file = snakemake.output[0]
    else:
        if len(sys.argv) != 3:
            print("Usage: python volcano_plot.py <annotated_results.csv> <volcano_plot.png>")
            sys.exit(1)
        in_file, out_file = sys.argv[1:3]

    df = pd.read_csv(in_file)

    # Basic volcano plot: log2FoldChange vs -log10(pvalue)
    df['-log10(pvalue)'] = -np.log10(df['pvalue'])
    plt.figure(figsize=(8,6))
    plt.scatter(df['log2FoldChange'], df['-log10(pvalue)'], c='grey', alpha=0.7, s=10)
    plt.xlabel('log2 Fold Change')
    plt.ylabel('-log10(p-value)')
    plt.title('Volcano Plot')
    plt.tight_layout()
    plt.savefig(out_file)
    print(f"Volcano plot saved to {out_file}")
