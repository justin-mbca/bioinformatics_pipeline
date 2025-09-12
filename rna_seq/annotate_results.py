"""
annotate_results.py
-------------------
Annotate DESeq2 results with gene information (e.g., gene symbols, descriptions).
Input: DESeq2 results CSV (with gene IDs)
Output: Annotated results CSV (with gene symbols, etc.)
"""

import pandas as pd
import sys

# Usage: python annotate_results.py deseq2_results.csv gene_annotation.csv annotated_results.csv

if __name__ == "__main__":
    # Use the global 'snakemake' variable if present (set by Snakemake), else fallback to CLI args
    if 'snakemake' in globals():
        de_file = snakemake.input[0]
        ann_file = snakemake.input[1]
        out_file = snakemake.output[0]
    else:
        if len(sys.argv) != 4:
            print("Usage: python annotate_results.py <deseq2_results.csv> <gene_annotation.csv> <annotated_results.csv>")
            sys.exit(1)
        de_file, ann_file, out_file = sys.argv[1:4]

    de_df = pd.read_csv(de_file, index_col=0)
    ann_df = pd.read_csv(ann_file)

    # Reset index to make gene IDs a column named 'gene_id'
    de_df = de_df.reset_index().rename(columns={de_df.index.name or 'index': 'gene_id'})

    # Merge on gene_id
    merged = pd.merge(de_df, ann_df, on='gene_id', how='left')
    merged.to_csv(out_file, index=False)
    print(f"Annotated results written to {out_file}")
