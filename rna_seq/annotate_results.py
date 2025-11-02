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

def annotate(de_file, ann_file, out_file, merge_key=None):
    """Annotate DE results.

    Args:
        de_file: path to DE results CSV (index may contain gene IDs)
        ann_file: path to annotation CSV
        out_file: path to write annotated CSV
        merge_key: optional explicit column name to merge on; if None, heuristics are used
    """

    # If an explicit merge_key is requested, read DE file without forcing the first
    # column to be the index so the key remains available as a column. Otherwise
    # read using the common pattern where gene IDs are in the index.
    if merge_key:
        de_df = pd.read_csv(de_file, index_col=None)
    else:
        de_df = pd.read_csv(de_file, index_col=0)
    ann_df = pd.read_csv(ann_file)

    # Normalize candidate key names and find a common merge key
    candidate_keys = ['gene_id', 'gene', 'Gene', 'gene_symbol', 'id']

    # Helper to ensure a dataframe has a 'gene_id' column
    def ensure_gene_id(df):
        # If already has gene_id, return
        if 'gene_id' in df.columns:
            return df
        # If index has a name or is meaningful, reset to column
        if df.index.name not in (None, ''):
            df = df.reset_index().rename(columns={df.index.name: 'gene_id'})
            return df
        # Otherwise if first column name is one of candidates, rename it
        first_col = df.columns[0]
        if first_col in candidate_keys:
            return df.rename(columns={first_col: 'gene_id'})
        # Otherwise, add an index-based gene_id column
        df = df.reset_index().rename(columns={df.index.name or 'index': 'gene_id'})
        return df

    de_df = ensure_gene_id(de_df)
    ann_df = ensure_gene_id(ann_df)

    # If explicit merge_key provided, prefer it (if present)
    if merge_key:
        if merge_key in de_df.columns and merge_key in ann_df.columns:
            # coerce merge columns to string to avoid object/int mismatches
            de_df[merge_key] = de_df[merge_key].astype(str)
            ann_df[merge_key] = ann_df[merge_key].astype(str)
            merged = pd.merge(de_df, ann_df, on=merge_key, how='left')
        else:
            # fallback to heuristics if explicit key isn't present in both
            merge_key = None

    if not merge_key:
        # If both have gene_id, merge; otherwise try to find any common column
        if 'gene_id' in de_df.columns and 'gene_id' in ann_df.columns:
            merged = pd.merge(de_df, ann_df, on='gene_id', how='left')
        else:
            common = set(de_df.columns).intersection(set(ann_df.columns))
            if common:
                key = list(common)[0]
                merged = pd.merge(de_df, ann_df, on=key, how='left')
            else:
                # Fallback: concatenate with annotation columns appended (no merge)
                merged = pd.concat([de_df.reset_index(drop=True), ann_df.reset_index(drop=True)], axis=1)

    merged.to_csv(out_file, index=False)
    print(f"Annotated results written to {out_file}")


if __name__ == "__main__":
    # Use the global 'snakemake' variable if present (set by Snakemake), else fallback to CLI args
    if 'snakemake' in globals():
        de_file = snakemake.input[0]
        ann_file = snakemake.input[1]
        out_file = snakemake.output[0]
        # allow config-driven merge key
        merge_key = None
        try:
            merge_key = snakemake.config.get('gene_annotation_key')
        except Exception:
            merge_key = None
    else:
        if len(sys.argv) not in (4, 5):
            print("Usage: python annotate_results.py <deseq2_results.csv> <gene_annotation.csv> <annotated_results.csv> [merge_key]")
            sys.exit(1)
        de_file, ann_file, out_file = sys.argv[1:4]
        merge_key = sys.argv[4] if len(sys.argv) == 5 else None

    annotate(de_file, ann_file, out_file, merge_key=merge_key)
