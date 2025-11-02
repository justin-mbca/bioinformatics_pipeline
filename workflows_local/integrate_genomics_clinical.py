#!/usr/bin/env python3
"""Local integrator script.
Usage: integrate_genomics_clinical.py genes.csv clinical_summary.csv combined_report.csv
"""
import sys
import pandas as pd


def main():
    genes = pd.read_csv(sys.argv[1])
    clinical = pd.read_csv(sys.argv[2])
    out = sys.argv[3]
    # Build standardized actionable report schema
    # Preferred columns: gene_id, symbol, log2FoldChange, padj, clinical_association, cohort_frequency, note
    report = pd.DataFrame()

    # gene_id - prefer human-readable name columns before numeric index-like columns
    gene_id_col = None
    for candidate in ['gene', 'Unnamed: 0', 'annotation', 'symbol', 'gene_id']:
        if candidate in genes.columns:
            gene_id_col = candidate
            break
    if gene_id_col is None:
        # pick the first non-numeric column if possible
        for c in genes.columns:
            if genes[c].dtype == object:
                gene_id_col = c
                break
    if gene_id_col is None:
        gene_id_col = genes.columns[0]

    report['gene_id'] = genes[gene_id_col].astype(str)

    # symbol - try genes, then annotation, then external mapping file
    if 'symbol' in genes.columns:
        report['symbol'] = genes['symbol']
    elif 'annotation' in genes.columns:
        report['symbol'] = genes['annotation']
    else:
        # attempt to look up gene symbol from local annotation file
        try:
            ann = pd.read_csv('workflows_local/gene_annotation.csv')
            # map columns 'gene' or 'gene_id' -> symbol/annotation
            key_col = None
            if 'gene' in ann.columns:
                key_col = 'gene'
            elif 'gene_id' in ann.columns:
                key_col = 'gene_id'
            if key_col:
                mapping = dict(zip(ann[key_col].astype(str), ann.get('annotation', ann.get('symbol', ann[key_col])).astype(str)))
                report['symbol'] = report['gene_id'].map(mapping).fillna(report['gene_id'])
            else:
                report['symbol'] = report['gene_id']
        except Exception:
            report['symbol'] = report['gene_id']

    # effect size and p-values
    report['log2FoldChange'] = genes.get('log2FoldChange', pd.NA)
    report['padj'] = genes.get('padj', pd.NA)

    # clinical association: try to pick a representative clinical column name
    clinical_assoc = None
    if clinical.shape[1] > 0:
        clinical_assoc = clinical.columns[0]
    report['clinical_association'] = clinical_assoc

    # cohort_frequency: if clinical provides a count column, use it; else placeholder
    if 'clinical_count' in genes.columns:
        # clinical_count might already be present after earlier merge
        report['cohort_frequency'] = genes['clinical_count']
    elif 'count' in clinical.columns:
        total = clinical['count'].sum() if clinical['count'].dtype.kind in 'biufc' else None
        if total and total > 0:
            report['cohort_frequency'] = clinical['count'] / total
        else:
            report['cohort_frequency'] = pd.NA
    else:
        report['cohort_frequency'] = pd.NA

    # note column for actionable recommendation (heuristics)
    report['note'] = pd.NA
    try:
        padj_vals = pd.to_numeric(report.get('padj', pd.Series([pd.NA]*len(report))), errors='coerce')
        l2fc_vals = pd.to_numeric(report.get('log2FoldChange', pd.Series([pd.NA]*len(report))), errors='coerce')
        notes = []
        for p, l in zip(padj_vals, l2fc_vals):
            if pd.notna(p) and p < 0.05:
                notes.append('validate by qPCR')
            elif pd.notna(l) and abs(l) > 1.0:
                notes.append('high priority')
            else:
                notes.append(pd.NA)
        report['note'] = notes
    except Exception:
        report['note'] = pd.NA

    report.to_csv(out, index=False)
    print(f"Wrote integrated report to {out}")


if __name__ == '__main__':
    main()
