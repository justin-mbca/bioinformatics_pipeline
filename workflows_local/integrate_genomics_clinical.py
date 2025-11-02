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

    # gene_id
    if 'gene_id' in genes.columns:
        report['gene_id'] = genes['gene_id']
    elif 'Unnamed: 0' in genes.columns:
        report['gene_id'] = genes['Unnamed: 0']
    elif 'annotation' in genes.columns:
        report['gene_id'] = genes['annotation']
    else:
        # fallback to first column
        report['gene_id'] = genes.iloc[:, 0]

    # symbol
    if 'symbol' in genes.columns:
        report['symbol'] = genes['symbol']
    elif 'annotation' in genes.columns:
        report['symbol'] = genes['annotation']
    else:
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

    # note column for actionable recommendation (left blank for now)
    report['note'] = pd.NA

    report.to_csv(out, index=False)
    print(f"Wrote integrated report to {out}")


if __name__ == '__main__':
    main()
