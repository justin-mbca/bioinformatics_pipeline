#!/usr/bin/env python3
"""Simple integrator: join annotated gene results with clinical summary on sample or cohort-level.
Usage: integrate_genomics_clinical.py genes.csv clinical_summary.csv combined_report.csv
"""
import sys
import pandas as pd


def main():
    genes = pd.read_csv(sys.argv[1])
    clinical = pd.read_csv(sys.argv[2])
    out = sys.argv[3]

    # Very simple integration: if genes contains a sample column, join; otherwise, produce a merged summary
    if 'sample' in genes.columns and 'sample' in clinical.columns:
        merged = genes.merge(clinical, on='sample', how='left')
    else:
        # Otherwise, attach clinical summary as metadata columns
        clinical_cols = {f"clinical_{c}":clinical[c].iloc[0] if c in clinical.columns else None for c in clinical.columns}
        for k,v in clinical_cols.items():
            genes[k] = v
        merged = genes

    merged.to_csv(out, index=False)
    print(f"Wrote integrated report to {out}")


if __name__ == '__main__':
    main()
