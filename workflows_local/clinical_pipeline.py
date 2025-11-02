#!/usr/bin/env python3
"""Local clinical pipeline script.
Usage: clinical_pipeline.py raw_clinical.csv cleaned_clinical.csv analysis_results.csv
"""
import sys
import pandas as pd


def main():
    raw = sys.argv[1]
    cleaned_out = sys.argv[2]
    analysis_out = sys.argv[3]

    df = pd.read_csv(raw)
    df_clean = df.dropna(how='all')
    df_clean.to_csv(cleaned_out, index=False)

    summary = df_clean.describe().transpose()
    summary.to_csv(analysis_out)
    print(f"Wrote cleaned clinical data to {cleaned_out} and analysis to {analysis_out}")


if __name__ == '__main__':
    main()
