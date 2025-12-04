"""
Calculate enrichment scores from barcode counts.
"""
import argparse
import pandas as pd

def enrichment(input_csv, output_csv):
    df = pd.read_csv(input_csv, index_col=0)
    df['enrichment'] = df['count'] / df['count'].sum()
    df.to_csv(output_csv)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    enrichment(args.input, args.output)
