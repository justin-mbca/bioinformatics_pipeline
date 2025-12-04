"""
Simulate barcode counting from FASTQ.
"""
import argparse
import pandas as pd

def count_barcodes(input_fastq, output_csv):
    # Placeholder: simulate random barcode counts
    barcodes = [f"BC{i:04d}" for i in range(100)]
    counts = pd.Series([abs(hash(bc)) % 1000 for bc in barcodes], index=barcodes)
    counts.to_csv(output_csv, header=["count"])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    count_barcodes(args.input, args.output)
