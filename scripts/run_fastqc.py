"""
Run FastQC on a FASTQ file.
"""
import argparse
import subprocess

def run_fastqc(input_fastq, output_dir):
    subprocess.run(["fastqc", input_fastq, "-o", output_dir], check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    run_fastqc(args.input, args.output)
