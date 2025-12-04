"""
Run featureCounts for quantification.
"""
import argparse
import subprocess

def run_featurecounts(input_bam, gtf, output_txt):
    subprocess.run([
        "featureCounts", "-a", gtf, "-o", output_txt, input_bam
    ], check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    run_featurecounts(args.input, args.gtf, args.output)
