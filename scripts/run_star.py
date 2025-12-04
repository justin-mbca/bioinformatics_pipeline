"""
Run STAR aligner on a FASTQ file.
"""
import argparse
import subprocess

def run_star(input_fastq, output_bam, star_index):
    subprocess.run([
        "STAR", "--genomeDir", star_index, "--readFilesIn", input_fastq,
        "--outSAMtype", "BAM", "SortedByCoordinate", "--outFileNamePrefix", output_bam.replace('.bam','')
    ], check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--index", required=True)
    args = parser.parse_args()
    run_star(args.input, args.output, args.index)
