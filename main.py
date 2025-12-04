"""
Main script to run all pipelines in order.
"""
import subprocess

pipelines = [
    "snakemake -s workflows/rna_seq.smk --cores 4",
    "snakemake -s workflows/variant_calling.smk --cores 4",
    "snakemake -s workflows/barcode_counting.smk --cores 4",
    "snakemake -s workflows/multiomics_integration.smk --cores 4"
]

for cmd in pipelines:
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
