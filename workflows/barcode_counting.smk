# Barcode/library-based screen pipeline rules
BARCODE_SAMPLES = ["sample1"]
rule all:
    input:
        "results/barcodes/enrichment.csv"


rule barcode_count:
    input: "data/example/{sample}_barcodes.fastq.gz"
    output: "results/barcodes/{sample}_counts.csv"
    wildcard_constraints:
        sample="|".join(BARCODE_SAMPLES)
    shell: "python scripts/barcode_count.py --input {input} --output {output}"

rule enrichment:
    input: expand("results/barcodes/{sample}_counts.csv", sample=BARCODE_SAMPLES)
    output: "results/barcodes/enrichment.csv"
    shell: "python scripts/enrichment_analysis.py --input {input} --output {output}"
