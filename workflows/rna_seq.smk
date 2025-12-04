# RNA-seq pipeline rules
SAMPLES = ["sample1", "sample2"]
rule all:
    input:
        "results/rnaseq/deseq2_results.csv",
        "results/rnaseq/pca.pdf",
        "results/rnaseq/volcano.pdf"

rule fastqc:
    input: "data/example/{sample}.fastq.gz"
    output: "results/fastqc/{sample}_fastqc.html"
    shell: "fastqc {input} -o results/fastqc/"

rule multiqc:
    input: expand("results/fastqc/{sample}_fastqc.html", sample=SAMPLES)
    output: "results/fastqc/multiqc_report.html"
    shell: "multiqc results/fastqc/ -o results/fastqc/"

# Add STAR/HISAT2, featureCounts/Salmon, DESeq2, plotting rules here
