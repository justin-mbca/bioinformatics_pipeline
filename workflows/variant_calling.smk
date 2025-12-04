# Variant calling pipeline rules
rule all:
    input:
        "results/variants/filtered_variants.csv"

rule bwa_mem:
    input: "data/example/{sample}.fastq.gz"
    output: "results/variants/{sample}.bam"
    shell: "bwa mem ref.fa {input} | samtools sort -o {output}"

rule gatk_haplotypecaller:
    input: "results/variants/{sample}.bam"
    output: "results/variants/{sample}.vcf"
    shell: "gatk HaplotypeCaller -I {input} -O {output} -R ref.fa"

# Add VEP/ANNOVAR annotation, filtering, and prioritization rules here
