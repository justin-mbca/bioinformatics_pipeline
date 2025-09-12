# Example Snakemake workflow for RNA-Seq analysis
# Save as Snakefile in the workflows/ directory

configfile: "config.yaml"

rule all:
    input:
        "../rna_seq/annotated_results.csv",
        "../rna_seq/volcano_plot.png"

rule preprocess_qc:
    input:
        counts="../rna_seq/sample_counts.csv"
    output:
        cleaned="../rna_seq/counts_for_deseq2.csv",
        plot="../rna_seq/library_sizes.png"
    script:
        "../rna_seq/preprocess_qc.py"

rule deseq2:
    input:
        counts="../rna_seq/counts_for_deseq2.csv"
    output:
        results="../rna_seq/deseq2_results.csv",
        maplot="../rna_seq/deseq2_MAplot.png"
    script:
        "../rna_seq/deseq2_analysis.R"

rule annotate:
    input:
        results="../rna_seq/deseq2_results.csv",
        annotation=config["gene_annotation"]
    output:
        annotated="../rna_seq/annotated_results.csv"
    script:
        "../rna_seq/annotate_results.py"

rule volcano_plot:
    input:
        annotated="../rna_seq/annotated_results.csv"
    output:
        plot="../rna_seq/volcano_plot.png"
    script:
        "../rna_seq/volcano_plot.py"
