# Multi-omics integration and ML pipeline rules
rule all:
    input:
        "results/integration/ml_metrics.csv"

rule integrate:
    input:
        rna="data/example/data_mrna_seq_v2_rsem.txt",
        variants="data/example/data_mutations.txt"
    output: "results/integration/integrated_features.csv"
    shell: "python scripts/multiomics_ml.py --rna {input.rna} --variants {input.variants} --output {output}"

rule ml:
    input:
        rna="data/example/data_mrna_seq_v2_rsem.txt",
        variants="data/example/data_mutations.txt",
        features="results/integration/integrated_features.csv"
    output: "results/integration/ml_metrics.csv"
    shell: "python scripts/multiomics_ml.py --rna {input.rna} --variants {input.variants} --output {input.features} --metrics {output} --run-ml"
