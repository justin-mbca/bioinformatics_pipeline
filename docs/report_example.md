# Example Integrated Actionable Report

This example shows a small, human-readable slice of the integrated biomarker report produced by the pipeline. It demonstrates the lightweight schema used to communicate prioritized findings.

| gene_id | symbol | log2FoldChange | padj | clinical_association | cohort_frequency | note |
|---|---:|---:|---:|---|---:|---|
| GeneA | GeneA | -0.0220 | 0.9936 | AGE | 0.02 | validate by qPCR |
| GeneB | GeneB | -0.0878 | 0.9936 | AGE | 0.01 | low priority |
| GeneC | GeneC | 0.0944 | 0.9936 | AGE | 0.03 | validate by qPCR |
| GeneD | GeneD | 0.0190 | 0.9936 | AGE | 0.00 | not significant |
| GeneE | GeneE | -0.0442 | 0.9936 | AGE | 0.05 | review in larger cohort |

Notes:
- `cohort_frequency` is a simple example metric; replace with study-specific prevalence or sample counts.
- `note` contains an actionable recommendation for experimental validation or follow-up.
