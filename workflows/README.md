# Workflow templates and usage

This folder contains example workflow definitions and helper scripts. The primary `Snakefile` demonstrates a typical bulk RNA-Seq workflow (preprocess -> DESeq2 -> annotate -> plots) and includes two additional rules for clinical processing and integration:

- `clinical_process`: cleans raw clinical CSV and produces a simple analysis summary (uses `workflows/clinical_pipeline.py`).
- `integrate_genomics_clinical`: merges annotated gene results with clinical summaries (`workflows/integrate_genomics_clinical.py`).

Run the Snakemake workflow from the repository root:

```bash
snakemake --snakefile workflows/Snakefile --cores 1
```

Customize source paths and tools as needed for your environment and data standards (CDISC) before running at scale.

