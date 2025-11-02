#!/usr/bin/env nextflow

// Example Nextflow workflow for RNA-Seq preprocessing

process PREPROCESS_QC {
    input:
    path counts_file
    output:
    path 'counts_for_deseq2.csv'
    path 'library_sizes.png'
    script:
    """
    python ../rna_seq/preprocess_qc.py
    """
}

workflow {
    PREPROCESS_QC(counts_file: file('../rna_seq/sample_counts.csv'))
}
