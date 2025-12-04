#!/usr/bin/env nextflow

// Nextflow template: Somatic variant calling (skeleton)

params.samplesheet = params.samplesheet ?: 'samplesheet.csv'
params.outdir = params.outdir ?: 'results'

workflow {
    Channel.fromPath(params.samplesheet)
    // ... Define channels and processes here
    println "This is a skeleton Nextflow pipeline. Fill with processes: QC, align, postprocessing, variant_calling, annotation"
}
