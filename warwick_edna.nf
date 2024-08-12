#!/usr/bin/env nextflow

/*
========================================================================================
Nextflow pipeline for analysing environmental DNA samples
========================================================================================

Author: Ash Sendell-Price

========================================================================================
Define initial files
========================================================================================
*/

params.version = "0.1.0"

// Initialize parameters from the command line
params.R1 = file(params.R1)
params.R2 = file(params.R2)
params.prefix = params.prefix

/*
========================================================================================
Create FASTQ Pairs Channel
========================================================================================
*/

Channel
    .of(tuple(params.prefix, [params.R1, params.R2]))
    .set { fastq_pairs }


/*
========================================================================================
Run FastQC on the raw FASTQ files
========================================================================================
*/

process fastqc_raw {
    tag "$prefix"

    input:
    tuple val(prefix), file(R1), file(R2)

    output:
    file("${prefix}_R1_fastqc.zip")
    file("${prefix}_R1_fastqc.html")
    file("${prefix}_R2_fastqc.zip")
    file("${prefix}_R2_fastqc.html")

    script:
    """
    fastqc \
    ${R1} ${R2}
    """
}


/*
========================================================================================
Workflow
========================================================================================
*/

workflow {
    // Check if required parameters are defined
    if (!params.R1 || !params.R2 || !params.prefix) {
        error "Missing required parameters. Ensure R1, R2, and prefix are provided."
    }

    // Create the channel for FASTQ pairs
    def fastq_pairs = Channel
        .of(tuple(params.prefix, [params.R1, params.R2]))

    // Run FastQC on the raw FASTQ files
    fastq_pairs
        | fastqc_raw

}