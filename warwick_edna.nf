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
    file R1 from params.R1
    file R2 from params.R2

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
Remove adaptors and trim low quality reads
========================================================================================
*/

process trim_galore {
    tag "$prefix"

    input:
    tuple val(prefix), file(fastqs)

    output:
    tuple val(prefix), file("${prefix}_R1_trimmed.fq.gz"), file("${prefix}_R2_trimmed.fq.gz")

    script:
    """
    trim_galore \
    --paired \
    --output_dir trimmed \
    --gzip \
    ${fastqs[0]} ${fastqs[1]}
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
    fastq_pairs
        // Run FastQC on the raw FASTQ files
        | fastqc_raw
        // Run Trim Galore to trim and clean FASTQ files
        //| trim_galore
}
