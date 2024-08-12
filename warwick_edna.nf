#!/usr/bin/env nextflow

/*
========================================================================================
Nextflow pipeline for analysing environmental DNA samples
========================================================================================

Author: Ash Sendell-Price

========================================================================================
Define initial parameters
========================================================================================
*/

params.version = "0.1.0"

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
    trim_galore --paired --output_dir trimmed --basename ${prefix} ${fastqs[0]} ${fastqs[1]}
    """
}


/*
========================================================================================
Workflow
========================================================================================
*/

workflow {
    fastq_pairs | trim_galore
}
