#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
Nextflow pipeline for analysing environmental DNA samples
Author: Ash Sendell-Price
Usage:  nextflow run warwick_edna.nf --input_csv test_samples.csv --blast_db /home/u2271009/eDNA/bin/blast_db/nt --outdir results --taxa 7742 --bin /home/u2271009/warwick_edna_nf/bin -process.cpus 16

========================================================================================

========================================================================================
Define input parameters
========================================================================================
*/

params.input_csv = null
params.blast_db = null
params.taxa = null
params.outdir = "results"

// Error checking for input file
if (params.input_csv == null) {
	error "Please provide the input CSV file with --input_csv parameter."
}
if (params.blast_db == null) {
	error "Please provide the path to blast database with --blast_db parameter"
}

// Print debug information
println "Input CSV: ${params.input_csv}"

/*
========================================================================================
Parse CSV file to extract sample information
========================================================================================
*/

ch_samples = Channel
	.fromPath(params.input_csv)
	.splitCsv(header: true)
	.map { row -> 
		tuple(row.sampleID, file(row.forward_reads), file(row.reverse_reads))
	}

/*
========================================================================================
Run FastQC on the raw FASTQ files
========================================================================================
*/

process fastqc_raw {
	publishDir "${params.outdir}/fastqc/raw", mode: 'copy'
	
	input:
	tuple val(sampleID), path(forward_reads), path(reverse_reads)
	
	output:
	file("${sampleID}*.{html,zip}")

	script:
	"""
	fastqc ${forward_reads} ${reverse_reads}
	"""
}

/*
========================================================================================
Summarise FastQC reports (raw reads) with MultiQC
========================================================================================
*/

process multiqc_raw {
	publishDir "${params.outdir}/multiqc/raw", mode: 'copy'
	
	input:
	path fastqc_reports

	output:
	file("multiqc_report.html")
	file("multiqc_data")

	script:
	"""
	multiqc ${fastqc_reports} -o .
	"""
}

/*
========================================================================================
Remove adaptors and quality trim 
========================================================================================
*/

process trim_adapters {
	publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

	input:
	tuple val(sampleID), path(forward_reads), path(reverse_reads)
	
	output:
	tuple val(sampleID), file("${sampleID}_val_1.fq.gz"), file("${sampleID}_val_2.fq.gz")

	script:
	"""
	trim_galore --paired --output_dir . --basename ${sampleID} ${forward_reads} ${reverse_reads}
	"""
}

/*
========================================================================================
Run FastQC on the trimmed FASTQ files
========================================================================================
*/

process fastqc_trimmed {
	publishDir "${params.outdir}/fastqc/trimmed", mode: 'copy'
	
	input:
	tuple val(sampleID), path(trimmed_forward), path(trimmed_reverse)
	
	output:
	file("${sampleID}_*.{html,zip}")

	script:
	"""
	fastqc -t 10 ${trimmed_forward} ${trimmed_reverse}
	"""
}

/*
========================================================================================
Summarise FastQC reports (trimmed reads) with MultiQC
========================================================================================
*/

process multiqc_trimmed {
	publishDir "${params.outdir}/multiqc/trimmed", mode: 'copy'
	
	input:
	path fastqc_reports

	output:
	file("multiqc_report.html")
	file("multiqc_data")

	script:
	"""
	multiqc ${fastqc_reports} -o .
	"""
}

/*
========================================================================================
Convert fastq files to fasta files and deduplicate fasta sequences
========================================================================================
*/

process seqtk_dedupe {
	publishDir "${params.outdir}/fastas", mode: 'copy'
	
	input:
	tuple val(sampleID), path(trimmed_forward), path(trimmed_reverse)
	
	output:
	tuple val(sampleID), file("${sampleID}.fasta")

	script:
	"""
	seqtk seq -a ${trimmed_forward} ${trimmed_reverse} | dedupe.sh in=stdin out=${sampleID}.fasta ac=f
	"""
}

/*
========================================================================================
Cluster fasta files with MMseqs2
========================================================================================
*/

process mmseqs_cluster {
	publishDir "${params.outdir}/clustering", mode: 'copy'
	
	input:
	tuple val(sampleID), path(fasta_file)
	
	output:
	tuple val(sampleID), file("${sampleID}_rep_seq.fasta")

	script:
	"""
	mmseqs easy-linclust ${fasta_file} ${sampleID} tmp --cov-mode 1 -c 0.9 --min-seq-id 0.97
	"""
}

/*
========================================================================================
Blast sequences against full database
========================================================================================
*/

process blastn {
	publishDir "${params.outdir}/blast", mode: 'copy'
	
	input:
	tuple val(sampleID), path(rep_seq_fasta)
	
	output:
	tuple val(sampleID), file("${sampleID}.blast.txt")

	script:
	"""
	blastn -db ${params.blast_db} \
	-query ${rep_seq_fasta} -out ${sampleID}.blast.txt \
	-max_target_seqs 500 -outfmt "6 std staxids" -num_threads 8
	"""
}

/*
========================================================================================
Blast sequences against selected database
========================================================================================
*/

process blastn_restricted_taxa {
	publishDir "${params.outdir}/blast", mode: 'copy'
	
	input:
	tuple val(sampleID), path(rep_seq_fasta)
	
	output:
	tuple val(sampleID), file("${sampleID}.${params.taxa}.blast.txt")

	script:
	"""
	blastn -db ${params.blast_db} \
	-query ${rep_seq_fasta} -out ${sampleID}.${params.taxa}.blast.txt \
	-max_target_seqs 500 -outfmt "6 std staxids" -num_threads 8 \
	-taxids ${params.taxa}
	"""
}

/*
========================================================================================
PIA analysis
========================================================================================
*/

process pia {
	publishDir "${params.outdir}/pia", mode: 'copy'
	
	input:
	tuple val(sampleID), path(rep_seq_fasta), path(blast_out)
	
	output:
	file("*.Full.txt")
	file("*.PIA_inner_logs.txt")
	file("*.Post-PIA.fasta")
	file("*.Summary_Basic.txt")
	file("*.Summary_Reads_MEGAN.txt")
	file("*.Summary_Reads.txt")

	script:
	"""
	cp -r /home/u2271009/warwick_edna_nf/bin/PIA/* .
	perl PIA.pl -f ${rep_seq_fasta} -b ${blast_out}
	mv *.PIA_output/* .
	"""
}

/*
========================================================================================
Merge pia "_basic.txt" files
========================================================================================
*/

process merge_pia {
	publishDir "${params.outdir}/pia", mode: 'copy'
	
	input:
	path pia_out_files from pia_ch
	
	output:
	file("taxaIDs.txt")
	file("taxa_info.txt")
	file("pia_merged_incl_taxa.txt")

	script:
	"""
	#!/usr/bin/env python

	#Import required modules
	import pandas as pd
	import os

	#Create a list of PIA "Summary_Basic.txt" files  
	keyword = '_Basic.txt'
	current_dir = os.getcwd()  # Get the current working directory
	files = [f for f in os.listdir(current_dir) if os.path.isfile(f) and keyword in f]

	#Load first PIA file in list as pandas df 
	sampleID = files[0].split(".")[0]
	df = pd.read_csv(files[0], sep='\t', skiprows=11)
	df = df.rename(columns={"Reads" : sampleID})

	#Load subsequent PIA files and merge
	for i in range(1,len(files)-1):
		sampleID = files[i].split(".")[0]
		df2 = pd.read_csv(files[i], sep='\t', skiprows=11)
		df2 = df2.rename(columns={"Reads" : sampleID})
		df = pd.merge(df, df2, how = 'outer')

	#Replae NAs with zero and rename column 1 to "taxa_ID"
	df = df.fillna(0)
	df = df.rename(columns={"# ID": "taxa_ID"})

	# Add taxonomic information
	df["taxa_ID"].to_csv('taxaIDs.txt', sep ='\t', index=False, header=False)
	os.system('taxaranks -i taxaIDs.txt -o taxa_info.txt')
	taxa = pd.read_csv('taxa_info.txt', sep='\t').rename(columns={"user_taxa": "taxa_ID"})
	df = pd.merge(df, taxa, how = 'outer')

	# Write dataframe to file
	df.to_csv("pia_merged_incl_taxa.txt", sep ='\t', index=False)

	"""
}

/*
========================================================================================
Workflow
========================================================================================
*/

workflow {
	// Define the input channel from the CSV
	ch_samples = Channel
		.fromPath(params.input_csv)
		.splitCsv(header: true)
		.map { row -> 
			tuple(row.sampleID, file(row.forward_reads), file(row.reverse_reads))
		}

	// Execute FASTQC process on raw reads and summarise with multiqc
	raw_fastqc_out = fastqc_raw(ch_samples).collect()
	multiqc_raw(raw_fastqc_out)

	// Execute TRIMGALORE process to trim adapters
	trimmed_reads_ch = trim_adapters(ch_samples)

	// Execute FASTQC process on trimmed reads and summarise with multiqc
	trimmed_fastqc_out = fastqc_trimmed(trimmed_reads_ch).collect()
	multiqc_trimmed(trimmed_fastqc_out)

	// Convert trimmed FASTQ files to FASTA and deduplicate sequences
	deduped_fasta_ch = seqtk_dedupe(trimmed_reads_ch)

	// Cluster deduplicated FASTA sequences with MMseqs2
	clustered_sequences_ch = mmseqs_cluster(deduped_fasta_ch)

	// Blast sequences against ncbi database
	if (params.taxa == null) {
		blast_output_ch = blastn(clustered_sequences_ch)
	} else {
		blast_output_ch = blastn_restricted_taxa(clustered_sequences_ch)
	}

	// Combine the outputs of mmseqs_cluster and blastn processes
	combined_ch = clustered_sequences_ch.join(blast_output_ch)

	// Run PIA (and collect output files)
	pia_ch = pia(combined_ch).concat()

	// Merge PIA files
	//merge_pia(pia_ch)
}
