#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
Define input parameters
========================================================================================
*/

params.input_csv = null
params.taxa = null
params.outdir = "results"

// Error checking for input file
if (params.input_csv == null) {
	error "Please provide the input CSV file with --input_csv parameter."
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
	tuple val(sampleID), 
		path("${sampleID}_rep_seq.Full.txt"),
		path("${sampleID}_rep_seq.PIA_inner_logs.txt"),
		path("${sampleID}_rep_seq.Summary_Basic.txt"),
		path("${sampleID}_rep_seq.Summary_Reads_MEGAN.txt"),
		path("${sampleID}_rep_seq.Summary_Reads.txt")

	script:
	"""
	cp -r $params.pia_dir/* .
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
	path summary_files

	output:
	file("taxaIDs.txt")
	file("taxa_info.txt")
	file("pia_merged_incl_taxa.txt")

	script:
	"""
	#!/usr/bin/env python3

	import pandas as pd
	import os

	# Create a list of PIA "Summary_Basic.txt" files from the input provided by Nextflow
	files = ${summary_files.collect { "'${it}'" }.join(',')}

	# Load the first PIA file as a pandas DataFrame
	first_file = files[0].strip("'")
	sampleID = os.path.basename(first_file).split(".")[0]
	df = pd.read_csv(first_file, sep='\\t', skiprows=11)
	df = df.rename(columns={"Reads": sampleID})

	# Iterate through remaining files and merge them
	for file in files[1:]:
		file = file.strip("'")
		sampleID = os.path.basename(file).split(".")[0]
		df2 = pd.read_csv(file, sep='\\t', skiprows=11)
		df2 = df2.rename(columns={"Reads": sampleID})
		df = pd.merge(df, df2, how="outer", on="# ID")

	# Replace NAs with zero and rename column 1 to "taxa_ID"
	df = df.fillna(0)
	df = df.rename(columns={"# ID": "taxa_ID"})

	# Save taxa_IDs for taxonomy enrichment
	df["taxa_ID"].to_csv('taxaIDs.txt', sep='\\t', index=False, header=False)

	# Run the taxaranks command to add taxonomic information
	os.system('taxaranks -i taxaIDs.txt -o taxa_info.txt')

	# Load taxonomy info and merge with the main dataframe
	taxa = pd.read_csv('taxa_info.txt', sep='\\t').rename(columns={"user_taxa": "taxa_ID"})
	df = pd.merge(df, taxa, how="outer", on="taxa_ID")

	# Write the final merged output to a file
	df.to_csv("pia_merged_incl_taxa.txt", sep='\\t', index=False)
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

	// Run PIA (plus collect output files and flatten into a single channel)
	pia_ch = pia(combined_ch)

	// Split PIA outputs to create `pia_merge_in_ch` for `_rep_seq.Summary_Basic.txt` files
	pia_merge_in_ch = pia_ch.map { sampleID, full_txt, inner_logs, summary_basic, summary_megan, summary_reads -> summary_basic}.collect()

	// Merge PIA results
	merge_pia(pia_merge_in_ch)
}
