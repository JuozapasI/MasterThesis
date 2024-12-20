#!/usr/bin/env nextflow

params.genome_fasta = "${projectDir}/data/genome/fasta/GRCh38.dna.primary_assembly.fa"
params.reference = "${projectDir}/data/genome/references/10x.gtf"

log.info """\
Running intergenic reads analysis.

Parameters:
Genome:
$params.genome_fasta

References:
$params.reference
"""

process STAR_INDEX {
	cpus = 14

	publishDir "${projectDir}", mode: 'copy'

	input:
	path transcriptomic_reference
	

	output:
	path "data/genome/indices/index_10x"

	script:
	"""
	mkdir -p data/genome/indices/index_10x
	echo STAR \
	--runMode genomeGenerate \
	--runThreadN $task.cpus \
	--genomeDir data/genome/indices/index_10x \
	--genomeFastaFiles ${params.genome_fasta} \
	--sjdbGTFfile $transcriptomic_reference >> data/genome/indices/index_10x/script.txt

	"""
}

workflow {

index_ch = STAR_INDEX(params.reference)
index_ch.view()
}

workflow.onComplete {
	log.info (workflow.success ? "\nDone" : "Something went wrong!")
}
