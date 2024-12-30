#!/usr/bin/env nextflow

params.genome_fasta = "${projectDir}/data/genome/fasta/GRCh38.dna.primary_assembly.fa"
params.references = "${projectDir}/data/genome/references/list.tsv"
params.barcodes_10x = "${projectDir}/data/barcode_whitelists/10x_v3.1_barcodes"

def references_list =  file(params.references)
def ref_names = []
references_list.eachLine { line ->
    def name = line.split('\t')[0]
    ref_names << name
}

log.info """\
Running intergenic reads analysis.

Parameters:
Genome:
$params.genome_fasta

References:
$ref_names

"""

process STAR_INDEX {
	cpus 14
	time '2d'
	memory '120 GB'

	publishDir "${projectDir}", mode: 'copy'

	input:
	tuple val(name), path(transcriptomic_reference)
	
	output:
	tuple val(name), path("data/genome/indices/index_${name}")

	script:
	"""
	mkdir -p data/genome/indices/${name}
	STAR \
	--runMode genomeGenerate \
	--runThreadN $task.cpus \
	--genomeDir data/genome/indices/index_${name} \
	--genomeFastaFiles ${params.genome_fasta} \
	--sjdbGTFfile $transcriptomic_reference

	"""
}

process STARSOLO {
	cpus 14
	memory '120 GB'
	time '2d'

	publishDir "${projectDir}", mode: 'copy'	

	input:
	tuple val(name), val(type), path(read1), path(read2), val(index_name), path(index)
	
	output:
	tuple val(name), path("data/datasets/${name}/solo_output_${index_name}")

	script:
	if (type == '10x_v3.1') {
	"""
	echo STAR 10x
	"""
	} else if (type == 'indrops2') {
	"""
	echo STAR indrops
	"""
	} else {
	throw new IllegalArgumentException("Sample type ${type} is not supported yet.")
	}
}

process SAMTOOLS_INDEX {
	input:
	path bam
	
	output:
	path "*.bai"
	
	script:
	"""
	samtools index ${bam}
	"""
}

workflow {

Channel
	.fromPath(params.references)
	.splitCsv(sep: '\t')
	.map { tuple(it[0], it[1]) }
	.set { references_ch }

references_ch.view()

indices_ch = STAR_INDEX(references_ch)

indices_ch.toSortedList { tuple1, tuple2 ->
        int idx1 = ref_names.indexOf(tuple1[0])
        int idx2 = ref_names.indexOf(tuple2[0])
        return idx1 <=> idx2 }

def ch = [Channel.of(1)]

// Loop to create new channels
for (int i = 1; i < 10; i++) {
    ch << ch[i - 1].flatMap { it + 1}
}
ch.view()
}

workflow.onComplete {
	log.info (workflow.success ? "\nDone" : "Something went wrong!")
}
