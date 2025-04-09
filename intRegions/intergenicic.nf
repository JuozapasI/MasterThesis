#!/usr/bin/env nextflow

params.bam = ""                          // Input BAM file
params.gtf = ""                          // Transcriptomic reference GTF file
params.cpm_threshold = 5                 // Custom cluster size threshold (in cpm)
params.threshold = null                  // Custom cluster size threshold
params.read_length = 1000                // Max allowed read length
params.unassigned_pattern = "GN:Z:-"     // Pattern to filter out 

process getReadCount {
    input:
    path bam_file

    output:
    stdout

    script:
    """
    samtools view -F 0x900 '$bam_file' | wc -l 
    """
}

process extractUnassignedReads {
    input:
    path bam_file

    output:
    path "unassigned.bam"

    script:
    """
    ( samtools view -H '$bam_file' ; \
    samtools view -F 0x900 '$bam_file' | grep -F "${params.unassigned_pattern}" ) |  samtools view -bS - > unassigned.bam
    """
}

process filterAgainstGTF {
    input:
    path unassigned_bam
    path gtf_file

    output:
    path "intergenic.bam"

    script:
    """
    bedtools intersect -s -v -abam ${unassigned_bam} -b <(sort -k1,1 -k4,4n ${gtf_file}) > intergenic.bam
    """
}

process clusterReads {
    input:
    path intergenic_bam
    val threshold

    output:
    path "clusters.bed"

    script:
    """
    bedtools bamtobed -i ${intergenic_bam} | \
    	awk -F '\t' -v len=${params.read_length} '\$3 - \$2 < len' | \
	sort -k1,1 -k2,2n | \
	bedtools merge -s -c 6 -o distinct,count | \
	sort -k1,1 -k2,2n | \
	awk -F'\t' 'BEGIN {OFS = FS} \$5 > ${threshold} {print \$1, \$2, \$3, ".", \$5, \$4}' > clusters.bed
    """
}

process annotateClusters {
    input:
    path clusters
    path gtf_file

    output:
    path "annotated_clusters.tsv"
    
    publishDir "${baseDir}", mode: 'copy'

    script:
    """
    
    bedtools closest -wa -t first -D a -id -s \
    -a ${clusters} -b <(sort -k1,1 -k4,4n ${gtf_file} | awk -F '\t' '\$3 == "gene"') | \
    awk -F '\t' 'BEGIN {OFS = FS} {if (\$0 ~ "gene_name") {split(\$0, a, "gene_name"); split(a[2], b, "\\""); gene = b[2]} \
    else if (\$0 ~ "gene_id") {split(\$0, a, "gene_id"); split(a[2], b, "\\""); gene = b[2]} \
    else {gene = "."} print \$1, \$2, \$3, \$4, \$5, \$6, gene, \$NF}' | sort -k1,1 -k2,2n > closest_upstream.bed
    
    bedtools closest -wa -t first -D a -iu -s \
    -a ${clusters} -b <(sort -k1,1 -k4,4n ${gtf_file} | awk -F '\t' '\$3 == "gene"') | \
    awk -F '\t' 'BEGIN {OFS = FS} {if (\$0 ~ "gene_name") {split(\$0, a, "gene_name"); split(a[2], b, "\\""); gene = b[2]} \
    else if (\$0 ~ "gene_id") {split(\$0, a, "gene_id"); split(a[2], b, "\\""); gene = b[2]} \
    else {gene = "."} print \$1, \$2, \$3, \$4, \$5, \$6, gene, \$NF}' | sort -k1,1 -k2,2n > closest_downstream.bed
    
    bedtools closest -wa -t first -D a -id -S \
    -a ${clusters} -b <(sort -k1,1 -k4,4n ${gtf_file} | awk -F '\t' '\$3 == "gene"') | \
    awk -F '\t' 'BEGIN {OFS = FS} {if (\$0 ~ "gene_name") {split(\$0, a, "gene_name"); split(a[2], b, "\\""); gene = b[2]} \
    else if (\$0 ~ "gene_id") {split(\$0, a, "gene_id"); split(a[2], b, "\\""); gene = b[2]} \
    else {gene = "."} print \$1, \$2, \$3, \$4, \$5, \$6, gene, \$NF}' | sort -k1,1 -k2,2n > closest_upstream_antisense.bed
    
    bedtools closest -wa -t first -D a -iu -S \
    -a ${clusters} -b <(sort -k1,1 -k4,4n ${gtf_file} | awk -F '\t' '\$3 == "gene"') | \
    awk -F '\t' 'BEGIN {OFS = FS} {if (\$0 ~ "gene_name") {split(\$0, a, "gene_name"); split(a[2], b, "\\""); gene = b[2]} \
    else if (\$0 ~ "gene_id") {split(\$0, a, "gene_id"); split(a[2], b, "\\""); gene = b[2]} \
    else {gene = "."} print \$1, \$2, \$3, \$4, \$5, \$6, gene, \$NF}' | sort -k1,1 -k2,2n > closest_downstream_antisense.bed
    
    paste closest_upstream.bed <(cut -f 7,8 closest_downstream.bed) <(cut -f 7,8 closest_upstream_antisense.bed)\
    <(cut -f 7,8 closest_downstream_antisense.bed) > annotated_clusters.tsv
    """
}

// Check if required parameters are provided
if (!params.bam) {
    error "The BAM file must be specified using --bam."
}

if (!params.gtf) {
    error "The GTF file must be specified using --gtf."
}

workflow {
    // Input BAM and GTF files
    Channel
        .fromPath(params.bam)
        .set { bam_file }

    Channel
        .fromPath(params.gtf)
        .set { gtf_file }


    // Defining threshold
    if (!params.threshold) {
    // Only get the read count if threshold is not manually set
    getReadCount(bam_file)
        .map { count -> count.toInteger() }
        .set { read_count }

    // Calculate threshold based on read count
    read_count
        .map { count ->
            (count / 1_000_000 * params.cpm_threshold) as int
        }
        .set { threshold }
        
    threshold
        .view { threshold_value -> println "Threshold value: $threshold_value" }

    } else {
        // Use the manually set threshold
        threshold = params.threshold as int
        println "Threshold value: $threshold"
    }


    // Continue with the rest of the processes
    extractUnassignedReads(bam_file)
        .set { unassigned_bam }

    filterAgainstGTF(unassigned_bam, gtf_file)
        .set { non_overlapping_bed }

    clusterReads(non_overlapping_bed, threshold)
        .set { clusters }

    annotateClusters(clusters, gtf_file)

}

