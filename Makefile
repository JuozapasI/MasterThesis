umi = 12
barcode = 16
clusterThreshold = 100
reference = data/gencode.chr12.gtf
fasta = data/12.fa

all: thesis.pdf

%.pdf: %.tex
	latexmk -pdf -silent -deps-out=.depend $*
	
%.bam.bai: %.bam
	samtools index $<
	
data/gene_ranges_sorted.bed: $(reference)
	awk 'BEGIN {FS = "\t"; OFS = FS} { \
		if ($$3 == "gene") { \
		if ($$0 ~ /transcript_id/) {print $$0} \
		else {print $$0" transcript_id \"\""} } }' $(reference) | \
	gtf2bed | \
	sort -k1,1 -k2,2n > $@
	
	
data/%/unassigned_reads/unassigned.bam: data/%/solo_output/Aligned.sortedByCoord.out.bam
	( samtools view -H $^  && samtools view -F 1024 -F 256 $^ | \
	grep "GN:Z:-" | \
	grep -P "\bNH:i:1\b" | \
	grep -E "CB:Z:[A-Z]{$(barcode)}" | \
	grep -E "UB:Z:[A-Z]{$(umi)}" | \
#	awk 'BEGIN {FS = "\t"; OFS=FS} { \
#		split($$0, a, "CB:Z:"); split(a[2], b, " "); cb = b[1]; \
#		split($$0, a, "UB:Z:"); split(a[2], b, " "); ub = b[1]; \
#		if (!seen[cb":"ub]++) {print $$0} }' ; ) | \
	samtools view -h -b - > $@
	
data/%/unassigned_reads/intergenic.bam: data/%/unassigned_reads/unassigned.bam data/gene_ranges_sorted.bed
	bedtools intersect -s -v -abam $< -b data/gene_ranges_sorted.bed > $@
	
data/%/unassigned_reads/intersecting.bam: data/%/unassigned_reads/unassigned.bam data/gene_ranges_sorted.bed
	bedtools intersect -s -u -wa -abam $< -b data/gene_ranges_sorted.bed > $@
	
data/%/unassigned_reads/intergenic.clusters.bed: data/%/unassigned_reads/intergenic.bam
	bedtools bamtobed -i $< | \
	sort -k1,1 -k2,2n | \
	bedtools merge -s -c 6 -o distinct,count > $@
	
data/%/unassigned_reads/intergenic.clusters.thrash.bed: data/%/unassigned_reads/intergenic.clusters.bed
	awk -F'\t' '{if ($$5 <= $(clusterThreshold)) {print $$0}}' $< > $@

data/%/unassigned_reads/intergenic.clusters.good.bed: data/%/unassigned_reads/intergenic.clusters.bed
	awk -F'\t' '{if ($$5 > $(clusterThreshold)) {print $$0}}' $< | sort -k5,5r > $@
	
data/%/unassigned_reads/intergenic.forward.bam: data/%/unassigned_reads/intergenic.bam
	samtools view -h -b -F 16 $< > $@
	
data/%/unassigned_reads/intergenic.backward.bam: data/%/unassigned_reads/intergenic.bam
	samtools view -h -b -f 16 $< > $@

data/%/unassigned_reads/intergenic.coverage.forward.bg: data/%/unassigned_reads/intergenic.forward.bam data/%/unassigned_reads/intergenic.forward.bam.bai
	bamCoverage -b $< -o $@ -of bedgraph -p max
	
data/%/unassigned_reads/intergenic.coverage.backward.bg: data/%/unassigned_reads/intergenic.backward.bam data/%/unassigned_reads/intergenic.backward.bam.bai
	bamCoverage -b $< -o $@ -of bedgraph -p max
	
data/%/unassigned_reads/intergenic.clusters.good.fine.bed: data/%/unassigned_reads/intergenic.clusters.good.bed data/%/unassigned_reads/intergenic.coverage.forward.bg data/%/unassigned_reads/intergenic.coverage.backward.bg
	# get bam coverage of those regions
	# take the largest coverage region of each
	# expand it to take adjecent values that are at least max/2
	# save the list of those 'fine' regions
	python scripts/python/intergenic_regions_tuning.py $^ $@
	
data/%/unassigned_reads/igv.snapshot.batch.txt: data/%/unassigned_reads/intergenic.clusters.good.fine.bed data/%/unassigned_reads/intergenic.bam
	mkdir data/$*/unassigned_reads/intergenic_snapshots/ 
	python scripts/python/igv_batch_script_generator.py $^ \
		  data/$*/unassigned_reads/intergenic_snapshots/ GRCh38.dna.primary_assembly.fa GRCh38.dna.primary_assembly.fa \
		  data/gencode.v47.sorted.gtf data/$*/unassigned_reads/igv.snapshot.batch.txt
	

data/%/unassigned_reads/intergenic.clusters.sequences.tsv: data/%/unassigned_reads/intergenic.clusters.good.fine.bed
	bedtools getfasta -fi $(fasta) -nameOnly -tab -bed $< | cut -f2 > $@
	
data/%/unassigned_reads/intergenic.clusters.GCcontent.tsv: data/%/unassigned_reads/intergenic.clusters.sequences.tsv 
	awk '{l=length($0); gc=gsub(/[GCgc]/,""); print (gc/l)*100 "%"}' $< > $@
		
clean:
	latexmk -c
	rm *.bbl
	rm *.bib.bak
