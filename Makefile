umi = 12
barcode = 16
clusterThreshold = 100
reference = data/gencode.chr12.gtf
fasta = data/12.fa
ATrichThreshold = 70

all: thesis.pdf data/chr12/unassigned_reads/intergenic.clusters.GCcontent.tsv

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
	
	
data/%/unassigned_reads/unassigned.bam: data/%/solo_output/Aligned.sortedByCoord.out.bam | data/%/unassigned_reads
	( samtools view -H $^  && samtools view -F 1024 -F 256 $^ | \
	grep "GN:Z:-" | \
	grep -P "\bNH:i:1\b" | \
	grep -E "CB:Z:[A-Z]{$(barcode)}" | \
	grep -E "UB:Z:[A-Z]{$(umi)}" ) | \
#	awk 'BEGIN {FS = "\t"; OFS=FS} { \
#		split($$0, a, "CB:Z:"); split(a[2], b, " "); cb = b[1]; \
#		split($$0, a, "UB:Z:"); split(a[2], b, " "); ub = b[1]; \
#		if (!seen[cb":"ub]++) {print $$0} }' ; ) | \
	samtools view -h -b - > $@
	
data/%/unassigned_reads:
	mkdir -p $@
	
data/%/unassigned_reads/intergenic.bam: data/%/unassigned_reads/unassigned.bam data/gene_ranges_sorted.bed
	bedtools intersect -s -v -abam $< -b data/gene_ranges_sorted.bed > $@
	
data/%/unassigned_reads/intersecting.bam: data/%/unassigned_reads/unassigned.bam data/gene_ranges_sorted.bed
	bedtools intersect -s -u -wa -abam $< -b data/gene_ranges_sorted.bed > $@
	
%.clusters.bed: %.bam
	bedtools bamtobed -i $< | \
	sort -k1,1 -k2,2n | \
	bedtools merge -s -c 6 -o distinct,count > $@
	
%.clusters.thrash.bed: %.clusters.bed
	awk -F'\t' '{if ($$5 <= $(clusterThreshold)) {print $$0}}' $< > $@
	
%.clusters.good.bed: %.clusters.bed
	awk -F'\t' '{if ($$5 > $(clusterThreshold)) {print $$0}}' $< | sort -rnk5,5r > $@
	
%.forward.bam: %.bam
	samtools view -h -b -F 16 $< > $@
	
%.backward.bam: %.bam
	samtools view -h -b -f 16 $< > $@

%.coverage.forward.bg: %.forward.bam %.forward.bam.bai
	bamCoverage -b $< -o $@ -of bedgraph -p max
	
%.coverage.backward.bg: %.backward.bam %.backward.bam.bai
	bamCoverage -b $< -o $@ -of bedgraph -p max
	
%.clusters.good.fine.bed: %.clusters.good.bed %.coverage.forward.bg %.coverage.backward.bg
	python scripts/python/intergenic_regions_tuning.py $^ $@
	
data/%/unassigned_reads/igv.intergenic.snapshot.batch.txt: data/%/unassigned_reads/intergenic.clusters.good.fine.bed data/%/unassigned_reads/intergenic.bam
	mkdir -p data/$*/unassigned_reads/intergenic_snapshots/ 
	python scripts/python/igv_batch_script_generator.py $^ \
		  data/$*/unassigned_reads/intergenic_snapshots/ GRCh38.dna.primary_assembly.fa GRCh38.dna.primary_assembly.fa \
		  data/gencode.v47.sorted.gtf data/$*/unassigned_reads/igv.snapshot.batch.txt
		  
data/%/unassigned_reads/igv.intersecting.snapshot.batch.txt: data/%/unassigned_reads/intersecting.clusters.good.fine.bed data/%/unassigned_reads/intersecting.bam
	mkdir -p data/$*/unassigned_reads/intersecting_snapshots/ 
	python scripts/python/igv_batch_script_generator.py $^ \
		  data/$*/unassigned_reads/intersecting_snapshots/ GRCh38.dna.primary_assembly.fa GRCh38.dna.primary_assembly.fa \
		  data/gencode.v47.sorted.gtf data/$*/unassigned_reads/intersecting.igv.snapshot.batch.txt
	

data/%/unassigned_reads/intergenic.clusters.sequences.tsv: data/%/unassigned_reads/intergenic.clusters.good.fine.bed
	bedtools getfasta -fi $(fasta) -nameOnly -tab -bed $< | cut -f2 > $@
	
data/%/unassigned_reads/intergenic.clusters.GCcontent.tsv: data/%/unassigned_reads/intergenic.clusters.sequences.tsv 
	awk '{l=length($0); gc=gsub(/[GCgc]/,""); print (gc/l)*100}' $< > $@
	
data/%/solo_intergenic_gencode/Aligned.sortedByCoord.out.bam: data/%/unassigned_reads/intergenic.bam data/index_gencode
	bash scripts/solo/starsolo_from_bam.sh data/$*/unassigned_reads/intergenic.bam data/index_gencode/ $@

data/%/solo_intergenic_gencode_ncbi/Aligned.sortedByCoord.out.bam: data/%/unassigned_reads/intergenic.bam data/index_gencode
	bash scripts/solo/starsolo_from_bam.sh data/$*/unassigned_reads/intergenic.bam data/index_gencode/ $@
	


data/%/unassigned_reads/stats.txt: data/%/unassigned_reads/unassigned.bam \
					data/%/unassigned_reads/intersecting.bam \
					data/%/unassigned_reads/intergenic.bam \
					data/%/unassigned_reads/intergenic.clusters.thrash.bed \
					data/%/unassigned_reads/intergenic.clusters.good.bed \
					data/%/unassigned_reads/intergenic.clusters.GCcontent.tsv
	echo -n "Total unassigned (and unique) reads: " > $@
	samtools view data/$*/unassigned_reads/unassigned.bam | wc -l >> $@
	echo "Of them:" >> $@
	echo -n "\tintersecting with genes: " >> $@
	samtools view data/$*/unassigned_reads/intersecting.bam | wc -l >> $@
	echo -n "\tintergenic: " >> $@
	samtools view data/$*/unassigned_reads/intergenic.bam | wc -l >> $@
	echo "Out of the intergenic: " >> $@
	echo -n "\tAre in clusters of size less than or equal to ${clusterThreshold}: " >> $@
	awk -F'\t' '{sum += $$5} END {print sum}' data/$*/unassigned_reads/intergenic.clusters.thrash.bed >> $@
	echo -n "\tAre in clusters of size greater than ${clusterThreshold}: " >> $@
	awk -F'\t' '{sum += $$5} END {print sum}' data/$*/unassigned_reads/intergenic.clusters.good.bed >> $@
	echo "From the cluster larger than ${clusterThreshold}:" >> $@
	echo -n "\tComes from the AT-rich regions (> ${ATrichThreshold}%): " >> $@
	awk -F'\t' 'BEGIN {sum = 0} FNR==NR {percentage[NR]=$$1; next} {if (percentage[FNR] < 100 - ${ATrichThreshold}) sum += $$5} END {print sum}' \
		data/$*/unassigned_reads/intergenic.clusters.GCcontent.tsv \
		data/$*/unassigned_reads/intergenic.clusters.good.bed >> $@
	echo -n "\tComes from non-AT-rich regions: " >> $@
	awk -F'\t' 'BEGIN {sum = 0} FNR==NR {percentage[NR]=$$1; next} {if (percentage[FNR] >= 100 - ${ATrichThreshold}) sum += $$5} END {print sum}' \
		data/$*/unassigned_reads/intergenic.clusters.GCcontent.tsv \
		data/$*/unassigned_reads/intergenic.clusters.good.bed >> $@
	echo "\tMapped to ncbi reference" >> $@

		
clean:
	latexmk -c
	rm *.bbl
	rm *.bib.bak
