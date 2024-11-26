umi = 12
barcode = 16
clusterThreshold = 100
fasta = data/12.fa
ATrichThreshold = 70
AT = 60

references = 10x gencode ncbi lnc
first =  $(firstword $(references))
last = $(lastword $(references))
without_first = $(wordlist 2, $(words $(references)), $(references))

order = $(first)
order += $(foreach n, $(without_first), $(lastword $(order)).$(n))

order_without_first = $(wordlist 2, $(words $(order)), $(order))

.SECONDARY:

all: thesis.pdf data/PBMC_10x/unassigned_reads/stats.txt

%.pdf: %.tex
	latexmk -pdf -silent -deps-out=.depend $*

%.bam.bai: %.bam
	samtools index $<

# Extracting unassigned reads
$(foreach f, $(order), $(eval \
data/%/unassigned_reads/$(f).unassigned.bam: data/%/solo_output.$(f)/Aligned.sortedByCoord.out.bam | data/%/unassigned_reads; \
bash scripts/bash/take_unassigned.sh $$< $$@ $(barcode) $(umi)))

# Generating gene location bed file from gtf
data/%.gene_ranges_sorted.bed: data/%.gtf
	awk 'BEGIN {FS = "\t"; OFS = FS} { \
		if ($$0 ~ /transcript_id/) {print $$0} \
		else {print $$0" transcript_id \"\""} }' $< | \
	gtf2bed | \
	sort -k1,1 -k2,2n > $@
	
data/%/unassigned_reads:
	mkdir -p $@

# Extracting intergenic reads
$(foreach f, $(order_without_first), $(eval \
data/%/unassigned_reads/$(f).intergenic.bam: data/%/unassigned_reads/$(f).unassigned.bam \
data/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gene_ranges_sorted.bed; \
bedtools intersect -s -v -abam $$< -b data/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gene_ranges_sorted.bed > $$@))
	
# Extracting intersecting reads
$(foreach f, $(order), $(eval \
data/%/unassigned_reads/$(f).intersecting.bam: data/%/unassigned_reads/$(f).unassigned.bam \
data/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gene_ranges_sorted.bed; \
bedtools intersect -s -u -wa -abam $$< -b data/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gene_ranges_sorted.bed > $$@))

# Separate cases for the first reference, as we additionally filter out AT rich sequences
data/%/unassigned_reads/$(first).intergenic.bam: data/%/unassigned_reads/$(first).unassigned.bam data/%/$(first).gene_ranges_sorted.bed
	bedtools intersect -s -v -abam $< -b data/$*/$(first).gene_ranges_sorted.bed | \
	( samtools view -H - ; samtools view - | grep -E -v "A{60}|T{60}" ; ) | samtools view -h -b - > $@
	bedtools intersect -s -v -abam $< -b data/$*/$(first).gene_ranges_sorted.bed | \
	( samtools view -H - ; samtools view - | grep -E "A{60}|T{$(AT)}" ; ) | samtools view -h -b - > AT_seq.bam


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
	
$(foreach f, $(order_without_first), $(eval \
data/%/solo_output.$(f)/Aligned.sortedByCoord.out.bam: data/%/unassigned_reads/$(basename $(f)).intergenic.bam; \
bash scripts/solo/starsolo_from_bam.sh $$< data/index_$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f)))/ data/$$*/solo_output.$(f)/))
	











# Summary
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
