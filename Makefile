umi = 12
barcode = 16
clusterThreshold = 100
fasta = data/GRCh38.dna.primary_assembly.fa
ATrichThreshold = 70
AT = 60
ends_dist = 1000
closeEndThreshold = 100
cpmThreshold = 5

seq_depth = 200000000

extensionCountThreshold = $(shell echo $$(($(cpmThreshold) * $(seq_depth) / 1000000 / 5)))
newGeneCountThreshold = $(shell echo $$(($(cpmThreshold) * $(seq_depth) / 1000000)))


references = 10x gencode ncbi lnc
first =  $(firstword $(references))
last = $(lastword $(references))
without_first = $(wordlist 2, $(words $(references)), $(references))

order = 10x 10x.gencode 10x.gencode.ncbi 10x.gencode.ncbi.lnc

order_without_first = $(wordlist 2, $(words $(order)), $(order))
order_last = $(lastword $(order))

debug:
	echo $(order)
	echo $(first)
	echo $(order_last)

.SECONDARY:

all: data/PBMC_10x/unassigned_reads/Summary.txt data/brain/unassigned_reads/Summary.txt \
data/brain/unassigned_reads/final.gtf data/brain/solo_output.final/Aligned.sortedByCoord.out.bam

%.pdf: %.tex
	latexmk -pdf -silent -deps-out=.depend $*

%.bam.bai: %.bam
	samtools index $<

# Extracting unassigned reads
$(foreach f, $(order), $(eval \
data/%/unassigned_reads/$(f).unassigned.bam: data/%/solo_output.$(f)/Aligned.sortedByCoord.out.bam | data/%/unassigned_reads; \
bash scripts/bash/take_unassigned.sh $$< $$@ $(barcode) $(umi)))

# Generating gene location bed file from gtf
%.gene_ranges_sorted.bed: %.gtf
	awk 'BEGIN {FS = "\t"; OFS = FS} { \
		if ($$3 == "gene") { \
		if ($$0 ~ /transcript_id/) {print $$0} \
		else {print $$0" transcript_id \"\""} } }' $< | \
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
data/%/unassigned_reads/$(first).intergenic.bam data/%/unassigned_reads/AT_seq.bam: data/%/unassigned_reads/$(first).unassigned.bam data/$(first).gene_ranges_sorted.bed
	bedtools intersect -s -v -abam $< -b data/$(first).gene_ranges_sorted.bed > tmp_bam
	( samtools view -H tmp_bam ; samtools view tmp_bam | grep -E -v "A{60}|T{$(AT)}" ; ) | \
	samtools view -h -b - > data/$*/unassigned_reads/$(first).intergenic.bam
	bedtools intersect -s -v -abam $< -b data/$(first).gene_ranges_sorted.bed > tmp_bam
	( samtools view -H tmp_bam ; samtools view tmp_bam | grep -E "A{60}|T{$(AT)}" ; ) | \
	samtools view -h -b - > data/$*/unassigned_reads/AT_seq.bam
	rm tmp_bam


%.clusters.bed: %.bam
	bedtools bamtobed -i $< | \
	sort -k1,1 -k2,2n | \
	bedtools merge -s -c 6 -o distinct,count | \
	sort -k1,1 -k2,2n > $@
	
%.clusters.thrash.bed: %.clusters.bed
	awk -F'\t' '{if ($$5 <= $(clusterThreshold)) {print $$0}}' $< > $@
	
%.clusters.good.bed: %.clusters.bed
	awk -F'\t' '{if ($$5 > $(clusterThreshold)) {print $$0}}' $< > $@

$(foreach f, $(order), $(eval \
data/%/unassigned_reads/$(f).intersecting_genes.bed: data/%/unassigned_reads/$(f).intersecting.clusters.good.fine.bed \
data/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gene_ranges_sorted.bed; \
bedtools intersect -s -wa -wb -abam $$< -b data/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gene_ranges_sorted.bed | \
sort -k1,1 -k2,2n > $$@))

%.intersecting_gene_list.tsv: %.intersecting_genes.bed
	cut -f 11 $< | sort -u > $@
	
$(foreach f, $(order), $(eval \
data/%/unassigned_reads/$(f).overlappers.csv: data/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gtf \
data/%/unassigned_reads/$(f).intersecting_gene_list.tsv; \
Rscript scripts/R/Overlappers.R $$^ $$@))

$(foreach f, $(order), $(eval \
data/%/unassigned_reads/$(f).overlaps_modified.gtf: data/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gtf \
data/%/unassigned_reads/$(f).overlappers.csv; \
Rscript scripts/R/ResolveOverlappers.R $$^ $(ends_dist) $$@))

# Assemble modified gtf for the first reference (i.e. with resolved overlaps)
data/%/unassigned_reads/$(first).modified.gtf: data/%/unassigned_reads/$(first).intersecting_gene_list.tsv \
data/$(first).gtf data/%/unassigned_reads/$(first).overlaps_modified.gtf
	( grep -v -f $< data/$(first).gtf; cat data/$*/unassigned_reads/$(first).overlaps_modified.gtf ; ) | \
	sort -k1,1 -k4,4n > $@

# For the second and on, we need to make sure that appended entries doesn't overlap with previous gtf
$(foreach f, $(order_without_first), $(eval \
data/%/unassigned_reads/$(f).new_entries.gtf: data/%/unassigned_reads/$(basename $(f)).modified.gtf \
data/%/unassigned_reads/$(f).overlaps_modified.gtf; \
Rscript scripts/R/ResolveGtfsOverlaps.R $$^ $$@))

# And then we can merge gtfs
$(foreach f, $(order_without_first), $(eval \
data/%/unassigned_reads/$(f).modified.gtf: data/%/unassigned_reads/$(basename $(f)).modified.gtf \
data/%/unassigned_reads/$(f).new_entries.gtf; \
( cat $$< ; cat data/$$*/unassigned_reads/$(f).new_entries.gtf ; ) | \
	sort -k1,1 -k4,4n > $$@))
	
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
	

data/%/unassigned_reads/$(order_last).intergenic.clusters.sequences.tsv: data/%/unassigned_reads/$(order_last).intergenic.clusters.good.fine.bed
	bedtools getfasta -fi $(fasta) -nameOnly -tab -bed $< | cut -f2 > $@
	
data/%/unassigned_reads/$(order_last).intergenic.clusters.GCcontent.tsv: data/%/unassigned_reads/$(order_last).intergenic.clusters.sequences.tsv 
	awk '{l=length($0); gc=gsub(/[GCgc]/,""); print (gc/l)*100}' $< > $@
	
$(foreach f, $(order_without_first), $(eval \
data/%/solo_output.$(f)/Aligned.sortedByCoord.out.bam: data/%/unassigned_reads/$(basename $(f)).intergenic.bam; \
bash scripts/solo/starsolo_from_bam.sh $$< data/index_$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f)))/ data/$$*/solo_output.$(f)/))
	
# Filter those clusters that came from AT rich regions:
%.filteredAT.intergenic.clusters.good.fine.bed: %.intergenic.clusters.good.fine.bed %.intergenic.clusters.GCcontent.tsv
	paste -d ',' $^ | awk -F',' '$$2 > 30 {print $$1}' > $@

#In final part, check if there are any clusters very close to genes 3' ends and extend those genes if there are
%.distances.clusters.good.fine.tsv: %.filteredAT.intergenic.clusters.good.fine.bed %.modified.gene_ranges_sorted.bed
	bedtools closest -a $< -b $*.modified.gene_ranges_sorted.bed -s -D a -id -fu > $@

# Split clusters into close to 3' ends and not
data/%/unassigned_reads/final.intergenic.clusters_near_genes.txt: data/%/unassigned_reads/$(order_last).distances.clusters.good.fine.tsv
	awk -F'\t' '{if((-1 * $$18) < $(closeEndThreshold) && $$11 != ".") {print $$0}}' $< > $@

data/%/unassigned_reads/final.intergenic.clusters_far_from_genes.txt: data/%/unassigned_reads/$(order_last).distances.clusters.good.fine.tsv
	awk -F'\t' '{if((-1 * $$18) >= $(closeEndThreshold)) {print $$0}}' $< > $@

# Make gene extension list
data/%/unassigned_reads/final.extension_candidates.csv: data/%/unassigned_reads/final.intergenic.clusters_near_genes.txt
	awk 'BEGIN {FS = "\t"; OFS = FS;} {if ($$7 >= $(extensionCountThreshold)) {print $$1, $$2, $$3, $$6, $$11}}' $< > $@

# Make new intergenic gene list:
data/%/unassigned_reads/final.new_gene_list.csv: data/%/unassigned_reads/final.intergenic.clusters_far_from_genes.txt
	awk 'BEGIN {FS = "\t"; OFS = FS;} {if ($$7 >= $(newGeneCountThreshold)) {print $$1, $$2, $$3, $$6}}' $< > $@

# Assemble final gtf
%/final.gtf: %/$(order_last).modified.gtf %/final.extension_candidates.csv %/final.new_gene_list.csv
	Rscript scripts/R/final_gtf.R $^ $@

# Make final genome index
data/%/index_final/: data/%/unassigned_reads/final.gtf 
	bash scripts/solo/genome_index.sh $< $(fasta) $@


# Run starsolo on final gtf
data/%/solo_output.final/Aligned.sortedByCoord.out.bam: data/%/solo_output.$(first)/Aligned.sortedByCoord.out.bam \
data/%/index_final/
	bash scripts/solo/starsolo_from_bam.sh $^ data/$*/solo_output.final/


# Make igv snapshots script for extension candidates
data/%/unassigned_reads/igv_snapshots/extension_candidates.batch.txt: data/%/unassigned_reads/final.extension_candidates.csv \
data/%/solo_output.10x/Aligned.sortedByCoord.out.bam
	mkdir -p data/$*/unassigned_reads/igv_snapshots/extension_candidates/ 
	python scripts/python/igv_batch_script_generator.py $^ \
		  data/$*/unassigned_reads/igv_snapshots/extension_candidates/ GRCh38.dna.primary_assembly.fa GRCh38.dna.primary_assembly.fa \
		  data/gencode.v47.sorted.gtf $@

# Make igv snaphots script for new gene candidates
data/%/unassigned_reads/igv_snapshots/new_gene_candidates.batch.txt: data/%/unassigned_reads/final.new_gene_list.csv \
data/%/solo_output.10x/Aligned.sortedByCoord.out.bam
	mkdir -p data/$*/unassigned_reads/igv_snapshots/new_gene_candidates/ 
	python scripts/python/igv_batch_script_generator.py $^ \
		  data/$*/unassigned_reads/igv_snapshots/new_gene_candidates/ GRCh38.dna.primary_assembly.fa GRCh38.dna.primary_assembly.fa \
		  data/gencode.v47.sorted.gtf $@
		  
# Make igv snaphots script for overlaps (first reference)
data/%/unassigned_reads/igv_snapshots/overlaps.batch.txt: data/%/unassigned_reads/$(first).intersecting.clusters.good.fine.bed \
data/%/solo_output.10x/Aligned.sortedByCoord.out.bam
	mkdir -p data/$*/unassigned_reads/igv_snapshots/overlaps/ 
	python scripts/python/igv_batch_script_generator.py $^ \
		  data/$*/unassigned_reads/igv_snapshots/overlaps/ GRCh38.dna.primary_assembly.fa GRCh38.dna.primary_assembly.fa \
		  data/gencode.v47.sorted.gtf $@
		  
		  
		  

# Summary
data/%/unassigned_reads/Summary.txt: data/%/solo_output.10x/Aligned.sortedByCoord.out.bam \
data/%/unassigned_reads/10x.unassigned.bam \
data/%/unassigned_reads/10x.intersecting.bam \
data/%/unassigned_reads/10x.intergenic.bam \
data/%/unassigned_reads/AT_seq.bam \
data/%/unassigned_reads/10x.gencode.ncbi.lnc.intergenic.clusters.good.bed \
data/%/unassigned_reads/10x.gencode.ncbi.lnc.intergenic.clusters.thrash.bed \
data/%/unassigned_reads/10x.gencode.unassigned.bam \
data/%/unassigned_reads/10x.gencode.intersecting.bam \
data/%/unassigned_reads/10x.gencode.intergenic.bam \
data/%/unassigned_reads/10x.gencode.ncbi.unassigned.bam \
data/%/unassigned_reads/10x.gencode.ncbi.intersecting.bam \
data/%/unassigned_reads/10x.gencode.ncbi.intergenic.bam \
data/%/unassigned_reads/10x.gencode.ncbi.lnc.unassigned.bam \
data/%/unassigned_reads/10x.gencode.ncbi.lnc.intersecting.bam \
data/%/unassigned_reads/10x.gencode.ncbi.lnc.intergenic.bam 
	bash scripts/bash/statistics.sh $^ > $@
	

# Captured gene list:
data/%/unassigned_reads/captured_genes_$(first).csv: data/%/solo_output.$(first)/Aligned.sortedByCoord.out.bam
	samtools view $< | grep -o -P "GN:Z:[^ \t]*" | cut -d':' -f3 | sort | uniq > $@
	
data/%/unassigned_reads/captured_genes_final.csv: data/%/solo_output.final/Aligned.sortedByCoord.out.bam
	samtools view $< | grep -o -P "GN:Z:[^ \t]*" | cut -d':' -f3 | sort | uniq > $@




		
clean:
	latexmk -c
	rm *.bbl
	rm *.bib.bak
