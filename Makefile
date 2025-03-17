export LC_ALL=C
SHELL := /bin/bash

umi = 12
barcode = 16
# for PBMC_indrops umi = 10 barcode = 28
# for PBMC_indrops_2 umi=8 barcode=17
# for eye umi=12 barcode=19
clusterThreshold = 100
fasta = data/genome/fasta/GRCh38.dna.primary_assembly.fa
ATrichThreshold = 70
AT = 60
ends_dist = 1000
closeEndThreshold = 100
cpmThreshold = 5

ref_dir = data/genome/references/
references_paths = $(ref_dir)10x.gene_ranges_sorted.bed $(ref_dir)gencode.gene_ranges_sorted.bed $(ref_dir)ncbi.gene_ranges_sorted.bed $(ref_dir)lnc.gene_ranges_sorted.bed

references = 10x gencode ncbi lnc
first =  $(firstword $(references))
last = $(lastword $(references))
without_first = $(wordlist 2, $(words $(references)), $(references))

order = 10x 10x.gencode 10x.gencode.ncbi 10x.gencode.ncbi.lnc

order_without_first = $(wordlist 2, $(words $(order)), $(order))
order_last = $(lastword $(order))

datasets = PBMC_10x.dataset brain.dataset PBMC_10x_2.dataset PBMC_indrops.dataset

#$(foreach dataset, $(datasets), $(eval $(dataset)_seq_depth = \
#$(shell samtools flagstat --threads 8 data/$(dataset)/solo_output.$(first)/Aligned.sortedByCoord.out.bam | head -2 | tail -1 | cut -d ' ' -f1)))
#
#$(foreach dataset, $(datasets), $(eval $(dataset)_extensionCountThreshold = \
#$(shell echo $$(($(cpmThreshold) * $($(dataset)_seq_depth) / 1000000 / 5)))))
#
#$(foreach dataset, $(datasets), $(eval $(dataset)_newGeneCountThreshold = \
#$(shell echo $$(($(cpmThreshold) * $($(dataset)_seq_depth) / 1000000)))))

# Defining tresholds (instead can be used above script to compute them automatically,
# however, in such case, these thresholds are recalculated every time you run make, which takes a bit of time)

PBMC_10x_seq_depth = 182330834
PBMC_10x_2_seq_depth = 496387931
PBMC_10x_3_seq_depth = 368640939
brain_seq_depth = 206360627
brain_2_seq_depth = 122556503
PBMC_indrops_seq_depth = 112932507
PBMC_indrops_2_seq_depth = 471705924
eye_seq_depth = 375397270
eye_2_seq_depth = 140981808
eye_3_seq_depth = 161261977
lung_2_seq_depth = 511080104
lung_5_seq_depth = 452105505
lung_7_seq_depth = 524095146
lung_8_seq_depth = 342092138
PBMC_diabetes_1_seq_depth = 107430657
PBMC_diabetes_2_seq_depth = 75410796

PBMC_10x_extensionCountThreshold = 182
PBMC_10x_2_extensionCountThreshold = 496
PBMC_10x_3_extensionCountThreshold = 368
brain_extensionCountThreshold = 206
brain_2_extensionCountThreshold = 122
PBMC_indrops_extensionCountThreshold = 112
PBMC_indrops_2_extensionCountThreshold = 472
eye_extensionCountThreshold = 375
eye_2_extensionCountThreshold = 141
eye_3_extensionCountThreshold = 161
lung_2_extensionCountThreshold = 511
lung_5_extensionCountThreshold = 452
lung_7_extensionCountThreshold = 524
lung_8_extensionCountThreshold = 342
PBMC_diabetes_1_extensionCountThreshold = 107
PBMC_diabetes_2_extensionCountThreshold = 75

PBMC_10x_newGeneCountThreshold = 911
PBMC_10x_2_newGeneCountThreshold = 2481
PBMC_10x_3_newGeneCountThreshold = 1843
brain_newGeneCountThreshold = 1031
brain_2_newGeneCountThreshold = 613
PBMC_indrops_newGeneCountThreshold = 564
PBMC_indrops_2_newGeneCountThreshold = 2358
eye_newGeneCountThreshold = 1876
eye_2_newGeneCountThreshold = 705
eye_3_newGeneCountThreshold = 806
lung_2_newGeneCountThreshold = 2555
lung_5_newGeneCountThreshold = 2261
lung_7_newGeneCountThreshold = 2620
lung_8_newGeneCountThreshold = 1710
PBMC_diabetes_1_newGeneCountThreshold = 535
PBMC_diabetes_2_newGeneCountThreshold = 375

debug:
	echo $(references_paths)

# Comment the line below out, if you want intermediate files to be cleaned afterwards automatically
.SECONDARY:

#.PHONY: $(datasets)

all: data/downstream/summaries/gene_summaries/intersecting_gene_summary.tex \
data/downstream/summaries/count_summaries/count_summary.tex \
data/downstream/summaries/gene_summaries/intersecting_gene_summary.tex


# rule to generate all outputs for one dataset
%.dataset: data/downstream/summaries/count_summaries/%.csv \
     data/downstream/summaries/gene_summaries/%.csv \
     data/downstream/summaries/captured_gene_summaries/%.final.csv \
     data/downstream/summaries/captured_gene_summaries/%.$(first).csv \
     data/downstream/matrices/%/10x/matrix.mtx.gz \
     data/downstream/matrices/%/10x/features.tsv.gz \
     data/downstream/matrices/%/10x/barcodes.tsv.gz \
     data/downstream/matrices/%/final/matrix.mtx.gz \
     data/downstream/matrices/%/final/features.tsv.gz \
     data/downstream/matrices/%/final/barcodes.tsv.gz \
     data/downstream/intergenic/%.csv
	@echo "Done with $* dataset"

# rule to combine outputs from different datasets
# important: should be runned only after all 'datasets', i.e. it doesn't check if outputs from all samples are present,
# combines only those that are present
downstream: data/downstream/intergenic/closest.bed \
     data/downstream/summaries/captured_gene_summaries/captured_gene_types_summary.tex \
     data/downstream/summaries/count_summaries/count_summary.tex \
     data/downstream/summaries/gene_summaries/intersecting_gene_summary.tex \
     data/downstream/igv/intergenic.batch.txt \
     data/downstream/intergenic/intergenic.gtf \
     data/downstream/intergenic/index/ \
     data/downstream/intergenic/augustus/predictions/ \
     data/downstream/intergenic/augustus/predictions.gtf

# rule to map unassigned reads using combined intergenic annotation
%.dataset2: data/downstream/matrices/intergenic/%/matrix.mtx.gz \
     data/downstream/matrices/intergenic/%/features.tsv.gz \
     data/downstream/matrices/intergenic/%/barcodes.tsv.gz
	@echo "$* done."

# rule for generating thesis pdf
%.pdf: %.tex
	latexmk -pdf -silent -deps-out=.depend $*

# indexing bam files
%.bam.bai: %.bam
	samtools index $<

# Extracting unassigned reads
$(foreach f, $(order), $(eval \
data/datasets/%/unassigned_reads/$(f).unassigned.bam: data/datasets/%/solo_output.$(f)/Aligned.sortedByCoord.out.bam | data/datasets/%/unassigned_reads; \
bash scripts/bash/take_unassigned.sh $$< $$@ $(barcode) $(umi)))

# Generating gene location bed file from gtf
%.gene_ranges_sorted.bed: %.gtf
	awk 'BEGIN {FS = "\t"; OFS = FS} { \
		if ($$3 == "gene") { \
		if ($$0 ~ /transcript_id/) {print $$0} \
		else {print $$0"; transcript_id \"\""} } }' $< | \
	gtf2bed --do-not-sort | \
	sort -k1,1 -k2,2n > $@

data/datasets/%/unassigned_reads:
	mkdir -p $@

# Extracting intergenic reads
$(foreach f, $(order_without_first), $(eval \
data/datasets/%/unassigned_reads/$(f).intergenic.bam: data/datasets/%/unassigned_reads/$(f).unassigned.bam \
data/genome/references/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gene_ranges_sorted.bed; \
bedtools intersect -s -v -abam $$< -b data/genome/references/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gene_ranges_sorted.bed > $$@))
	
# Extracting intersecting reads
$(foreach f, $(order), $(eval \
data/datasets/%/unassigned_reads/$(f).intersecting.bam: data/datasets/%/unassigned_reads/$(f).unassigned.bam \
data/genome/references/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gene_ranges_sorted.bed; \
bedtools intersect -s -u -wa -abam $$< -b data/genome/references/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gene_ranges_sorted.bed > $$@))

# Separate cases for the first reference, as we additionally filter out AT rich sequences
data/datasets/%/unassigned_reads/$(first).intergenic.bam data/datasets/%/unassigned_reads/AT_seq.bam: data/datasets/%/unassigned_reads/$(first).unassigned.bam data/genome/references/$(first).gene_ranges_sorted.bed
	bedtools intersect -s -v -abam $< -b data/genome/references/$(first).gene_ranges_sorted.bed > tmp_bam_$*
	( samtools view -H tmp_bam_$* ; samtools view tmp_bam_$* | grep -E -v "A{60}|T{$(AT)}" ; ) | \
	samtools view -h -b - > data/datasets/$*/unassigned_reads/$(first).intergenic.bam
	bedtools intersect -s -v -abam $< -b data/genome/references/$(first).gene_ranges_sorted.bed > tmp_bam_$*
	( samtools view -H tmp_bam_$* ; samtools view tmp_bam_$* | grep -E "A{60}|T{$(AT)}" ; ) | \
	samtools view -h -b - > data/datasets/$*/unassigned_reads/AT_seq.bam
	rm tmp_bam_$*

# rule for clustering reads
%.clusters.bed: %.bam
	bedtools bamtobed -i $< | \
	sort -k1,1 -k2,2n | \
	bedtools merge -s -c 6 -o distinct,count | \
	sort -k1,1 -k2,2n | \
	awk -F'\t' 'BEGIN {OFS = FS} {print $$1, $$2, $$3, ".", $$5, $$4}' > $@

# rule for taking big enough clusters of reads
# initially we take clusters of at least 100 reads, later, when constructing intergenic region list,
# more strict filtering (5 CPM) is applied.
%.clusters.good.bed: %.clusters.bed
	awk -F'\t' '{if ($$5 > $(clusterThreshold)) {print $$0}}' $< > $@

# getting list of genes intersecting with unassigned reads
$(foreach f, $(order), $(eval \
data/datasets/%/unassigned_reads/$(f).intersecting_genes.bed: data/datasets/%/unassigned_reads/$(f).intersecting.clusters.good.fine.bed \
data/genome/references/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gene_ranges_sorted.bed; \
bedtools intersect -s -wa -wb -abam $$< -b data/genome/references/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gene_ranges_sorted.bed | \
sort -k1,1 -k2,2n > $$@))

# Extracting intersecting genes
%.intersecting_gene_list.tsv: %.intersecting_genes.bed
	cut -f 11 $< | sort -u > $@

# Getting list of overlapping genes that contain unassigned reads
$(foreach f, $(order), $(eval \
data/datasets/%/unassigned_reads/$(f).overlappers.csv: data/genome/references/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gtf \
data/datasets/%/unassigned_reads/$(f).intersecting_gene_list.tsv; \
Rscript scripts/R/Overlappers.R $$^ $$@))

# Constructing GTF of the adjusted overlapping genes
$(foreach f, $(order), $(eval \
data/datasets/%/unassigned_reads/$(f).overlaps_modified.gtf: data/genome/references/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gtf \
data/datasets/%/unassigned_reads/$(f).overlappers.csv; \
Rscript scripts/R/ResolveOverlappers.R $$^ $(ends_dist) $$@))

# Assemble modified gtf for the first reference (i.e. with resolved overlaps)
data/datasets/%/unassigned_reads/$(first).modified.gtf: data/datasets/%/unassigned_reads/$(first).intersecting_gene_list.tsv \
data/genome/references/$(first).gtf data/datasets/%/unassigned_reads/$(first).overlaps_modified.gtf
	( grep -v -f $< data/genome/references/$(first).gtf; cat data/datasets/$*/unassigned_reads/$(first).overlaps_modified.gtf ; ) | \
	sort -k1,1 -k4,4n > $@

# Assembling adjusted gtf for second and on
# In this case, we need to make sure that appended entries doesn't overlap with previous gtf
$(foreach f, $(order_without_first), $(eval \
data/datasets/%/unassigned_reads/$(f).new_entries.gtf: data/datasets/%/unassigned_reads/$(basename $(f)).modified.gtf \
data/datasets/%/unassigned_reads/$(f).overlaps_modified.gtf; \
Rscript scripts/R/ResolveGtfsOverlaps.R $$^ $$@))

# Merging GTFs (previously generated + new entries)
$(foreach f, $(order_without_first), $(eval \
data/datasets/%/unassigned_reads/$(f).modified.gtf: data/datasets/%/unassigned_reads/$(basename $(f)).modified.gtf \
data/datasets/%/unassigned_reads/$(f).new_entries.gtf; \
( cat $$< ; cat data/datasets/$$*/unassigned_reads/$(f).new_entries.gtf ; ) | \
	sort -k1,1 -k4,4n > $$@))

# rule for taking reads only on the positive strand
%.forward.bam: %.bam
	samtools view -h -b -F 16 $< > $@

# rule for taking reads only on the negative strand
%.backward.bam: %.bam
	samtools view -h -b -f 16 $< > $@

# rule for compute coverage of reads on the positive strand
%.coverage.forward.bg: %.forward.bam %.forward.bam.bai
	bamCoverage -b $< -o $@ -of bedgraph -p max
	
# rule for compute coverage of reads on the negative strand
%.coverage.backward.bg: %.backward.bam %.backward.bam.bai
	bamCoverage -b $< -o $@ -of bedgraph -p max

# rule for adjusting intergenic regions boundaries
%.clusters.good.fine.bed: %.clusters.good.bed %.coverage.forward.bg %.coverage.backward.bg
	python scripts/python/intergenic_regions_tuning.py $^ $@
	
# rule to extract sequences of intergenic regions
data/datasets/%/unassigned_reads/$(order_last).intergenic.clusters.sequences.tsv: data/datasets/%/unassigned_reads/$(order_last).intergenic.clusters.good.fine.bed
	bedtools getfasta -fi $(fasta) -nameOnly -tab -bed $< | cut -f2 > $@

# rule to compute GC content of intergenic regions
data/datasets/%/unassigned_reads/$(order_last).intergenic.clusters.GCcontent.tsv: data/datasets/%/unassigned_reads/$(order_last).intergenic.clusters.sequences.tsv 
	awk '{l=length($0); gc=gsub(/[GCgc]/,""); print (gc/l)*100}' $< > $@

# rule to map unassigned reads not intersecting with genes from current reference using the next reference
$(foreach f, $(order_without_first), $(eval \
data/datasets/%/solo_output.$(f)/Aligned.sortedByCoord.out.bam: data/datasets/%/unassigned_reads/$(basename $(f)).intergenic.bam \
data/genome/indices/index_$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f)))/; \
bash scripts/solo/starsolo_from_bam.sh $$< data/genome/indices/index_$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f)))/ data/datasets/$$*/solo_output.$(f)/ $(umi)))
	
# rule to filter those clusters that came from AT rich regions:
%.filteredAT.intergenic.clusters.good.fine.bed: %.intergenic.clusters.good.fine.bed %.intergenic.clusters.GCcontent.tsv
	paste -d ',' $^ | awk -F',' '$$2 > 30 {print $$1}' > $@

# rule to find distances from clusters to closest 3' ends of genes
%.distances.clusters.good.fine.tsv: %.filteredAT.intergenic.clusters.good.fine.bed %.modified.gene_ranges_sorted.bed
	bedtools closest -t first -a $< -b $*.modified.gene_ranges_sorted.bed -s -D a -id -fu > $@

# rule to split clusters into close to 3' ends of genes and not
data/datasets/%/unassigned_reads/final.intergenic.clusters_near_genes.txt: data/datasets/%/unassigned_reads/$(order_last).distances.clusters.good.fine.tsv
	awk -F'\t' '{if((-1 * $$18) < $(closeEndThreshold) && $$11 != ".") {print $$0}}' $< > $@

# rule to take those clusters that are not very close to 3' ends
data/datasets/%/unassigned_reads/final.intergenic.clusters_far_from_genes.txt: data/datasets/%/unassigned_reads/$(order_last).distances.clusters.good.fine.tsv
	awk -F'\t' '{if((-1 * $$18) >= $(closeEndThreshold)) {print $$0}}' $< > $@

# rule to make gene extension list
data/datasets/%/unassigned_reads/final.extension_candidates.csv: data/datasets/%/unassigned_reads/final.intergenic.clusters_near_genes.txt
	awk 'BEGIN {FS = "\t"; OFS = FS;} {if ($$7 >= $($*_extensionCountThreshold)) {print $$1, $$2, $$3, $$6, $$11}}' $< > $@

# rule to make new intergenic gene list:
data/datasets/%/unassigned_reads/final.new_gene_list.bed: data/datasets/%/unassigned_reads/final.intergenic.clusters_far_from_genes.txt
	awk -v depth="$($*_seq_depth)" -v threshold="$($*_newGeneCountThreshold)" 'BEGIN {FS = "\t"; OFS = FS;} \
	{if ($$7 >= threshold) {i += 1; print $$1, $$2, $$3, "INTERGENIC" i, $$7 * 1000000 / depth, $$6}}' $< > $@

# rule to assemble final gtf (modified entries + intergenic genes)
# %/final.gtf: %/$(order_last).modified.gtf %/final.extension_candidates.csv %/final.new_gene_list.bed
%/final.gtf: %/$(order_last).modified.gtf %/final.extension_candidates.csv
	Rscript scripts/R/final_gtf.R $^ $@

# rule to make final genome index
data/datasets/%/index_final/: data/datasets/%/unassigned_reads/final.gtf 
	bash scripts/solo/genome_index.sh $< $(fasta) $@

# rule to make STAR indices from gtf files
$(foreach reference, $(references), $(eval data/genome/indices/index_$(reference)/: data/genome/references/$(reference).gtf ; bash scripts/solo/genome_index.sh $$< $(fasta) $$@ ))

# rule to run starsolo on final gtf
data/datasets/%/solo_output.final/Aligned.sortedByCoord.out.bam: data/datasets/%/solo_output.$(first)/Aligned.sortedByCoord.out.bam \
data/datasets/%/index_final/
	bash scripts/solo/starsolo_from_bam.sh $^ data/datasets/$*/solo_output.final/ $(umi)


# Make igv snapshots script for extension candidates
#data/datasets/%/unassigned_reads/igv_snapshots/extension_candidates.batch.txt: data/datasets/%/unassigned_reads/final.extension_candidates.csv \
#data/datasets/%/solo_output.10x/Aligned.sortedByCoord.out.bam
#	mkdir -p data/datasets/$*/unassigned_reads/igv_snapshots/extension_candidates/ 
#	python scripts/python/igv_batch_script_generator.py $^ \
#		  data/datasets/$*/unassigned_reads/igv_snapshots/extension_candidates/ \
#		  GRCh38.dna.primary_assembly.fa GRCh38.dna.primary_assembly.fa \
#		  data/gencode.v47.sorted.gtf $@

# Make igv snaphots script for new gene candidates
data/downstream/igv/intergenic.batch.txt: data/downstream/intergenic/filtered.bed 
	mkdir -p data/downstream/igv/
	python scripts/python/igv_batch_script_generator.py $^ \
		  bam_file \
		  images/igv/intergenic/ \
		  GRCh38.dna.primary_assembly.fa GRCh38.dna.primary_assembly.fa \
		  data/gencode.v47.sorted.gtf $@
		  
# Make igv snaphots script for overlaps (first reference)
#data/datasets/%/unassigned_reads/igv_snapshots/overlaps.batch.txt: data/datasets/%/unassigned_reads/$(first).intersecting.clusters.good.fine.bed \
#data/datasets/%/solo_output.10x/Aligned.sortedByCoord.out.bam
#	mkdir -p data/datasets/$*/unassigned_reads/igv_snapshots/overlaps/ 
#	python scripts/python/igv_batch_script_generator.py $^ \
#		  data/datasets/$*/unassigned_reads/igv_snapshots/overlaps/ GRCh38.dna.primary_assembly.fa GRCh38.dna.primary_assembly.fa \
#		  data/gencode.v47.sorted.gtf $@
		  
		  
		  

# Gene count summary
data/datasets/%/unassigned_reads/Summary.txt: data/datasets/%/solo_output.10x/Aligned.sortedByCoord.out.bam \
data/datasets/%/unassigned_reads/10x.unassigned.bam \
data/datasets/%/unassigned_reads/10x.intersecting.bam \
data/datasets/%/unassigned_reads/10x.intergenic.bam \
data/datasets/%/unassigned_reads/AT_seq.bam \
data/datasets/%/unassigned_reads/10x.gencode.ncbi.lnc.intergenic.clusters.good.bed \
data/datasets/%/unassigned_reads/10x.gencode.unassigned.bam \
data/datasets/%/unassigned_reads/10x.gencode.intersecting.bam \
data/datasets/%/unassigned_reads/10x.gencode.intergenic.bam \
data/datasets/%/unassigned_reads/10x.gencode.ncbi.unassigned.bam \
data/datasets/%/unassigned_reads/10x.gencode.ncbi.intersecting.bam \
data/datasets/%/unassigned_reads/10x.gencode.ncbi.intergenic.bam \
data/datasets/%/unassigned_reads/10x.gencode.ncbi.lnc.unassigned.bam \
data/datasets/%/unassigned_reads/10x.gencode.ncbi.lnc.intersecting.bam \
data/datasets/%/unassigned_reads/10x.gencode.ncbi.lnc.intergenic.bam 
	bash scripts/bash/statistics.sh $^ $($*_newGeneCountThreshold) > $@

# rule to compute conservation scores for intergenic regions
%/conservation_scores.csv: %/final.new_gene_list.bed data/genome/conservation/hg38.phastCons100way.bed
	bedtools map -a $< -b data/genome/conservation/hg38.phastCons100way.bed -c 5 -o mean > $@

# rule to make captured gene list:
data/datasets/%/unassigned_reads/captured_genes_$(first).csv: data/datasets/%/solo_output.$(first)/Aligned.sortedByCoord.out.bam
	samtools view -F 1024 -F 256 $< | grep -o -P "GX:Z:[^ \t]*" | cut -d':' -f3 | sort | uniq > $@
# rule to make captured gene list (with final annotation):	
data/datasets/%/unassigned_reads/captured_genes_final.csv: data/datasets/%/solo_output.final/Aligned.sortedByCoord.out.bam
	samtools view -F 1024 -F 256 $< | grep -o -P "GX:Z:[^ \t]*" | cut -d':' -f3 | sort | uniq > $@

# rule to extract gene types from gtf:
%.gene_types.tsv: %.gtf
	bash scripts/bash/extract_gene_types.sh  $< > $@

# rule to make captured gene list:
data/datasets/%/unassigned_reads/captured_gene_types.$(first).csv: data/datasets/%/unassigned_reads/captured_genes_$(first).csv \
data/genome/references/$(first).gene_types.tsv
	grep -f data/datasets/$*/unassigned_reads/captured_genes_$(first).csv data/genome/references/$(first).gene_types.tsv | \
	cut -f 2 | sort | uniq -c | sort -nrk1,1 > $@
# rule to make captured gene list (final annotation):
data/datasets/%/unassigned_reads/captured_gene_types.final.csv: data/datasets/%/unassigned_reads/captured_genes_final.csv \
data/datasets/%/unassigned_reads/final.gene_types.tsv
	grep -f data/datasets/$*/unassigned_reads/captured_genes_final.csv data/datasets/$*/unassigned_reads/final.gene_types.tsv | \
	cut -f 2 | sort | uniq -c | sort -nrk1,1 > $@

# rule to make list of genes that were not detected using initial references, but are detected using final reference
data/datasets/%/unassigned_reads/additional_gene_summary.txt: data/datasets/%/unassigned_reads/$(first).intergenic.clusters.good.bed $(references_paths)
	bedtools intersect -wb -s -a $< -b $(references_paths) -names $(references) | \
	awk -F '\t' '{split($$17, a, "gene_type \""); split(a[2], b, "\""); print $$7, b[1];}' | \
	sort | uniq -c | sort -k2,2 -k1,1nr > $@



### Data arrangement for downstream analyses ###
# Collect and compress matrices
data/downstream/matrices/%/10x/matrix.mtx.gz: data/datasets/%/solo_output.10x/Aligned.sortedByCoord.out.bam
	mkdir -p data/downstream/matrices/$*/10x/
	cp data/datasets/$*/solo_output.10x/Solo.out/GeneFull/raw/matrix.mtx data/downstream/matrices/$*/10x/matrix.mtx
	gzip -f data/downstream/matrices/$*/10x/matrix.mtx
	
data/downstream/matrices/%/10x/features.tsv.gz: data/datasets/%/solo_output.10x/Aligned.sortedByCoord.out.bam
	mkdir -p data/downstream/matrices/$*/10x/
	cp data/datasets/$*/solo_output.10x/Solo.out/GeneFull/raw/features.tsv data/downstream/matrices/$*/10x/features.tsv
	gzip -f data/downstream/matrices/$*/10x/features.tsv
	
data/downstream/matrices/%/10x/barcodes.tsv.gz: data/datasets/%/solo_output.10x/Aligned.sortedByCoord.out.bam
	mkdir -p data/downstream/matrices/$*/10x/
	cp data/datasets/$*/solo_output.10x/Solo.out/GeneFull/raw/barcodes.tsv data/downstream/matrices/$*/10x/barcodes.tsv
	gzip -f data/downstream/matrices/$*/10x/barcodes.tsv
	
data/downstream/matrices/%/final/matrix.mtx.gz: data/datasets/%/solo_output.final/Aligned.sortedByCoord.out.bam
	mkdir -p data/downstream/matrices/$*/final/
	cp data/datasets/$*/solo_output.final/Solo.out/GeneFull/raw/matrix.mtx data/downstream/matrices/$*/final/matrix.mtx
	gzip -f data/downstream/matrices/$*/final/matrix.mtx
	
data/downstream/matrices/%/final/features.tsv.gz: data/datasets/%/solo_output.final/Aligned.sortedByCoord.out.bam
	mkdir -p data/downstream/matrices/$*/final/
	cp data/datasets/$*/solo_output.final/Solo.out/GeneFull/raw/features.tsv data/downstream/matrices/$*/final/features.tsv
	gzip -f data/downstream/matrices/$*/final/features.tsv
	
data/downstream/matrices/%/final/barcodes.tsv.gz: data/datasets/%/solo_output.final/Aligned.sortedByCoord.out.bam
	mkdir -p data/downstream/matrices/$*/final/
	cp data/datasets/$*/solo_output.final/Solo.out/GeneFull/raw/barcodes.tsv data/downstream/matrices/$*/final/barcodes.tsv
	gzip -f data/downstream/matrices/$*/final/barcodes.tsv

# Collect summaries:
data/downstream/summaries/count_summaries/%.csv: data/datasets/%/unassigned_reads/Summary.txt
	mkdir -p data/downstream/summaries/count_summaries/
	cp $< $@
	
data/downstream/summaries/gene_summaries/%.csv: data/datasets/%/unassigned_reads/additional_gene_summary.txt
	mkdir -p data/downstream/summaries/gene_summaries/
	cp $< $@

data/downstream/summaries/captured_gene_summaries/%.final.csv: data/datasets/%/unassigned_reads/captured_gene_types.final.csv
	mkdir -p data/downstream/summaries/captured_gene_summaries/
	cp $< $@
	
data/downstream/summaries/captured_gene_summaries/%.$(first).csv: data/datasets/%/unassigned_reads/captured_gene_types.$(first).csv
	mkdir -p data/downstream/summaries/captured_gene_summaries/
	cp $< $@
	
# Collect intergenic regions info
data/downstream/intergenic/%.csv: data/datasets/%/unassigned_reads/conservation_scores.csv
	mkdir -p data/downstream/intergenic/
	cp $< $@
	


# Combine count summaries into one latex table
data/downstream/summaries/count_summaries/count_summary.tex: data/downstream/summaries/count_summaries/
	python scripts/python/combine_count_summaries.py $< $@
	
# Combine intersecting genes summaries into one latex table
data/downstream/summaries/gene_summaries/intersecting_gene_summary.tex: data/downstream/summaries/gene_summaries/
	python scripts/python/combine_intersecting_gene_summaries.py $< $@
	
# Combine captured genes summaries into one latex table
data/downstream/summaries/captured_gene_summaries/captured_gene_types_summary.tex: data/downstream/summaries/captured_gene_summaries/
	python scripts/python/combine_captured_gene_summaries.py $< $@

# Combine intergenic regions
data/downstream/intergenic/combined.bed:
	awk -F '\t' 'BEGIN {OFS = FS} !seen[$$1, $$2, $$3, FILENAME]++ {print $$1, $$2, $$3, substr(FILENAME, 28, length(FILENAME) - 31) "." FNR, $$5, $$6, $$7}' data/downstream/intergenic/*.csv | \
	sort -k1,1 -k2,2n | \
	awk -F'\t' 'BEGIN {OFS = FS} {if($$7 == ".") {print $$1, $$2, $$3, $$4, $$5, $$6, 0} else {print $$0}}' | \
	bedtools merge -s -c 4,5,6,6,7 -o collapse,mean,distinct,count,mean -i - | \
	sort -k7,7nr -k5,5nr > $@

# Check if intergenic regions overlap with prediction tools
data/downstream/intergenic/predictions.bed: data/downstream/intergenic/combined.bed
	bedtools intersect -s -loj -wa -wb -a $< -b data/genome/predictions/*.bed -names Augustus Geneid Gescan SIB SPG | \
	awk -F'\t' 'BEGIN {OFS = FS} !seen[$$4 ":" $$9]++ {print $$1, $$2, $$3, $$4, $$5, $$6, $$7, $$8, $$9}' | \
	awk -F'\t' '{if($$4 in data) {cols[$$4] = cols[$$4] "," $$9;} \
	else {data[$$4] = $$1 "\t" $$2 "\t" $$3 "\t" $$4 "\t" $$5 "\t" $$6 "\t" $$7 "\t" $$8; cols[$$4] = $$9;}} \
	END {for (key in data) {print data[key] "\t" cols[key];}}' | sort -k7,7nr -k5,5nr > $@
	
# More comprehensive intersection list, if needed
data/downstream/intergenic/predictions_full_prediction_entries.bed: data/downstream/intergenic/combined.bed
	bedtools intersect -s -wa -wb -a $< -b data/genome/predictions/*.bed -names Augustus Geneid Gescan SIB SPG > $@
	
# Combine all those intergenic regions from various samples into one gtf annotaion
data/downstream/intergenic/intergenic.gtf: data/downstream/intergenic/filtered.bed
	awk -F '\t' 'BEGIN {OFS=FS} {print $$1, "scRNAseqData", "gene", $$2, $$3, ".", $$6, ".", \
	"gene_id \"" $$4 "\"; gene_name \"" $$4 "\"; no_samples " $$7 "; mean_cpm " $$5 "; aliases \"" $$14 "\"; conservation_score " $$8 \
	"; predictions \"" $$9 "\";"; \
	print $$1, "scRNAseqData", "exon", $$2, $$3, ".", $$6, ".", \
	"gene_id \"" $$4 "\"; gene_name \"" $$4 "\"; no_samples " $$7 "; mean_cpm " $$5 "; aliases \"" $$14 "\"; conservation_score " $$8 \
	"; predictions \"" $$9 "\";";}' $< > $@

# Make solo index for this annotation
data/downstream/intergenic/index/: data/downstream/intergenic/intergenic.gtf
	bash scripts/solo/genome_index.sh $< $(fasta) $@
	
# Run starsolo with combined intergenic annotations for all samples (use only unassigned sequences)
data/datasets/%/solo_output.intergenic/Aligned.sortedByCoord.out.bam: data/datasets/%/unassigned_reads/$(order_last).intergenic.bam \
     data/downstream/intergenic/index/
	bash scripts/solo/starsolo_from_bam.sh $^ data/datasets/$*/solo_output.intergenic/ $(umi)
	
# Colect intergenic matrices
data/downstream/matrices/intergenic/%/matrix.mtx.gz: data/datasets/%/solo_output.intergenic/Aligned.sortedByCoord.out.bam
	mkdir -p data/downstream/matrices/intergenic/$*/
	cp data/datasets/$*/solo_output.intergenic/Solo.out/GeneFull/raw/matrix.mtx data/downstream/matrices/intergenic/$*/matrix.mtx
	gzip -f data/downstream/matrices/intergenic/$*/matrix.mtx
	
data/downstream/matrices/intergenic/%/features.tsv.gz: data/datasets/%/solo_output.intergenic/Aligned.sortedByCoord.out.bam
	mkdir -p data/downstream/matrices/intergenic/$*/
	cp data/datasets/$*/solo_output.intergenic/Solo.out/GeneFull/raw/features.tsv data/downstream/matrices/intergenic/$*/features.tsv
	gzip -f data/downstream/matrices/intergenic/$*/features.tsv
	
data/downstream/matrices/intergenic/%/barcodes.tsv.gz: data/datasets/%/solo_output.intergenic/Aligned.sortedByCoord.out.bam
	mkdir -p data/downstream/matrices/intergenic/$*/
	cp data/datasets/$*/solo_output.intergenic/Solo.out/GeneFull/raw/barcodes.tsv data/downstream/matrices/intergenic/$*/barcodes.tsv
	gzip -f data/downstream/matrices/intergenic/$*/barcodes.tsv
	
	
# Update combined intergenic gene list to contain distances to the closest genes
data/downstream/intergenic/closest.bed: data/downstream/intergenic/predictions.bed
	sort -k1,1 -k2,2n $< | \
	bedtools closest -mdb all -t first -d -s -a stdin -b data/genome/references/gencode.gene_ranges_sorted.bed \
	data/genome/references/ncbi.gene_ranges_sorted.bed data/genome/references/lnc.gene_ranges_sorted.bed -names gencode ncbi lnc | \
	awk -F '\t' 'BEGIN {OFS=FS} {split($$20, a, "gene_name"); split (a[2], b, "\""); \
	print $$1, $$2, $$3, $$4, $$5, $$6, $$7, $$8, $$9, $$10 ":" b[2], $$21}' | \
	bedtools closest -mdb all -t first -d -S -a stdin -b data/genome/references/gencode.gene_ranges_sorted.bed \
	data/genome/references/ncbi.gene_ranges_sorted.bed data/genome/references/lnc.gene_ranges_sorted.bed -names gencode ncbi lnc | \
	awk -F '\t' 'BEGIN {OFS=FS} {split($$22, a, "gene_name"); split (a[2], b, "\""); \
	print $$1, $$2, $$3, $$4, $$5, $$6, $$7, $$8, $$9, $$10, $$11, $$12 ":" b[2], $$23}' | \
	sort -k7,7nr -k5,5rn | \
	awk -F'\t' 'BEGIN {OFS=FS} {print $$1, $$2, $$3, "INT" NR, $$5, $$6, $$7, $$8, $$9, $$10, $$11, $$12, $$13, "aliases:" $$4}' > $@
	
# extend intergenic regions for the augustus prediction tool
data/downstream/intergenic/extended_intervals.bed: data/downstream/intergenic/predictions.bed
	awk -F '\t' 'BEGIN{OFS=FS} {if ($$6 == "+") {start = $$2 - 10000; end = $$3 + 1000} else {start = $$2 - 1000; end = $$3 + 10000}; \
	print $$1, start, end, "INT" NR, ".", $$6}' $< > $@

# get sequences of intergenic regions
data/downstream/intergenic/augustus/sequences/: data/downstream/intergenic/extended_intervals.bed
	mkdir -p $@
	bedtools getfasta -name -fi $(fasta) -bed $< | cut -d ':' -f 1 > data/downstream/intergenic/augustus/sequences.fa
	awk '/^>/ {if (out) close(out); gene = substr($$1,2);\
	out="data/downstream/intergenic/augustus/sequences/" gene ".fa"} {print > out}' data/downstream/intergenic/augustus/sequences.fa

# compute augustus predictions for neighborhoods of intergenic regions
data/downstream/intergenic/augustus/predictions/: data/downstream/intergenic/closest.bed data/downstream/intergenic/augustus/sequences/
	mkdir -p $@
	awk -F '\t' '{if ($$6 == "+") {strand = "forward"} else {strand = "backward"}; \
	cmd = "augustus --singlestrand=true --genemodel=partial --strand=" strand \
	" --uniqueGeneId=true --species=human data/downstream/intergenic/augustus/sequences/" $$4 \
	".fa > data/downstream/intergenic/augustus/predictions/" $$4; system(cmd)}' $<

# combine augustus predicted genes into one gtf
data/downstream/intergenic/augustus/predictions.gtf: data/downstream/intergenic/augustus/predictions/ \
      data/downstream/intergenic/extended_intervals.bed
	bash scripts/bash/combine_augustus.sh $^ $@

# update intergenic regions list to have a column showing if they intersect with local augustus predictions
data/downstream/intergenic/augustus.bed: data/downstream/intergenic/closest.bed \
     data/downstream/intergenic/augustus/predictions.gene_ranges_sorted.bed
		bedtools intersect -s -wa -loj -a $< -b data/downstream/intergenic/augustus/predictions.gene_ranges_sorted.bed | \
		cut -f 1-15 | \
		awk -F '\t' 'BEGIN {OFS=FS} !seen[$$4]++ {if ($$15 == ".") {a = 0} else {a = 1}; \
		print $$1,$$2,$$3,$$4,$$5,$$6,$$7,$$8,$$9,$$10,$$11,$$12,$$13,$$14,a}' > $@

# filter intergenic regions to contain only intergenic regions found in sufficient number of samples
data/downstream/intergenic/filtered.bed: data/downstream/intergenic/augustus.bed
	awk '{if (gsub(/lung/, "&") == 4 || gsub(/eye/, "&") == 3 || gsub(/brain/, "&") == 2 || \
	gsub(/PBMC_indrops/, "&") == 2 || gsub(/PBMC_10x/, "&") == 3) {print $$0}}' $< > $@
	
data/downstream/intergenic/isolated/isolated.bed: data/downstream/intergenic/filtered.bed
	mkdir -p data/downstream/intergenic/isolated/
	awk -F '\t' '$$13 != 0' $< > $@
	
data/downstream/intergenic/antisense/antisense.bed: data/downstream/intergenic/filtered.bed
	mkdir -p data/downstream/intergenic/antisense/
	awk -F '\t' '$$13 == 0' $< > $@
	
data/downstream/intergenic/isolated_sequences.tsv: data/downstream/intergenic/isolated.bed
	awk -F '\t' 'BEGIN {OFS = FS} {if($$6 == "+") {print $$1, $$2 - 10000, $$3, $$4} else {print $$1, $$2, $$3 + 10000, $$4}}' $< | \
	bedtools getfasta -nameOnly -tab -fi $(fasta) -bed stdin > $@

data/downstream/intergenic/AT_distances_isolated.tsv: data/downstream/intergenic/isolated_sequences.tsv data/downstream/intergenic/isolated.bed
	python scripts/python/find_ATrich_upstream.py $^ $@
	
data/downstream/intergenic/ATAC_distances_isolated.bed: data/downstream/intergenic/isolated.bed data/genome/ATAC/PBMC/ATAC.bed
	bedtools closest -t first -D a -id -a <(sort -k1,1 -k2,2n $<) \
	-b <( grep -f <( cut -f 1 data/downstream/intergenic/isolated.bed | sort | uniq | awk '{print "^" $$1 "\t"}' ) \
	data/genome/ATAC/PBMC/ATAC.bed) | awk -F '\t' 'BEGIN {OFS = FS} {print $$4, $$22 * -1}' > $@
	
data/downstream/intergenic/gene_distances_upstream_isolated.bed: data/downstream/intergenic/isolated.bed
	bedtools closest -mdb all -t first -D a -id -a <(sort -k1,1 -k2,2n $<) \
	-b data/genome/references/gencode.gene_ranges_sorted.bed \
	data/genome/references/ncbi.gene_ranges_sorted.bed data/genome/references/lnc.gene_ranges_sorted.bed -names gencode ncbi lnc | \
	cut -f 4,16,20,22,27 | awk -F '\t' 'BEGIN{OFS=FS} {print $$1, $$2, $$3, $$4, $$5 * -1}' > $@
	
# rule to find distances from 3' ends of random 10000 genes to open chromatin upstream
data/downstream/intergenic/ATrich_distances_from_random_genes.bed: data/genome/references/gencode.gene_ranges_sorted.bed 
	shuf $< | head -10000 | cut -f 1-6 > $@_tmp
	awk -F '\t' 'BEGIN {OFS = FS} function max(a, b) { return (a > b) ? a : b } {if($$6 == "+") {print $$1, max(0, $$2 - 10000), $$3, $$4} \
	else {print $$1, $$2, $$3 + 10000, $$4}}' $@_tmp | \
	bedtools getfasta -nameOnly -tab -fi $(fasta) -bed stdin > $@_sequences
	python scripts/python/find_ATrich_upstream.py $@_sequences $@_tmp $@
	rm $@_sequences $@_tmp
	
# rule to find distances from 10000 random locations to open chromatin upstream
data/downstream/intergenic/ATrich_distances_from_random_locs.bed:
	bash scripts/bash/generate_random_locations.sh $@_random_loc 
	awk -F '\t' 'BEGIN {OFS = FS} {print $$1, $$2, $$3, "RAND" NR, ".", (rand() < 0.5 ? "-" : "+")}' $@_random_loc > $@_rand_bed
	awk -F '\t' 'BEGIN {OFS = FS} function max(a, b) { return (a > b) ? a : b } {if($$6 == "+") {print $$1, max(0, $$2 - 10000), $$3, $$4} \
	else {print $$1, $$2, $$3 + 10000, $$4}}' $@_rand_bed | \
	bedtools getfasta -nameOnly -tab -fi $(fasta) -bed stdin > $@_sequences
	python scripts/python/find_ATrich_upstream.py $@_sequences $@_rand_bed $@
	rm $@_sequences $@_rand_bed $@_random_loc
	
data/downstream/intergenic/isolated_distance_combined.tsv: data/downstream/intergenic/gene_distances_upstream_isolated.bed \
data/downstream/intergenic/ATAC_distances_isolated.bed \
data/downstream/intergenic/AT_distances_isolated.tsv \
data/downstream/intergenic/isolated.bed
	paste <(sort -k1,1 $<) <(sort -k1,1 data/downstream/intergenic/ATAC_distances_isolated.bed) \
	<(sort -k4,4 data/downstream/intergenic/isolated.bed) | \
	cut -f 1,5,7-13 | awk -F '\t' '$$2 > $$3 + 2000' > $@
	
data/downstream/intergenic/isolated_distance_combined_closest_ATAC.bed: data/downstream/intergenic/isolated_distance_combined.tsv
	cut -f 4- $< | sort -k1,1 -k2,2n | bedtools closest -wa -wb -D a -id -t first -a stdin \
	-b <( grep -f <( cut -f 4 $< | uniq | awk '{print "^" $$1 "\t"}' ) data/genome/ATAC/PBMC/ATAC.bed) > $@
	
data/downstream/intergenic/isolated_distance_combined_closest_ATAC_sequences.tsv: data/downstream/intergenic/isolated_distance_combined_closest_ATAC.bed
	awk -F '\t' 'BEGIN {OFS = FS} {print $$7, $$8, $$9, $$4}' $< | \
	bedtools getfasta -tab -nameOnly -fi $(fasta) -bed stdin > $@
	
data/downstream/intergenic/isolated_distance_combined_closest_ATAC_sequences_ATcheck.tsv: data/downstream/intergenic/isolated_distance_combined_closest_ATAC_sequences.tsv
	python scripts/python/check_if_AT_rich.py $< > $@

data/downstream/intergenic/sequences_downstream.tsv: data/downstream/intergenic/filtered.bed
	bedtools getfasta -nameOnly -tab -fi $(fasta) -bed <(awk -F '\t' 'BEGIN {OFS = FS} {if ($$6 == "+") {print $$1, $$3, $$3 + 100, $$4} \
	else {print $$1, $$2 - 100, $$2, $$4}}' $< ) > $@
	
data/downstream/intergenic/sequences_upstream.tsv: data/downstream/intergenic/filtered.bed
	bedtools getfasta -nameOnly -tab -fi $(fasta) -bed <(awk -F '\t' 'BEGIN {OFS = FS} {if ($$6 == "+") {print $$1, $$2 - 100, $$2, $$4} \
	else {print $$1, $$3, $$3 + 100, $$4}}' $< ) > $@
	
data/downstream/intergenic/sequences.tsv: data/downstream/intergenic/filtered.bed
	bedtools getfasta -nameOnly -tab -fi $(fasta) -bed $< > $@
	
data/downstream/intergenic/sequences_downstream_ATcheck.tsv: data/downstream/intergenic/sequences_downstream.tsv
	python scripts/python/check_if_AT_rich.py $< > $@
	
data/downstream/intergenic/sequences_upstream_ATcheck.tsv: data/downstream/intergenic/sequences_upstream.tsv
	python scripts/python/check_if_AT_rich.py $< > $@
	
data/downstream/intergenic/sequences_ATcheck.tsv: data/downstream/intergenic/sequences.tsv
	python scripts/python/check_if_AT_rich.py $< > $@

data/downstream/intergenic/closest_antisense_not_intersecting.bed: data/downstream/intergenic/filtered.bed
	sort -k1,1 -k2,2n $< | \
	bedtools closest -mdb all -t first -d -S -io -a stdin -b data/genome/references/gencode.gene_ranges_sorted.bed \
	data/genome/references/ncbi.gene_ranges_sorted.bed data/genome/references/lnc.gene_ranges_sorted.bed -names gencode ncbi lnc > $@

data/downstream/intergenic/sequences_extended_with_strands.tsv: data/downstream/intergenic/filtered.bed
	bedtools getfasta -nameOnly -tab -fi $(fasta) \
	-bed <(awk -F '\t' 'BEGIN {OFS=FS} {print $$1, $$2 - 1000, $$3 + 1000, $$4, ".", $$6}' $< | sort -k1,1 -k2,2n) > $@.tmp
	paste <(sort -k1,1 $@.tmp) <(cut -f 4,6 $< | sort -k1,1 | cut -f 2) | awk -F'\t' 'BEGIN {OFS = FS} {print $$1, $$3, $$2}' > $@
	rm $@.tmp
	
	
data/downstream/intergenic/TATA_AT.tsv: data/downstream/intergenic/sequences_extended_with_strands.tsv
	python scripts/python/TATA_AT.py $< $@
	
data/downstream/intergenic/isolated/TATA_AT.tsv: data/downstream/intergenic/TATA_AT.tsv data/downstream/intergenic/isolated/isolated.bed
	awk -F '\t' 'NR==FNR {keys[$$4]; next} $$1 in keys' data/downstream/intergenic/isolated/isolated.bed $< > $@
	
data/downstream/intergenic/antisense/TATA_AT.tsv: data/downstream/intergenic/TATA_AT.tsv data/downstream/intergenic/antisense/antisense.bed
	awk -F '\t' 'NR==FNR {keys[$$4]; next} $$1 in keys' data/downstream/intergenic/antisense/antisense.bed $< > $@
	
data/downstream/intergenic/isolated/ATAC_distances.bed: data/downstream/intergenic/isolated/isolated.bed
	bedtools closest -D a -id -a $< -b data/genome/ATAC/PBMC/ATAC.bed > $@
	
# rule to clean working directory (mainly from the latex intermediates)
clean:
	latexmk -c
	rm -f *.bbl
	rm -f *.bib.bak$
