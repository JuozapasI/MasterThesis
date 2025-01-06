umi = 12
barcode = 16
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

PBMC_10x_seq_depth = 182330834
PBMC_10x_2_seq_depth = 496387931
brain_seq_depth = 206360627
PBMC_indrops_seq_depth = 112932507

PBMC_10x_extensionCountThreshold = 182
PBMC_10x_2_extensionCountThreshold = 496
brain_extensionCountThreshold = 206
PBMC_indrops_extensionCountThreshold = 112

PBMC_10x_newGeneCountThreshold = 911
PBMC_10x_2_newGeneCountThreshold = 2481
brain_newGeneCountThreshold = 1031
PBMC_indrops_newGeneCountThreshold = 564

debug:
	echo $(references_paths)

.SECONDARY:

#.PHONY: $(datasets)

all: $(datasets)

%.dataset: data/downstream/summaries/count_summaries/%.csv \
     data/downstream/summaries/gene_summaries/%.csv \
     data/downstream/matrices/%/10x/matrix.mtx.gz \
     data/downstream/matrices/%/10x/features.tsv.gz \
     data/downstream/matrices/%/10x/barcodes.tsv.gz \
     data/downstream/matrices/%/final/matrix.mtx.gz \
     data/downstream/matrices/%/final/features.tsv.gz \
     data/downstream/matrices/%/final/barcodes.tsv.gz \
     data/downstream/intergenic/%.csv
	@echo "Building $* dataset"

%.pdf: %.tex
	latexmk -pdf -silent -deps-out=.depend $*

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
	gtf2bed | \
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


%.clusters.bed: %.bam
	bedtools bamtobed -i $< | \
	sort -k1,1 -k2,2n | \
	bedtools merge -s -c 6 -o distinct,count | \
	sort -k1,1 -k2,2n | \
	awk -F'\t' 'BEGIN {OFS = FS} {print $$1, $$2, $$3, ".", $$5, $$4}' > $@
	
%.clusters.thrash.bed: %.clusters.bed
	awk -F'\t' '{if ($$5 <= $(clusterThreshold)) {print $$0}}' $< > $@
	
%.clusters.good.bed: %.clusters.bed
	awk -F'\t' '{if ($$5 > $(clusterThreshold)) {print $$0}}' $< > $@

$(foreach f, $(order), $(eval \
data/datasets/%/unassigned_reads/$(f).intersecting_genes.bed: data/datasets/%/unassigned_reads/$(f).intersecting.clusters.good.fine.bed \
data/genome/references/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gene_ranges_sorted.bed; \
bedtools intersect -s -wa -wb -abam $$< -b data/genome/references/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gene_ranges_sorted.bed | \
sort -k1,1 -k2,2n > $$@))

%.intersecting_gene_list.tsv: %.intersecting_genes.bed
	cut -f 11 $< | sort -u > $@
	
$(foreach f, $(order), $(eval \
data/datasets/%/unassigned_reads/$(f).overlappers.csv: data/genome/references/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gtf \
data/datasets/%/unassigned_reads/$(f).intersecting_gene_list.tsv; \
Rscript scripts/R/Overlappers.R $$^ $$@))

$(foreach f, $(order), $(eval \
data/datasets/%/unassigned_reads/$(f).overlaps_modified.gtf: data/genome/references/$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f))).gtf \
data/datasets/%/unassigned_reads/$(f).overlappers.csv; \
Rscript scripts/R/ResolveOverlappers.R $$^ $(ends_dist) $$@))

# Assemble modified gtf for the first reference (i.e. with resolved overlaps)
data/datasets/%/unassigned_reads/$(first).modified.gtf: data/datasets/%/unassigned_reads/$(first).intersecting_gene_list.tsv \
data/genome/references/$(first).gtf data/datasets/%/unassigned_reads/$(first).overlaps_modified.gtf
	( grep -v -f $< data/genome/references/$(first).gtf; cat data/datasets/$*/unassigned_reads/$(first).overlaps_modified.gtf ; ) | \
	sort -k1,1 -k4,4n > $@

# For the second and on, we need to make sure that appended entries doesn't overlap with previous gtf
$(foreach f, $(order_without_first), $(eval \
data/datasets/%/unassigned_reads/$(f).new_entries.gtf: data/datasets/%/unassigned_reads/$(basename $(f)).modified.gtf \
data/datasets/%/unassigned_reads/$(f).overlaps_modified.gtf; \
Rscript scripts/R/ResolveGtfsOverlaps.R $$^ $$@))

# And then we can merge gtfs
$(foreach f, $(order_without_first), $(eval \
data/datasets/%/unassigned_reads/$(f).modified.gtf: data/datasets/%/unassigned_reads/$(basename $(f)).modified.gtf \
data/datasets/%/unassigned_reads/$(f).new_entries.gtf; \
( cat $$< ; cat data/datasets/$$*/unassigned_reads/$(f).new_entries.gtf ; ) | \
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
	

data/datasets/%/unassigned_reads/$(order_last).intergenic.clusters.sequences.tsv: data/datasets/%/unassigned_reads/$(order_last).intergenic.clusters.good.fine.bed
	bedtools getfasta -fi $(fasta) -nameOnly -tab -bed $< | cut -f2 > $@
	
data/datasets/%/unassigned_reads/$(order_last).intergenic.clusters.GCcontent.tsv: data/datasets/%/unassigned_reads/$(order_last).intergenic.clusters.sequences.tsv 
	awk '{l=length($0); gc=gsub(/[GCgc]/,""); print (gc/l)*100}' $< > $@
	
$(foreach f, $(order_without_first), $(eval \
data/datasets/%/solo_output.$(f)/Aligned.sortedByCoord.out.bam: data/datasets/%/unassigned_reads/$(basename $(f)).intergenic.bam \
data/genome/indices/index_$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f)))/; \
bash scripts/solo/starsolo_from_bam.sh $$< data/genome/indices/index_$(word $(words $(subst ., ,$(f))), $(subst ., ,$(f)))/ data/datasets/$$*/solo_output.$(f)/ $(umi)))
	
# Filter those clusters that came from AT rich regions:
%.filteredAT.intergenic.clusters.good.fine.bed: %.intergenic.clusters.good.fine.bed %.intergenic.clusters.GCcontent.tsv
	paste -d ',' $^ | awk -F',' '$$2 > 30 {print $$1}' > $@

#In final part, check if there are any clusters very close to genes 3' ends and extend those genes if there are
%.distances.clusters.good.fine.tsv: %.filteredAT.intergenic.clusters.good.fine.bed %.modified.gene_ranges_sorted.bed
	bedtools closest -a $< -b $*.modified.gene_ranges_sorted.bed -s -D a -id -fu > $@

#Check if everything is ok here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#After changing bed format, if required column is still $18 etc.
# Split clusters into close to 3' ends and not
data/datasets/%/unassigned_reads/final.intergenic.clusters_near_genes.txt: data/datasets/%/unassigned_reads/$(order_last).distances.clusters.good.fine.tsv
	awk -F'\t' '{if((-1 * $$18) < $(closeEndThreshold) && $$11 != ".") {print $$0}}' $< > $@

data/datasets/%/unassigned_reads/final.intergenic.clusters_far_from_genes.txt: data/datasets/%/unassigned_reads/$(order_last).distances.clusters.good.fine.tsv
	awk -F'\t' '{if((-1 * $$18) >= $(closeEndThreshold)) {print $$0}}' $< > $@

# Make gene extension list
data/datasets/%/unassigned_reads/final.extension_candidates.csv: data/datasets/%/unassigned_reads/final.intergenic.clusters_near_genes.txt
	awk 'BEGIN {FS = "\t"; OFS = FS;} {if ($$7 >= $($*_extensionCountThreshold)) {print $$1, $$2, $$3, $$6, $$11}}' $< > $@

# Make new intergenic gene list:
data/datasets/%/unassigned_reads/final.new_gene_list.bed: data/datasets/%/unassigned_reads/final.intergenic.clusters_far_from_genes.txt
	awk -v depth="$($*_seq_depth)" -v threshold="$($*_newGeneCountThreshold)" 'BEGIN {FS = "\t"; OFS = FS;} \
	{if ($$7 >= threshold) {i += 1; print $$1, $$2, $$3, "INTERGENIC" i, $$7 * 1000000 / depth, $$6}}' $< > $@

# Assemble final gtf
%/final.gtf: %/$(order_last).modified.gtf %/final.extension_candidates.csv %/final.new_gene_list.bed
	Rscript scripts/R/final_gtf.R $^ $@

# Make final genome index
data/datasets/%/index_final/: data/datasets/%/unassigned_reads/final.gtf 
	bash scripts/solo/genome_index.sh $< $(fasta) $@
	
$(foreach reference, $(references), $(eval data/genome/indices/index_$(reference)/: data/genome/references/$(reference).gtf ; bash scripts/solo/genome_index.sh $$< $(fasta) $$@ ))


# Run starsolo on final gtf
data/datasets/%/solo_output.final/Aligned.sortedByCoord.out.bam: data/datasets/%/solo_output.$(first)/Aligned.sortedByCoord.out.bam \
data/datasets/%/index_final/
	bash scripts/solo/starsolo_from_bam.sh $^ data/datasets/$*/solo_output.final/ $(umi)


# Make igv snapshots script for extension candidates
data/datasets/%/unassigned_reads/igv_snapshots/extension_candidates.batch.txt: data/datasets/%/unassigned_reads/final.extension_candidates.csv \
data/datasets/%/solo_output.10x/Aligned.sortedByCoord.out.bam
	mkdir -p data/datasets/$*/unassigned_reads/igv_snapshots/extension_candidates/ 
	python scripts/python/igv_batch_script_generator.py $^ \
		  data/datasets/$*/unassigned_reads/igv_snapshots/extension_candidates/ GRCh38.dna.primary_assembly.fa GRCh38.dna.primary_assembly.fa \
		  data/gencode.v47.sorted.gtf $@

# Make igv snaphots script for new gene candidates
data/datasets/%/unassigned_reads/igv_snapshots/new_gene_candidates.batch.txt: data/datasets/%/unassigned_reads/final.new_gene_list.bed \
data/datasets/%/solo_output.10x/Aligned.sortedByCoord.out.bam
	mkdir -p data/datasets/$*/unassigned_reads/igv_snapshots/new_gene_candidates/ 
	python scripts/python/igv_batch_script_generator.py $^ \
		  data/datasets/$*/unassigned_reads/igv_snapshots/new_gene_candidates/ GRCh38.dna.primary_assembly.fa GRCh38.dna.primary_assembly.fa \
		  data/gencode.v47.sorted.gtf $@
		  
# Make igv snaphots script for overlaps (first reference)
data/datasets/%/unassigned_reads/igv_snapshots/overlaps.batch.txt: data/datasets/%/unassigned_reads/$(first).intersecting.clusters.good.fine.bed \
data/datasets/%/solo_output.10x/Aligned.sortedByCoord.out.bam
	mkdir -p data/datasets/$*/unassigned_reads/igv_snapshots/overlaps/ 
	python scripts/python/igv_batch_script_generator.py $^ \
		  data/datasets/$*/unassigned_reads/igv_snapshots/overlaps/ GRCh38.dna.primary_assembly.fa GRCh38.dna.primary_assembly.fa \
		  data/gencode.v47.sorted.gtf $@
		  
		  
		  

# Summary
data/datasets/%/unassigned_reads/Summary.txt: data/datasets/%/solo_output.10x/Aligned.sortedByCoord.out.bam \
data/datasets/%/unassigned_reads/10x.unassigned.bam \
data/datasets/%/unassigned_reads/10x.intersecting.bam \
data/datasets/%/unassigned_reads/10x.intergenic.bam \
data/datasets/%/unassigned_reads/AT_seq.bam \
data/datasets/%/unassigned_reads/10x.gencode.ncbi.lnc.intergenic.clusters.good.bed \
data/datasets/%/unassigned_reads/10x.gencode.ncbi.lnc.intergenic.clusters.thrash.bed \
data/datasets/%/unassigned_reads/10x.gencode.unassigned.bam \
data/datasets/%/unassigned_reads/10x.gencode.intersecting.bam \
data/datasets/%/unassigned_reads/10x.gencode.intergenic.bam \
data/datasets/%/unassigned_reads/10x.gencode.ncbi.unassigned.bam \
data/datasets/%/unassigned_reads/10x.gencode.ncbi.intersecting.bam \
data/datasets/%/unassigned_reads/10x.gencode.ncbi.intergenic.bam \
data/datasets/%/unassigned_reads/10x.gencode.ncbi.lnc.unassigned.bam \
data/datasets/%/unassigned_reads/10x.gencode.ncbi.lnc.intersecting.bam \
data/datasets/%/unassigned_reads/10x.gencode.ncbi.lnc.intergenic.bam 
	bash scripts/bash/statistics.sh $^ > $@

# Computing conservation scores for intergenic genes
%/conservation_scores.csv: %/final.new_gene_list.bed data/genome/conservation/hg38.phastCons100way.bed
	bedtools map -a $< -b data/genome/conservation/hg38.phastCons100way.bed -c 5 -o mean > $@

# Captured gene list:
data/datasets/%/unassigned_reads/captured_genes_$(first).csv: data/datasets/%/solo_output.$(first)/Aligned.sortedByCoord.out.bam
	samtools view $< | grep -o -P "GN:Z:[^ \t]*" | cut -d':' -f3 | sort | uniq > $@
	
data/datasets/%/unassigned_reads/captured_genes_final.csv: data/datasets/%/solo_output.final/Aligned.sortedByCoord.out.bam
	samtools view $< | grep -o -P "GN:Z:[^ \t]*" | cut -d':' -f3 | sort | uniq > $@


data/datasets/%/unassigned_reads/additional_gene_summary.txt: data/datasets/%/unassigned_reads/$(first).intergenic.clusters.good.bed $(references_paths)
	bedtools intersect -wb -s -a $< -b $(references_paths) -names $(references) | \
	awk -F '\t' '{split($$17, a, "gene_type \""); split(a[2], b, "\""); print $$7, b[1];}' | \
	sort | uniq -c | sort -k2,2 -k1,1nr > $@




# data arrangement for downstream analyses

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
	

		
clean:
	latexmk -c
	rm *.bbl
	rm *.bib.bak
