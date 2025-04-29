#!/bin/bash

region=$1
strand=$2
gene=$3

cd /tmp/Mazutislab-out/Juozapas/Thesis/data/imageGeneration/
mkdir -p ${gene}
cd ${gene}

for sample in brain PBMC_10x PBMC_10x_2 PBMC_10x_3 brain_2 lung_8 lung_5 lung_2 lung_7 eye_2 eye_3 PBMC_indrops PBMC_indrops_2 eye; do
	if [[ ! -f "../../datasets/${sample}/solo_output.10x/Aligned.sortedByCoord.out.bam.bai" ]]; then
		samtools index ../../datasets/${sample}/solo_output.10x/Aligned.sortedByCoord.out.bam
	fi
	if [[ "$strand" == "+" ]]; then
	bamCoverage -b ../../datasets/${sample}/solo_output.10x/Aligned.sortedByCoord.out.bam \
	-o ${sample}.bedgraph -of bedgraph -bs 1 -r ${region} --samFlagExclude 16
	else
	bamCoverage -b ../../datasets/${sample}/solo_output.10x/Aligned.sortedByCoord.out.bam \
	-o ${sample}.bedgraph -of bedgraph -bs 1 -r ${region} --samFlagInclude 16
	fi
done
