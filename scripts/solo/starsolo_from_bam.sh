#!/bin/bash

bam=$1
index=$2
outdir=$3

STAR \
    --genomeDir $index \
    --readFilesIn  $bam --readFilesType SAM SE \
    --soloCBwhitelist /tmp/Mazutislab-out/Juozapas/Thesis/data/PBMC_10x/3M-february-2018.txt \
    --runThreadN 14 \
    --outFileNamePrefix $outdir \
    --readFilesCommand samtools view -F 0x100 \
    --readFilesSAMattrKeep None \
    --runDirPerm All_RWX \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS CR UR CB UB sS sQ sM GX GN \
    --outSAMunmapped Within \
    --soloMultiMappers PropUnique \
    --soloFeatures GeneFull \
    --soloType CB_UMI_Simple \
    --soloCBmatchWLtype 1MM \
    --soloUMIdedup Exact \
    --soloUMIlen 12
    
