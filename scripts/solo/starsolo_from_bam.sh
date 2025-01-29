#!/bin/bash

bam=$1
index=$2
outdir=$3
umi=$4

# Shuffle bam first for more efficient mapping
samtools collate -u -o ${bam}.shuffled.bam $1

STAR \
    --genomeDir $index \
    --readFilesIn  ${bam}.shuffled.bam --readFilesType SAM SE \
    --soloInputSAMattrBarcodeSeq CR UR \
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
    --soloCBtype String \
    --soloCBwhitelist None \
    --soloBarcodeReadLength 0 \
    --soloUMIlen ${umi}
    
rm ${bam}.shuffled.bam
    
