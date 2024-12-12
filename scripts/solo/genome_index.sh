#!/bin/bash

gtf=$1
fasta=$2
directory=$3

# run STAR
STAR \
    --runMode genomeGenerate \
    --runThreadN 14 \
    --genomeDir ${directory} \
    --genomeFastaFiles ${fasta} \
    --sjdbGTFfile ${gtf}

