#!/bin/bash
#SBATCH -J solo_J
#SBATCH -o /tmp/Mazutislab-out/Juozapas/Thesis/data/index_10x/slurm.log
#SBATCH --partition Cluster-public
#SBATCH --nodes=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=120G
#SBATCH --time=48:00:00


# run STAR
STAR \
    --runMode genomeGenerate \
    --runThreadN 14 \
    --genomeDir /tmp/Mazutislab-out/Juozapas/Thesis/data/index_10x/ \
    --genomeFastaFiles /tmp/Mazutislab-out/Juozapas/Thesis/data/GRCh38.dna.primary_assembly.fa \
    --sjdbGTFfile /tmp/Mazutislab-out/Juozapas/Thesis/data/genes_10x.gtf

