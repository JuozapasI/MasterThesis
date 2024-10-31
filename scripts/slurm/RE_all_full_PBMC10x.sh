#!/bin/bash
#SBATCH -J RE
#SBATCH -o /tmp/Mazutislab-out/Juozapas/Thesis/data/PBMC_10x/RefEnh_full_reference/slurm.log
#SBATCH --partition Cluster-public
#SBATCH --nodes=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=120G
#SBATCH --time=48:00:00


# activate conda-env, if the version from env needed
source activate /home/MazutisLab/software/pkg/miniconda3/envs/RefEnh
export PATH="/home/MazutisLab/software/pkg/miniconda3/envs/star-J/bin:$PATH"
    
Rscript /tmp/Mazutislab-out/Juozapas/Thesis/scripts/RefEnh/main1.R \
  /tmp/Mazutislab-out/Juozapas/Thesis/data/PBMC_10x/RefEnh_full_reference/ \
  /tmp/Mazutislab-out/Juozapas/Thesis/data/gencode.v47.gtf \
  /tmp/Mazutislab-out/Juozapas/Thesis/data/PBMC_10x/solo_output_full_reference/Aligned.sortedByCoord.out.bam \
  /tmp/Mazutislab-out/Juozapas/Thesis/data/PBMC_10x/solo_output_full_reference/Aligned.sortedByCoord.out.bam.bai \
  /tmp/Mazutislab-out/Juozapas/Thesis/data/GRCh38.dna.primary_assembly.fa \
  /tmp/Mazutislab-out/Juozapas/Thesis/data/GRCh38.dna.primary_assembly.fa.fai \
  /tmp/Mazutislab-out/Juozapas/Thesis/data/chr_lenghts.tsv \
  10000 \
  182330834 \
  1000 \
  5

Rscript /tmp/Mazutislab-out/Juozapas/Thesis/scripts/RefEnh/main2.R \
  /tmp/Mazutislab-out/Juozapas/Thesis/data/PBMC_10x/RefEnh_full_reference/ \
  /tmp/Mazutislab-out/Juozapas/Thesis/data/chr_lenghts.tsv 

