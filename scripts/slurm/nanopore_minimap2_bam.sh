#!/bin/bash
#SBATCH -J minimap2
#SBATCH -o /tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/PBMC_10x/slurm_minimap2_bam.log
#SBATCH --partition Cluster-public
#SBATCH --nodes=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=120G
#SBATCH --time=48:00:00


# activate conda-env, if the version from env needed
source activate /home/MazutisLab/software/pkg/miniconda3/envs/explore-J
export PATH="/home/MazutisLab/software/pkg/miniconda3/envs/explore-J/bin:$PATH"

cd /tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/PBMC_10x/nanopore/

ref="/tmp/Mazutislab-out/Juozapas/Thesis/data/genome/fasta/GRCh38.dna.primary_assembly.fa"

minimap2 -t 14 -ax splice ${ref} SC3pv3_GEX_Human_PBMC_ONT.fastq.gz | samtools view -@ 4 -b -o aln.bam
