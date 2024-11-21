#!/bin/bash
#SBATCH -J make-j
#SBATCH -o /tmp/Mazutislab-out/Juozapas/Thesis/slurm_make.log
#SBATCH --partition Cluster-public
#SBATCH --nodes=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=120G
#SBATCH --time=48:00:00

source activate /home/MazutisLab/software/pkg/miniconda3/envs/RefEnh
export PATH="/home/MazutisLab/software/pkg/miniconda3/envs/RefEnh/bin:$PATH"

cd /tmp/Mazutislab-out/Juozapas/Thesis/
make data/PBMC_10x/solo_intergenic_gencode/Aligned.sortedByCoord.out.bam

