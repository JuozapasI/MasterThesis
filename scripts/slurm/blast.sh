#!/bin/bash
#SBATCH -J blast-j
#SBATCH -o /tmp/Mazutislab-out/Juozapas/Thesis/slurm_make.log
#SBATCH --partition Cluster-public
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60G
#SBATCH --time=72:00:00

source activate /home/MazutisLab/software/pkg/miniconda3/envs/RE
export PATH="/home/MazutisLab/software/pkg/miniconda3/envs/RE/bin:$PATH"

sample=$1

cd /tmp/Mazutislab-out/Juozapas/Thesis/
make data/blast/viral/${sample}.txt
