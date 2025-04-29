#!/bin/bash
#SBATCH -J make-j
#SBATCH -o /tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/brain_2/slurm_make.log
#SBATCH --partition Cluster-public
#SBATCH --nodes=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=120G
#SBATCH --time=72:00:00

source activate /home/MazutisLab/software/pkg/miniconda3/envs/RE
export PATH="/home/MazutisLab/software/pkg/miniconda3/envs/RE/bin:$PATH"

sample=$1
umi=$2
barcode=$3

cd /tmp/Mazutislab-out/Juozapas/Thesis/
make ${sample}.dataset umi=${umi} barcode=${barcode}
