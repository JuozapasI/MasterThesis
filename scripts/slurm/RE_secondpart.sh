#!/bin/bash
#SBATCH -J RE
#SBATCH -o /tmp/Mazutislab-out/Juozapas/REv3/slurm.log
#SBATCH --partition Cluster-public
#SBATCH --nodes=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=120G
#SBATCH --time=48:00:00


# activate conda-env, if the version from env needed
source activate /home/MazutisLab/software/pkg/miniconda3/envs/RefEnh
export PATH="/home/MazutisLab/software/pkg/miniconda3/envs/star-J/bin:$PATH"
    
Rscript /tmp/Mazutislab-out/Juozapas/Thesis/scripts/RefEnh/main2.R \
  /tmp/Mazutislab-out/Juozapas/Thesis/data/PBMC_10x/RefEnh/ \
  /tmp/Mazutislab-out/Juozapas/Thesis/data/chr_lenghts.tsv 


