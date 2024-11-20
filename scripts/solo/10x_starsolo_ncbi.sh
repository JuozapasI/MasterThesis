#!/bin/bash
#SBATCH -J solo_J
#SBATCH -o /tmp/Mazutislab-out/Juozapas/Thesis/data/PBMC_ncbi/slurm.log
#SBATCH --partition Cluster-public
#SBATCH --nodes=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=120G
#SBATCH --time=48:00:00


# activate conda-env, if the version from env needed
# source activate /home/MazutisLab/software/pkg/miniconda3/envs/star-J
# export PATH="/home/MazutisLab/software/pkg/miniconda3/envs/star-J/bin:$PATH"

fastq_dir="/tmp/Mazutislab-out/Juozapas/Thesis/data/PBMC_10x/Chromium_3p_GEX_Human_PBMC_fastqs/"

# run STAR with comperhensive annotations
STAR \
    --genomeDir /tmp/Mazutislab-out/Juozapas/Thesis/data/index_ncbi/ \
    --readFilesIn  "${fastq_dir}Chromium_3p_GEX_Human_PBMC_S1_L001_R2_001.fastq.gz" "${fastq_dir}Chromium_3p_GEX_Human_PBMC_S1_L001_R1_001.fastq.gz" \
    --soloCBwhitelist /tmp/Mazutislab-out/Juozapas/Thesis/data/PBMC_10x/3M-february-2018.txt \
    --runThreadN 14 \
    --outFileNamePrefix /tmp/Mazutislab-out/Juozapas/Thesis/data/PBMC_10x/solo_output_ncbi/ \
    --readFilesCommand zcat \
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
    
