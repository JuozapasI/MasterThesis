#!/bin/bash
#SBATCH -J solo_J
#SBATCH -o /tmp/Mazutislab-out/Juozapas/Thesis/data/PBMC_indrops/solo_output.10x/slurm.log
#SBATCH --partition Cluster-public
#SBATCH --nodes=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=120G
#SBATCH --time=48:00:00


# activate conda-env, if the version from env needed
# source activate /home/MazutisLab/software/pkg/miniconda3/envs/star-J
# export PATH="/home/MazutisLab/software/pkg/miniconda3/envs/star-J/bin:$PATH"

fastq_dir="/tmp/simonasj-out/samples/indrop_tso_ivt/raw/Rev/241023_VH00558_635_AAG32JKM5_out/24_BEAD_KG/merged_fastq/"
cdna="${fastq_dir}24_BEAD_KG_06_S4_cdna_001.fastq.gz"
bc="${fastq_dir}24_BEAD_KG_06_S4_bc_001.fastq.gz"
whitelist_dir="/tmp/Mazutislab-out/Juozapas/Thesis/data/indrops_barcodes/"

# run STAR with comperhensive annotations
STAR \
    --genomeDir /tmp/Mazutislab-out/Juozapas/Thesis/data/index_10x/ \
    --readFilesIn  ${cdna} ${bc} \
    --soloCBwhitelist  ${whitelist_dir}bc{1,2,3}_list.txt \
    --runThreadN 14 \
    --outFileNamePrefix /tmp/Mazutislab-out/Juozapas/Thesis/data/PBMC_indrops/solo_output.10x/ \
    --readFilesCommand zcat \
    --runDirPerm All_RWX \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS CR UR CB UB sS sQ sM GX GN \
    --outSAMunmapped Within \
    --soloMultiMappers PropUnique \
    --soloFeatures GeneFull \
    --soloType CB_UMI_Complex \
    --soloCBposition 0_0_0_7 0_8_0_17 0_18_0_25 \
    --soloUMIposition 0_26_0_35 \
    --soloCBmatchWLtype 1MM \
    --soloUMIdedup Exact 
    
