#!/bin/bash
#SBATCH -J solo_J
#SBATCH -o /tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/PBMC_10x_3/solo_output.10x/slurm.log
#SBATCH --partition Cluster-public
#SBATCH --nodes=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=120G
#SBATCH --time=48:00:00


# activate conda-env, if the version from env needed
# source activate /home/MazutisLab/software/pkg/miniconda3/envs/star-J
# export PATH="/home/MazutisLab/software/pkg/miniconda3/envs/star-J/bin:$PATH"

fastq_dir="/tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/PBMC_10x_3/fastq/10k_PBMC_3p_nextgem_Chromium_Controller_fastqs/"

read1_files=""
read2_files=""

for lane in {1..4}; do
    # Concatenate read 1 files
    read1_files+="${fastq_dir}10k_PBMC_3p_nextgem_Chromium_Controller_S2_L00${lane}_R1_001.fastq.gz,"
    # Concatenate read 2 files
    read2_files+="${fastq_dir}10k_PBMC_3p_nextgem_Chromium_Controller_S2_L00${lane}_R2_001.fastq.gz,"
done

# Remove the trailing comma from the concatenated file names
read1_files="${read1_files%,}"
read2_files="${read2_files%,}"


# run STAR with comperhensive annotations
STAR \
    --genomeDir /tmp/Mazutislab-out/Juozapas/Thesis/data/genome/indices/index_10x/ \
    --readFilesIn  "$read2_files" "$read1_files" \
    --soloCBwhitelist /tmp/Mazutislab-out/Juozapas/Thesis/data/barcode_whitelists/10x_v3.1_barcodes/3M-february-2018.txt \
    --runThreadN 14 \
    --outFileNamePrefix /tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/PBMC_10x_3/solo_output.10x/ \
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
    
