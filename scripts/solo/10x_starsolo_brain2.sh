#!/bin/bash
#SBATCH -J solo_J
#SBATCH -o /tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/brain_2/solo_output.10x/slurm.log
#SBATCH --partition Cluster-public
#SBATCH --nodes=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=120G
#SBATCH --time=48:00:00


# activate conda-env, if the version from env needed
source activate /home/MazutisLab/software/pkg/miniconda3/envs/RE
export PATH="/home/MazutisLab/software/pkg/miniconda3/envs/RE/bin:$PATH"

fastq_dir="/tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/brain_2/fastq/"

r1="${fastq_dir}10X146-2_S12_L001/10X146-2_S12_L001_R1_001.fastq.gz,${fastq_dir}10X146-2_S12_L002/10X146-2_S12_L002_R1_001.fastq.gz"
r2="${fastq_dir}10X146-2_S12_L001/10X146-2_S12_L001_R2_001.fastq.gz,${fastq_dir}10X146-2_S12_L002/10X146-2_S12_L002_R2_001.fastq.gz"

# run STAR with comperhensive annotations
STAR \
    --genomeDir /tmp/Mazutislab-out/Juozapas/Thesis/data/genome/indices/index_10x/ \
    --readFilesIn  ${r2} ${r1} \
    --soloCBwhitelist /tmp/Mazutislab-out/Juozapas/Thesis/data/barcode_whitelists/10x_v3.1_barcodes/3M-february-2018.txt \
    --runThreadN 14 \
    --outFileNamePrefix /tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/brain_2/solo_output.10x/ \
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
    
