#!/bin/bash
#SBATCH -J solo_J
#SBATCH -o /tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/PBMC_indrops_2/solo_output.10x/slurm.log
#SBATCH --partition Cluster-public
#SBATCH --nodes=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=120G
#SBATCH --time=48:00:00


# activate conda-env, if the version from env needed
source activate /home/MazutisLab/software/pkg/miniconda3/envs/RE
export PATH="/home/MazutisLab/software/pkg/miniconda3/envs/RE/bin:$PATH"

fastq_dir="/tmp/simonasj-out/samples/indrop_tso_ivt/raw/PBMC/PBMC_V4_fresh_InDrop2_rep/"
cdna="${fastq_dir}genomic/PBMC_fresh_InDrop2_rep_IGO_11244_7_S2_L001_R2_001.fastq.gz,${fastq_dir}genomic/PBMC_fresh_InDrop2_rep_IGO_11244_7_S2_L002_R2_001.fastq.gz"
bc="${fastq_dir}barcode/PBMC_fresh_InDrop2_rep_IGO_11244_7_S2_L001_R1_001.fastq.gz,${fastq_dir}barcode/PBMC_fresh_InDrop2_rep_IGO_11244_7_S2_L002_R1_001.fastq.gz"

# run STAR with comperhensive annotations
STAR \
    --genomeDir /tmp/Mazutislab-out/Juozapas/Thesis/data/genome/indices/index_10x/ \
    --readFilesIn  ${cdna} ${bc} \
    --soloCBwhitelist  /home/simonasj/references/barcodes/in_drop_v4/flat/v4_cb1.txt /home/simonasj/references/barcodes/in_drop_v4/flat/v4_cb2.txt \
    --runThreadN 14 \
    --outFileNamePrefix /tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/PBMC_indrops_2/solo_output.10x/ \
    --readFilesCommand zcat \
    --runDirPerm All_RWX \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS CR UR CB UB sS sQ sM GX GN \
    --outSAMunmapped Within \
    --soloMultiMappers PropUnique \
    --soloFeatures GeneFull \
    --soloType CB_UMI_Complex \
    --soloCBposition 0_0_0_7 0_12_0_19  \
    --soloUMIposition 0_20_0_27 \
    --soloCBmatchWLtype 1MM \
    --soloUMIdedup Exact 
    
