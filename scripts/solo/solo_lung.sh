#!/bin/bash
#SBATCH -J solo_J
#SBATCH -o /tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/lung/slurm_solo.log
#SBATCH --partition Cluster-public
#SBATCH --nodes=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=160G
#SBATCH --time=48:00:00

source activate /home/MazutisLab/software/pkg/miniconda3/envs/RE
export PATH="/home/MazutisLab/software/pkg/miniconda3/envs/RE/bin:$PATH"

read1_files=$1
read2_files=$2
i=$3

STAR \
	    --genomeDir /tmp/Mazutislab-out/Juozapas/Thesis/data/genome/indices/index_10x/ \
	    --readFilesIn  "$read2_files" "$read1_files" \
	    --soloCBwhitelist /tmp/Mazutislab-out/Juozapas/Thesis/data/barcode_whitelists/10x_v3.1_barcodes/3M-february-2018.txt \
	    --runThreadN 14 \
	    --outFileNamePrefix /tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/lung/solo_output.$i/ \
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
	    --soloUMIlen 12 \
	    --limitBAMsortRAM 126958727012
