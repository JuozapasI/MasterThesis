#!/bin/bash
#SBATCH -J solo_J
#SBATCH -o /tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/eye2/slurm_solo.log
#SBATCH --partition Cluster-public
#SBATCH --nodes=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=120G
#SBATCH --time=48:00:00


# activate conda-env, if the version from env needed
source activate /home/MazutisLab/software/pkg/miniconda3/envs/RE
export PATH="/home/MazutisLab/software/pkg/miniconda3/envs/RE/bin:$PATH"

samples_dir="/tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/"

for i in {1..3}; do
	for type in "fovea" "peripheral"; do
		read1_files=""
		read2_files=""
		
		for dir in ${samples_dir}eye2_d$i/fastq/${type}/*; do
			# Find R1 and R2 files and store them in variables
			r1=$(find "$dir" -type f -name "*1.fastq.gz" | sort)
			r2=$(find "$dir" -type f -name "*2.fastq.gz" | sort)

			# Append to arrays
			read1_files+=($r1)
			read2_files+=($r2)
		done
		
		# Merge filenames into strings separated by comma
		read1_files=$(IFS=,; echo "${read1_files[*]}" | cut -c 2-)
		read2_files=$(IFS=,; echo "${read2_files[*]}" | cut -c 2-)
		
		# run STAR with comperhensive annotations
		STAR \
		    --genomeDir /tmp/Mazutislab-out/Juozapas/Thesis/data/genome/indices/index_10x/ \
		    --readFilesIn  "$read2_files" "$read1_files" \
		    --soloCBwhitelist /tmp/Mazutislab-out/Juozapas/Thesis/data/barcode_whitelists/10x_v3.1_barcodes/3M-february-2018.txt \
		    --runThreadN 14 \
		    --outFileNamePrefix /tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/eye2_d$i/solo_output.${type}/ \
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
		    --soloBarcodeReadLength 0
		    
		unset read1_files read2_files
	done
done
