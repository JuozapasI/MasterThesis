#!/bin/bash
#SBATCH -J solo_J
#SBATCH -o /tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/lung/slurm_solo.log
#SBATCH --partition Cluster-public
#SBATCH --nodes=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=120G
#SBATCH --time=48:00:00


# activate conda-env, if the version from env needed
#source activate /home/MazutisLab/software/pkg/miniconda3/envs/RE
#export PATH="/home/MazutisLab/software/pkg/miniconda3/envs/RE/bin:$PATH"

fastq_dir="/tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/lung/fastq/"

for i in 2 5 7 8; do

	read1_files=""
	read2_files=""
	
	for dir in ${fastq_dir}$i/*; do
		# Find R1 and R2 files and store them in variables
		r1=$(find "$dir" -type f -name "merged_R1_shufled_shufled.fastq.gz" | sort)
		r2=$(find "$dir" -type f -name "merged_R2_shufled_shufled.fastq.gz" | sort)

		# Append to arrays
		read1_files+=($r1)
		read2_files+=($r2)
	done
	
	# Merge filenames into strings separated by comma
	read1_files=$(IFS=,; echo "${read1_files[*]}" | cut -c 2-)
	read2_files=$(IFS=,; echo "${read2_files[*]}" | cut -c 2-)
	
	# run STAR with comperhensive annotations
	sbatch --output /tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/lung/solo_output_${i}.log \
	/tmp/Mazutislab-out/Juozapas/Thesis/scripts/solo/solo_lung.sh $read1_files $read2_files $i
	    
	unset read1_files read2_files
done
