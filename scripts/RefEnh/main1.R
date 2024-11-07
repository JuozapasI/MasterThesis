library(rtracklayer)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 14) {
    stop(paste("Error: Expected 14 arguments, but got", length(args)))
}

workingDirectory = args[1]
original_gtf_file = args[2]
bam_file = args[3]
index_file = args[4]
fasta_file = args[5]
genome_index = args[6]
chr_lengths = args[7]
ends_dist = as.integer(args[8])
seq_depth = as.integer(args[9])
extensionThreshold = as.integer(args[10])
cpmThreshold = as.integer(args[11])
barc_length = as.integer(args[12])
umi_length = as.integer(args[13])
max_read_length = as.integer(args[14])

script_location = "/tmp/Mazutislab-out/Juozapas/Thesis/scripts/RefEnh/"
gtf_file = "unique_names.gtf"

setwd(workingDirectory)
options(scipen = 999)

source(paste0(script_location, "0.supplementaryFunctions.R"))
source(paste0(script_location, "1.uniqueNames.R"))
source(paste0(script_location, "2.geneLocationBED.R"))
source(paste0(script_location, "3.intergenicReads.R"))
source(paste0(script_location, "4.clusterIntergenicReads.R"))
source(paste0(script_location, "5.extensionCandidates.R"))
source(paste0(script_location, "6.identifyOverlappingGenes.R"))
source(paste0(script_location, "7.resolveOverlaps.R"))

print("Please do manual curration and then run main2.R")
