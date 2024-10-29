library(rtracklayer)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

workingDirectory = args[1]
original_gtf_file = args[2]
bam_file = args[3]
index_file = args[4]
fasta_file = args[5]
genome_index = args[6]
chr_lengths = args[7]
end_dist = as.integer(args[8])
seq_depth = as.integer(args[9])
extensionThreshold = as.integer(args[10])
cpmThreshold = as.integer(args[11])
barc_umi_length = as.integer(args[12])

gtf_file = "unique_names.gtf"

setwd(workingDirectory)
options(scipen = 999)

source("0.supplementaryFunctions.R")
source("1.uniqueNames.R")
source("2.geneLocationBED.R")
source("3.intergenicReads.R")
source("4.extensionCandidates.R")
source("5.clusterIntergenicReads.R")
source("6.identifyOverlappingGenes.R")
source("7.resolveOverlaps.R")

print("Please do manual curration and then run main2.R")
