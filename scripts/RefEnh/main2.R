library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)

workingDirectory = args[1]
chr_lengths = args[2]


gtf_file = "unique_names.gtf"

setwd(workingDirectory)
options(scipen = 999)

source("0.supplementaryFunctions.R")
source("8.finalGTF.R")

print("Everything is done.")