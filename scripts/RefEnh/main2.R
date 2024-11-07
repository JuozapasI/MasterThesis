library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop(paste("Error: Expected 2 arguments, but got", length(args)))
}

workingDirectory = args[1]
chr_lengths = args[2]


gtf_file = "unique_names.gtf"
script_location = "/tmp/Mazutislab-out/Juozapas/Thesis/scripts/RefEnh/"

setwd(workingDirectory)
options(scipen = 999)

source(paste0(script_location, "0.supplementaryFunctions.R"))
source(paste0(script_location, "8.finalGTF.R"))

print("Everything is done.")
