## EPITOPOLOGY
## GOAL: label epitopes on the viral topologyace based on topology information.

library(dplyr, quietly = T, warn.conflicts = F)
library(tidyr, quietly = T, warn.conflicts = F)
library(stringr, quietly = T, warn.conflicts = F)
library(argparser, quietly = T, warn.conflicts = F)

### SCRIPT ARGUMENTS
# Create a parser
p <- arg_parser("EPITOPOLOGY from MANUAL")

# Add command line arguments
p <- add_argument(p, "--start_pos", help= "Numbers sep by commas with the initial positions of the extracellular domains", type="character", default = "1")
p <- add_argument(p, "--end_pos", help= "Numbers sep by commas with the initial positions of the extracellular domains", type="character", default = "5")
p <- add_argument(p, "--input_epitopes", help= "Path to epitope dataset resulting from Epimerger.R (.csv format)", type="character", default = ".")
p <- add_argument(p, "--sample", help = "Sample name to label output files", type = "character", default = "topology_extracted_manual.csv")
p <- add_argument(p, "--outdir", help = "Path to output files", type = "character", default = "E_epitopology")
#p <- add_argument(p, "--save_rdata_dir", help = "Path to save rData image", type = "character", default = "E_epitoplogy/")

# Parse the command line arguments
argv <- parse_args(p)

## MANUAL DOMAIN UPLOAD (info extracted from Uniprot)
## ADD ONLY EXTRAVIRAL DOMAINS
topology <- data.frame(extracellular.start = c(argv$start_positions), extracellular.end = c(argv$end_positions))
#topology

## CHARACTER TO NUMERIC
topology$extracellular.start <- as.numeric(topology$extracellular.start)
topology$extracellular.end <- as.numeric(topology$extracellular.end)

### IMPORT DATAFRAME from EPIMERGER.PY
merged <- read.csv(file = argv$input_epitopes, header = T, sep = ";")

## LABEL EXTRACELLULAR DOMAINS (each extracellular correspond to one extracellular domain. Add or remove according to the extraviral domains of the protein.)
# extracellular1 <- c()
# for (x in 1:length(merged$Sequence)){
#   y <- ifelse(merged$Start[x] >= topology$extracellular.start[1] & merged$End[x] <= topology$extracellular.end[1], "Extracellular", "Non-extracellular")
#   extracellular1 <- c(extracellular1, y)
#   print(extracellular1)
# }

extracellular_epitope <- c()
extracellular_label <- c()
for (x in 1:length(merged$Sequence)){
  for (z in 1:length(topology$extracellular.start)){
    y <- ifelse(merged$Start[x] >= topology$extracellular.start[z] & merged$End[x] <= topology$extracellular.end[z], "Extracellular", "Non-extracel")
    extracellular_epitope <- c(extracellular_epitope, y)
  }
  extracellular_label <- c(extracellular_label, extracellular_epitope)
  extracellular_epitope <- c()
  #print(extracellular_label)
}

### LABEL EPITOPES AS EXTRACELLULAR
merged <- cbind(merged, extracellular_label)

## EXPORT RESULTS
write.table(merged, sep = ";", file = paste0(argv$outdir, "/", argv$sample, sep = ""), quote = F, row.names = F)
print(paste("Find your output file at: ", argv$outdir, "/", argv$sample, sep = ""))

### EXPORT RDATA
#p <- add_argument(p, "--save_rdata_dir", help = "Path to save rData image", type = "character", default = ".")
#dir.create(paste0(argv$save_rdata_dir, "/rdata", sep = ""))
#save.image(paste0(argv$save_rdata_dir, "/rdata/", argv$sample, ".rdata", sep = ""))
