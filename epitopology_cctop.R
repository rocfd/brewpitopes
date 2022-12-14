## EPITOPOLOGY
## GOAL: label epitopes on the viral topologyace based on topology information.

library(dplyr, quietly = T, warn.conflicts = F)
library(tidyr, quietly = T, warn.conflicts = F)
library(stringr, quietly = T, warn.conflicts = F)
library(argparser, quietly = T, warn.conflicts = F)

### SCRIPT ARGUMENTS
# Create a parser
p <- arg_parser("EPITOPOLOGY from CCTOP")

# Add command line arguments
p <- add_argument(p, "--input_CCTOP", help= "Path to CCTOP output dataset from XML parser (.csv format)", type="character", default = "E_epitopology/CCTOP/cctop_domains.csv")
p <- add_argument(p, "--input_epitopes", help= "Path to epitope dataset resulting from Epimerger.py (.csv format)", type="character", default = "D_epimerger/merged.csv")
p <- add_argument(p, "--sample", help = "Sample name to label output files", type = "character", default = "topology_extracted.csv")
p <- add_argument(p, "--outdir", help = "Path to output files", type = "character", default = "E_epitopology")
# p <- add_argument(p, "--save_rdata_dir", help = "Path to save rData image", type = "character", default = "E_epitoplogy")

# Parse the command line arguments
argv <- parse_args(p)

## UPLOAD TOPOLOGY DATA
## OPTION 1
## IMPORT DATAFRAME from CCTOP PREDICTIONS
topology <- read.csv(file = argv$input_CCTOP, header = T)

### IMPORT DATAFRAME from EPIMERGER.PY
merged <- read.csv(file = argv$input_epitopes, header = T, sep = ";")

### FILTER EXTRACELLULAR DOMAINS
topology <- filter(topology, loc == "O" | loc == "S")

## CHARACTER TO NUMBER
topology$extracellular.start <- as.numeric(topology$extracellular.start)
topology$extracellular.end <- as.numeric(topology$extracellular.end)

### LOOP TO LABEL EXTRACELLULARITY
extracellular_list <- list()
extracellular_epitope <- c()
extracellular_label <- c()
for (x in 1:length(merged$Sequence)){
  for (z in 1:length(topology$extracellular.start)){
    y <- ifelse(merged$Start[x] >= topology$extracellular.start[z] & merged$End[x] <= topology$extracellular.end[z], "Extracellular", "Non-extracel")
    extracellular_epitope <- c(extracellular_epitope, y)
    extracellular_label <- c(extracellular_label, extracellular_epitope)
    extracellular_epitope <- c()
  }
  extracellular_list[[x]] <- extracellular_label
  extracellular_label <- c()
}

### ADD LIST TO DATAFRAME
merged$Extracellular <- extracellular_list

### TRANSFORM LIST TO CHARACTER
merged$Extracellular <- as.character(merged$Extracellular)

## EXPORT RESULTS
write.table(merged, sep = ";", file = paste0(argv$outdir, "/", argv$sample, sep = ""), quote = F, row.names = F)
print(paste("Find your output file at: ", argv$outdir, "/", argv$sample, sep = ""))

### EXPORT RDATA
# #p <- add_argument(p, "--save_rdata_dir", help = "Path to save rData image", type = "character", default = ".")
# dir.create(paste0(argv$save_rdata_dir, "rdata", sep = ""))
# save.image(paste0(argv$save_rdata_dir, "rdata/", argv$sample, ".rdata", sep = ""))
