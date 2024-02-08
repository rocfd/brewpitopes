## EPITOPOLOGY
## GOAL: label epitopes on the viral topologyace based on topology information.

# Created by: Roc Farriol
#    XX XX 2022
# 
# Modified by: Victor Montal
#    7 Feb 2024

# --
## Libraries
# --
rm(list = ls())
library(dplyr, quietly = T, warn.conflicts = F)
library(tidyr, quietly = T, warn.conflicts = F)
library(stringr, quietly = T, warn.conflicts = F)
library(argparser, quietly = T, warn.conflicts = F)


# --
## SCRIPT PARSER
# --
p <- arg_parser("EPITOPOLOGY from CCTOP")
p <- add_argument(p, "--path", help= "Path to the brewpitope project", type="character")
argv <- parse_args(p)

# --
## Define defaults 
# --
icctop <- paste0(argv$path,"/E_epitopology/CCTOP/cctop_domains.csv")
imerged <- paste0(argv$path,"/D_epimerger/merged.csv")
ofile <- paste0(argv$path,"/E_epitopology/topology_extracted.csv")
  
# --
## Load topology and merged df
# --
topology <- read.csv(file=icctop, header = T)
topology <- filter(topology, loc == "O" | loc == "S")      # filter segment 
topology <- rename(topology, extracellular.start = "from")
topology <- rename(topology, extracellular.end = "to")
topology$extracellular.start <- as.numeric(topology$extracellular.start)
topology$extracellular.end <- as.numeric(topology$extracellular.end)

merged <- read.csv(file=imerged, header = T, sep = ";")

# --
## Label Epitopes based on extracellularity
# --
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

merged$Extracellular <- extracellular_list
merged$Extracellular <- as.character(merged$Extracellular)


# --
## Export results
# --
write.table(merged, sep=";", file=ofile, quote = F, row.names = F)
print(paste("Find your output file at: ", ofile))
