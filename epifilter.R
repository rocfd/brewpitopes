## EPIFILTER

# Created by: Roc Farriol
#    XX XX 2022
# 
# Modified by: Victor Montal
#    8 Feb 2024


## GOAL: IMPLEMENT THE FILTERS TO SELECT THE FINAL CANDIDATES
# TOPOLOGY = EXTRACELLULAR
# PTM = NON-GLYCOSILATED and NON-PHOSPHO
# ACCESSIBILITY = ACCESSIBLE
# LENGTH >= 5

# --
## Libraries
# --
rm(list = ls())
library(dplyr, quietly = T, warn.conflicts = F)
library(tidyr, quietly = T, warn.conflicts = F)
library(argparser, quietly = T, warn.conflicts = F)

# --
## SCRIPT PARSER
# --
p <- arg_parser("EPITOPE FILTER")
p <- add_argument(p, "--path", help= "Path to the brewpitope project", type="character")
argv <- parse_args(p)

# --
## Define defaults 
# --
iextract <- paste0(argv$path,"/G_episurf/access_extracted.csv")
ounfilt <- paste0(argv$path,"/H_epifilter/unfilter.csv")
ofilt <- paste0(argv$path,"/I_final_candidates/brewpitopes_results_df.csv")

# --
## Load data
# --
data <- read.csv(file=iextract)
write.table(data, sep = ";", row.names=F, quote=F, file = ounfilt)

# --
## Apply filters
# --
data_filt <- subset(data, grepl('Extracellular', data$Extracellular) &
                          PTM == "No-PTM" &
                          accessibility == "Accessible" &
                          Length > 5)

# --
## Dump results
# --
write.table(data_filt, sep=";", row.names=F, quote=F, file=ofilt)
print(paste("Find your filtered candidates at: ", ofilt, sep = ""))
