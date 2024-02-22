## EPIPTM
## GOAL: LABEL THE EPITOPES DEPENDING ON THEIR POST-TRANSLATIONAL MODIFICATIONS STATUS

# Created by: Roc Farriol
#    XX XX 2022
# 
# Modified by: Victor Montal
#    8 Feb 2024

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
p <- arg_parser("EPITOPE POST-TRANSLATIONAL MODIFICATIONS LABELER")
p <- add_argument(p, "--path", help= "Path to the brewpitope project", type="character")
argv <- parse_args(p)


# --
## Define defaults 
# --
iptm <- paste0(argv$path,"/F_epiptm/ptm_positions.csv")
itopol <- paste0(argv$path,"/E_epitopology/topology_extracted.csv")
ofile <- paste0(argv$path,"/F_epiptm/ptm_extracted.csv")

# --
## Load dataFrames
# --
ptm <- read.csv( file=iptm, sep=";")
topol <- read.csv( file=itopol, sep=";")

outdf <- topol
outdf$PTM <- ""

# --
## Label Epitopes
# --
# Loop each epitope
for (idx_epi in 1:nrow(outdf)){
  cpos_epi = outdf[idx_epi,"Positions"]
  cpos_epi <- as.integer(str_split(cpos_epi,",")[[1]])
  
  # Loop each PTM
  cptm <- c()
  for (idx_ptm in 1:1:nrow(ptm)){
    cpos_ptm <- ptm[idx_ptm,"ptm_pos"]
    ctype_ptm <- ptm[idx_ptm,"Type"]
    
    if (cpos_ptm %in% cpos_epi){
      cptm <- c(cptm,ctype_ptm)
    }
  }
  if (is_empty(cptm))
  {
    cptm <- "No-PTM"
  }else{
  cptm <- paste(cptm,collapse=",")
  }
  outdf[idx_epi,"PTM"] <- cptm
}


# --
## Dump tables
# --  
write.csv2(outdf, file=ofile, row.names=F, quote=F)
print(paste("Find your output files at: ", ofile, sep = ""))
