## EPIGLYCAN EXTRACTOR
## GOAL: extract Post-translational modifications positions from:
##    Net-N-Glyc, Net-O-glyc, Net-C-glyc, NetPhos

# Created by: Roc Farriol
#    XX XX 2022
# 
# Modified by: Victor Montal
#    8 Feb 2024

## N-GLYCOSILATIONS AT:
## https://services.healthtech.dtu.dk/service.php?NetNGlyc-1.0
## SAVE THE DATAFRAME HEADED SeqName	Position	Potential	Jury_agreement	NGlyc_result	Prediction
## AS CSV

## O-GLYCOSILATIONS AT:
## https://services.healthtech.dtu.dk/service.php?NetOGlyc-4.0
## SAVE THE DATAFRAME HEADED seqname	source	feature	start	end	score	strand	frame	comment
## AS CSV

## C-GLYCOSILATIONS AT:
## https://services.healthtech.dtu.dk/services/NetCGlyc-1.0/
## SAVE THE TABLE WITHOUT HEADER (i.e #) AS CSV

## PHOSPHORYLATION AT:
## https://services.healthtech.dtu.dk/services/NetPhos-3.1/
## SAVE THE TABLE WITHOUT HEADER (i.e #) AS CSV


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
p <- arg_parser("EPI POST-TRANSLATIONAL MODIFICATIONS POSITION EXTRACTOR")
p <- add_argument(p, "--path", help= "Path to the brewpitope project", type="character")
argv <- parse_args(p)


# --
## Define defaults 
# --
ioglyc <- paste0(argv$path,"/F_epiptm/netoglyc/oglyc_predict.csv")
inglyc <- paste0(argv$path,"/F_epiptm/netnglyc/nglyc_predict.csv") 
icglyc <- paste0(argv$path,"/F_epiptm/netcglyc/cglyc_predict.csv")
iphos <- paste0(argv$path,"/F_epiptm/netphos/phospho_predict.csv")
ofile <- paste0(argv$path,"/F_epiptm/ptm_positions.csv")

# --
## Load dataFrames, filter
# --
oglyc <- read.csv(header=F, col.names =c("seqname","source","feature","start","end","score","strand","frame","comment"), file=ioglyc, sep="\t")
nglyc <- read.csv(header=F, file=inglyc, col.names=c("SeqName", "Position", "Potential", "Jury_agreement", "N-Glyc_result", "Prediction"), sep="")
cglyc <- read.csv(header=F, file=icglyc, col.names=c("seqname","source","feature","start","end","score","ratio","Prediction"), sep="")
phosp <- read.csv(header=F, file=iphos, col.names=c("seqname","source","Phospho_type","start","end","score","ratio","oratio","Prediction"), sep="")


# Glycosilation
oglyc$Type <- "O-Glyc"
oglyc <- filter(oglyc, comment == "#POSITIVE")
oglyc <- rename(oglyc, ptm_pos = "start")

nglyc$Type <- "N-Glyc"
nglyc <- filter(nglyc, Prediction == "++" | Prediction == "+" | Prediction == "+++")
nglyc <- rename(nglyc, ptm_pos = "Position")

cglyc$Type <- "C-Glyc"
cglyc <- filter(cglyc, Prediction == "W")
cglyc <- rename(cglyc, ptm_pos = "start")

# Phosphorylation
phosp <- filter(phosp, Prediction == "YES")

phosp_update <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(phosp_update) <- c("ptm_pos","Type")

idx = 1
for (cresid in unique(phosp$start)){
  # subset residue
  tmpsubset <- subset(phosp, start == cresid)
  # concat all phospho-types
  phosp_type = c()
  for (cphospho in unique(tmpsubset$Phospho_type) )
  {
    phosp_type = c(phosp_type, cphospho)
  }
  phosp_type = paste(phosp_type, collapse=",")
  phosp_update[idx,"ptm_pos"] <- cresid
  phosp_update[idx,"Type"] <- phosp_type
  idx = idx+1
}

# --
## Merge tables
# --
oglyc <- oglyc[,c("ptm_pos","Type")]
nglyc <- nglyc[,c("ptm_pos","Type")]
cglyc <- cglyc[,c("ptm_pos","Type")]

merged_ptm <- bind_rows(oglyc,nglyc,cglyc,phosp_update)

# --
## Dump tables
# --  
write.csv2(merged_ptm, file=ofile, row.names=F, quote=F)
print(paste("Find your output files at: ", ofile, sep = ""))
