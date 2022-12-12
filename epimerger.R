### EPIMERGER

# GOAL: MERGE LINEAR EPITOPES FROM bepiPRED AND ABCPRED AND DISCOTOPE

## LIBRARIES
library(tidyr, quietly = T, warn.conflicts = F)
library(tibble, quietly = T, warn.conflicts = F)
library(dplyr, quietly = T, warn.conflicts = F)
library(data.table, quietly = T, warn.conflicts = F)
library(argparser, quietly = T, warn.conflicts = F)

### SCRIPT ARGUMENTS
# Create a parser
p <- arg_parser("EPIMERGER")

# Add command line arguments
p <- add_argument(p, "--abcpred", help= "Path to abcpred extracted epitopes (.csv format)", type="character", default = "./brewpitopes/C_epixtractor/abcpred_extracted.csv")
p <- add_argument(p, "--bepipred", help= "Path to bepipred extracted epitopes (.csv format)", type="character", default = "./brewpitopes/C_epixtractor/bepipred_extracted.csv")
p <- add_argument(p, "--discotope", help= "Path to discotope extracted epitopes (.csv format)", type="character", default = "./brewpitopes/C_epixtractor/discotope_extracted.csv")
p <- add_argument(p, "--sample", help = "Sample name to label output files", type = "character", default = "merged")
p <- add_argument(p, "--outdir", help = "Path to output files", type = "character", default = "./brewpitopes/D_epimerger")
# p <- add_argument(p, "--save_rdata_dir", help = "Path to save rData image", type = "character", default = "./brewpitopes/D_epimerger")

# Parse the command line arguments
argv <- parse_args(p)

## IMPORT RESULTS
## IMPORT ABCPRED RESULTS EXTRACTED USING EPIXTRACTOR
abc <- read.csv(file = argv$abcpred,sep = ";")
## IMPORT BEPIPRED RESULTS EXTRACTED USING EPIXTRACTOR
bepi <- read.csv(file = argv$bepipred, sep = ";")

## ADD EXTRA COLUMNS TO bepi
bepi$Sequence <- as.character(bepi$Sequence)
bepi$Tool <- "Bepipred 2.0"

## CALCULATE SEQUENCE LENGTH
len <- c()
for (x in 1:length(bepi$Sequence)){
  x <- nchar(bepi$Sequence[x])
  len <- c(len, x)
}

## ADD LENGTH TO RESULTS
bepi$Length <- len
#glimpse(bepi)

## RENAME SCORE COLUMNS
abc <- rename(abc, Score = "ABCscore")
bepi <- rename(bepi, Score = "BebiScore")

## BIND DATAFRAMES INTO LINEAR RESULTS
linear <- rbind(bepi, abc)

## IMPORT STRUCTURAL RESULTS
disc <- read.csv(file = argv$discotope, sep = ";")

## ADD VARIABLES
disc$Sequence <- as.character(disc$Sequence)

## CALCULATE EPITOPE LENGTH
len2 <- c()
for (x in 1:length(disc$Sequence)){
  x <- nchar(disc$Sequence[x])
  len2 <- c(len2, x)
}

disc$Length <- len2

## ADD PREDICTION TOOL
disc$Tool <- "Discotope 2.0"

## BIND DATAFRAMES
merged <- rbind(linear, disc) 

## FILTER LENGTH 0 FROM ABCPRED (optional)
merged <- filter(merged, Length != 0)

## EXPORT DATAFRAME
write.table(merged, file = paste0(argv$outdir, "/", argv$sample, ".csv"), row.names = F, quote = F, sep = ";")
print(paste("Find your output file at:", argv$outdir, sep = " ")

### EXPORT RDATA
#p <- add_argument(p, "--save_rdata_dir", help = "Path to save rData image", type = "character", default = ".")
# dir.create(paste0(argv$save_rdata_dir, "/rdata"))
# save.image(paste0(argv$save_rdata_dir, "/rdata/", argv$sample, ".rdata"))
