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

### IF EMPTY
if(dim(bepi)[1] == 0){
  bepi[1,] <- NA
  bepi[,"Tool"] <- "Bepipred 2.0"
  bepi <- bepi[0,]
} else {
  bepi$Tool <- "Bepipred 2.0"
}

## CALCULATE SEQUENCE LENGTH
len <- c()

### IF EMPTY
if(dim(bepi)[1] == 0){
  bepi[1,] <- NA
  bepi[,"Length"] <- NA
  bepi <- bepi[0,]
} else {
  for (x in 1:length(bepi$Sequence)){
  x <- nchar(bepi$Sequence[x])
  len <- c(len, x)
  }
bepi$Length <- len
}

## ADD LENGTH TO RESULTS
#bepi$Length <- len
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

### IF EMPTY
if(dim(disc)[1] == 0){
  disc[1,] <- NA
  disc[,"Tool"] <- "Discotope 2.0"
  disc <- disc[0,]
} else {
  disc$Tool <- "Discotope 2.0"
}

## CALCULATE SEQUENCE LENGTH
len2 <- c()

### IF EMPTY
if(dim(disc)[1] == 0){
  disc[1,] <- NA
  disc[,"Length"] <- NA
  disc <- disc[0,]
} else {
  for (x in 1:length(disc$Sequence)){
  x <- nchar(disc$Sequence[x])
  len2 <- c(len2, x)
  }
disc$Length <- len2
}

## BIND DATAFRAMES
merged <- rbind(linear, disc) 

## FILTER LENGTH 0 FROM ABCPRED (optional)
merged <- filter(merged, Length != 0)

## EXPORT DATAFRAME
write.table(merged, file = paste0(argv$outdir, "/", argv$sample, ".csv"), row.names = F, quote = F, sep = ";")

## FINAL PRINT
print(paste("Find your output merged file at: ", argv$outdir, "/", argv$sample, ".csv", sep = ""))

## CHECK IF EPITOPES WERE PREDICTED BY ANY TOOL
### CONTAINS EPITOPES PRINT
if(dim(merged)[1] == 0){
  print("None of the three predictors used (ABCpred, Bepipred 2.0 and Discotope 2.0) was able to predict a single epitope in your target protein. You should STOP THE PIPELINE and look for other epitope predictors.")
} else {
  print("One or more tools predicted epitopes in your target protein. Go ahead!")
}

### EXPORT RDATA
#p <- add_argument(p, "--save_rdata_dir", help = "Path to save rData image", type = "character", default = ".")
# dir.create(paste0(argv$save_rdata_dir, "/rdata"))
# save.image(paste0(argv$save_rdata_dir, "/rdata/", argv$sample, ".rdata"))
