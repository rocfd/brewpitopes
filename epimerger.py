### EPIMERGER

# GOAL: MERGE LINEAR EPITOPES FROM BEBIPRED AND ABCPRED AND DISCOTOPE

## LIBRARIES
library(tidyr)
library(tibble)
library(dplyr)
library(data.table)

## IMPORT RESULTS
## IMPORT ABCPRED RESULTS EXTRACTED USING EPIXTRACTOR
abc <- read.csv("path/to/abcpred/results_allmers.csv",sep = ";")
## IMPORT BEPIPRED RESULTS EXTRACTED USING EPIXTRACTOR
bebi <- read.csv(file = "path/to/bepipred/results_Out.csv", sep = ";")
colnames(abc)
colnames(bebi)

## ADD EXTRA COLUMNS TO BEBI
bebi$Sequence <- as.character(bebi$Sequence)
bebi$Tool <- "Bepipred 2.0"

## CALCULATE SEQUENCE LENGTH
len <- c()
for (x in 1:length(bebi$Sequence)){
  x <- nchar(bebi$Sequence[x])
  len <- c(len, x)
}

## ADD LENGTH TO RESULTS
bebi$Length <- len
glimpse(bebi)

## RENAME SCORE COLUMNS
abc <- rename(abc, Score = "ABCscore")
bebi <- rename(bebi, Score = "BebiScore")

## BIND DATAFRAMES INTO LINEAR RESULTS
linear <- rbind(bebi, abc)
glimpse(linear)

## IMPORT STRUCTURAL RESULTS
disc <- read.csv("path/to/C_epixtractor/discotope_Out.csv", sep = ";")

### ADD MORE RESULTS IF MORE STRUCTURES WERE ANALYZED
'''
disc2 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210920_R1A_brewpitopes/C_epixtractor/7cz4_Out.csv", sep = ";")
disc3 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210920_R1A_brewpitopes/C_epixtractor/7d47_Out.csv", sep = ";")
disc4 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210920_R1A_brewpitopes/C_epixtractor/7exm_Out.csv", sep = ";")
#disc5 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210914_R1AB_brewpitopes/C_epixtractor/7kag_819-929_Out.csv", sep = ";")
#disc6 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210914_R1AB_brewpitopes/C_epixtractor/6zlw_1-180_Out.csv", sep = ";")
#disc7 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210914_R1AB_brewpitopes/C_epixtractor/6zct_4236-4384_Out.csv", sep = ";")
#disc8 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210914_R1AB_brewpitopes/C_epixtractor/6wuu_1562-1879_Out.csv", sep = ";")
#disc9 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210914_R1AB_brewpitopes/C_epixtractor/6w4b_4141-4253_Out.csv", sep = ";")
#disc10 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210914_R1AB_brewpitopes/C_epixtractor/5s71_6543-6798_Out.csv", sep = ";")
'''

### MERGE THE DIFFERENT STRUCTURAL RESULTS
disc_all <- rbind(disc)#, disc2, disc4, disc4)
colnames(disc_all)
glimpse(disc_all)

## ADD VARIABLES
disc_all$Sequence <- as.character(disc_all$Sequence)
#disc_all$Length = apply(disc,1, function(x) nchar(disc_all$Sequence[x]))

## CALCULATE EPITOPE LENGTH
len2 <- c()
for (x in 1:length(disc_all$Sequence)){
  x <- nchar(disc_all$Sequence[x])
  len2 <- c(len2, x)
}

disc_all$Length <- len2

## ADD PREDICTION TOOL
disc_all$Tool <- "Discotope 2.0"
glimpse(disc_all)
glimpse(linear)

## BIND DATAFRAMES
merged <- rbind(linear, disc_all) #since no discotope predictions
glimpse(merged)

## FILTER LENGTH 0 FROM ABCPRED (optional)
merged <- filter(merged, Length != 0)

## EXPORT DATAFRAME
write.table(merged, "D_epimerger/XXX_linear_structural.csv", row.names = F, quote = F, sep = ";")
