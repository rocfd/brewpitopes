### EPIMERGER

# GOAL: MERGE LINEAR EPITOPES FROM BEBIPRED AND ABCPRED AND DISCOTOPE

## LIBRARIES
library(tidyr)
library(tibble)
library(dplyr)

## IMPORT RESULTS
abc <- read.csv( "/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210112_brewpitopes/A_epixtractor/epixtractor_linear_abcpred/eprotein_abcpred_allmers_extract.csv",sep = ";")
bebi <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210902_brewpitopes_setup/test_data/eprotein_bebipred_results_Out.csv", sep = ";")
colnames(abc)
colnames(bebi)

## ADD EXTRA COLUMNS TO BEBI
bebi$Sequence <- as.character(bebi$Sequence)
bebi$Length = apply(bebi,1, function(x) nchar(bebi$Sequence))
bebi$Tool <- "Bebipred 2.0"

## RENAME SCORE COLUMNS
abc <- rename(abc, Score = "ABCscore")
bebi <- rename(bebi, Score = "BebiScore")

## BIND DATAFRAMES
linear <- rbind(bebi, abc)
glimpse(linear)

## IMPORT STRUCTURAL RESULTS
disc <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210902_brewpitopes_setup/test_data/6VXX_Discotope_chainABC_Out.csv", sep = ";")
colnames(disc)

## ADD VARIABLES
disc$Tool <- "Discotope 2.0"
disc$Sequence <- as.character(disc$Sequence)
disc$Length = apply(disc,1, function(x) nchar(disc$Sequence))

## BIND DATAFRAMES
merged <- rbind(linear, disc)
glimpse(merged)

## EXPORT DATAFRAME
write.table(merged, "/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210112_brewpitopes/B_epimerger/epimerger_linear_struct.csv", row.names = F, quote = F, sep = ";")
