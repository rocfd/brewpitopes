## EPIXTRACTOR LINEAR ABCPRED
# GOAL: extract epitopes from ABCpred results.
# SUBGOAL: merge result files fro ABCpred results

library(dplyr)
library(purrr)

## IMPORT CSV - ABCPred RESULTS
abc_10 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210902_brewpitopes_setup/test_data/abcpred/P0DTC9_w10_abcpred_result.csv")
abc_12 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210902_brewpitopes_setup/test_data/abcpred/P0DTC9_w12_abcpred_result.csv")
abc_14 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210902_brewpitopes_setup/test_data/abcpred/P0DTC9_w14_abcpred_result.csv")
abc_16 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210902_brewpitopes_setup/test_data/abcpred/P0DTC9_w16_abcpred_result.csv")
abc_18 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210902_brewpitopes_setup/test_data/abcpred/P0DTC9_w18_abcpred_result.csv")
abc_20 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210902_brewpitopes_setup/test_data/abcpred/P0DTC9_w20_abcpred_result.csv")

## MERGE RESULT FILES
abc_all <- rbind(abc_10, abc_12, abc_14, abc_16, abc_18, abc_20)

## FILTER by SCORE
abc_all <- filter(abc_all, Score >= 0.5)

## Extract END position
colnames(abc_all)
abc_all <- rename(abc_all, Length = "Window")
abc_all <- rename(abc_all, Start = "Start.position")
abc_all <- rename(abc_all, ABCscore = "Score")

End <- c()
for (x in 1:length(abc_all$Start)){
  z <- abc_all$Start[x] + abc_all$Length[x] -1
  End <- c(End, z)
  print(End)
}

abc_all <- cbind(abc_all, End)

## Extract epitope POSITIONS
abc_all$Positions = apply(abc_all,1, function(x) paste(x['Start']:x['End'],collapse=','))
#df=   data.frame ( a=rep(1, 4), b = rep(10, 4))
#df$c = apply(df,1, function(x) paste(x['a']:x['b'],collapse=','))

## Create tool column
abc_all$Tool <- "ABCpred"

## Export csv
write.table(abc_all, "/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210112_brewpitopes/A_epixtractor/epixtractor_linear_abcpred/eprotein_abcpred_allmers_extract.csv", quote = F, row.names = F, sep = ";")



