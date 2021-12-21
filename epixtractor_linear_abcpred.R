## EPIXTRACTOR LINEAR ABCPRED
# GOAL: extract epitopes from ABCpred results.
# SUBGOAL: merge result files from ABCpred results

## EPITOPE PREDICTION
# https://webs.iiitd.edu.in/raghava/abcpred/ABC_submission.html
# Predict epitopes using all the epitope length windows (10-20mers).
# Save the output results as csv (for instance: copying the output into excel and save as csv)
# Proceed with the following script.


library(dplyr)
library(purrr)

## IMPORT CSV - ABCPred RESULTS
abc_10 <- read.csv("path/to/10mers.csv")
abc_12 <- read.csv("path/to/12mers.csv")
abc_14 <- read.csv("path/to/14mers.csv")
abc_16 <- read.csv("path/to/16mers.csv")
abc_18 <- read.csv("path/to/18mers.csv")
abc_20 <- read.csv("path/to/20mers.csv")

## REMOVE COLUMN X
abc_10 <- select(abc_10, -X)
abc_12 <- select(abc_12, -X)
abc_14 <- select(abc_14, -X)
abc_16 <- select(abc_16, -X)
abc_18 <- select(abc_18, -X)
abc_20 <- select(abc_20, -X)

## MERGE RESULT FILES
abc_all <- rbind(abc_10, abc_12, abc_14, abc_16, abc_18, abc_20)

## FILTER by SCORE
abc_all <- filter(abc_all, Score >= 0.8)

## RENAME COLUMNS
colnames(abc_all)
#abc_all <- rename(abc_all, Length = "Window")
abc_all <- rename(abc_all, Start = "Start.position")
abc_all <- rename(abc_all, ABCscore = "Score")
#abc_all <- select(abc_all, -X)


## EXTRACT LENGTH
abc_all$Sequence <- as.character(abc_all$Sequence)
abc_all$Sequence
Length <- c()
for (x in 1:length(abc_all$Sequence)){
  z <- nchar(abc_all$Sequence[x])
  Length <- c(Length, z)
  print(Length)
}

abc_all <- cbind(abc_all, Length)

## EXTRACT END POSITION
End <- c()
for (x in 1:length(abc_all$Start)){
  z <- abc_all$Start[x] + abc_all$Length[x] -1
  End <- c(End, z)
  print(End)
}

abc_all <- cbind(abc_all, End)

## EXTRACT EPITOPE POSITIONS
abc_all$Positions = apply(abc_all,1, function(x) paste(x['Start']:x['End'],collapse=','))

## CREATE TOOL COLUMN
abc_all$Tool <- "ABCpred"

## EXPORT RESULTS AS CSV
write.table(abc_all, "path/to/abcpred/XXX_abcpred_allmers.csv", quote = F, row.names = F, sep = ";")
