## EPIXTRACTOR LINEAR ABCPRED
# GOAL: extract epitopes from ABCpred results.
# SUBGOAL: merge result files from ABCpred results

## EPITOPE PREDICTION
# https://webs.iiitd.edu.in/raghava/abcpred/ABC_submission.html
# Predict epitopes using all the epitope length windows (10-20mers).
# Save the output results as csv (for instance: copying the output into excel and save as csv)
# Proceed with the following script.

library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(purrr, quietly = TRUE, warn.conflicts = FALSE)
library(argparser, quietly = TRUE, warn.conflicts = FALSE)

### SCRIPT ARGUMENTS
# Create a parser
p <- arg_parser("EPIXTRACTOR LINEAR ABCPRED")

# Add command line arguments
p <- add_argument(p, "--input_allmers", help= "Path to ABCPRED all-mers output dataset (.csv format)", type="character", default = "A_linear_predictions/abcpred/abcpred.tsv")
p <- add_argument(p, "--sample", help = "Sample name to label output files", type = "character", default = "abcpred_results_extracted")
p <- add_argument(p, "--outdir", help = "Path to output files", type = "character", default = "brewpitopes/C_epixtractor")
#p <- add_argument(p, "--rdata", help = "Path to save rData image", type = "character", default = "brewpitopes/A_linear_predictions/abcpred")

# Parse the command line arguments
argv <- parse_args(p)

## IMPORT CSV - ABCPred RESULTS
abc_all <- read.table(argv$input_allmers, sep = "\t", header = TRUE, row.names = NULL)  #, fill=TRUE, row.names = NULL)
#abc_10 <- read.table(argv$input_10mers, sep = "\t")
#abc_12 <- read.table(argv$input_12mers, sep = "\t")
#abc_14 <- read.table(argv$input_14mers, sep = "\t")
#abc_16 <- read.table(argv$input_16mers, sep = "\t")
#abc_18 <- read.table(argv$input_18mers, sep = "\t")
#abc_20 <- read.table(argv$input_20mers, sep = "\t")

### RENAME COLUMNS ACCORDINGLY
colnames(abc_all) <- colnames(abc_all)[2:ncol(abc_all)]  

## REMOVE NA COLUMN
abc_all <- abc_all[ , - ncol(abc_all)]

## REMOVE COLUMN X
#abc_10 <- select(abc_10, -X)
#abc_12 <- select(abc_12, -X)
#abc_14 <- select(abc_14, -X)
#abc_16 <- select(abc_16, -X)
#abc_18 <- select(abc_18, -X)
#abc_20 <- select(abc_20, -X)

## MERGE RESULT FILES
#abc_all <- rbind(abc_10, abc_12, abc_14, abc_16, abc_18, abc_20)

## FILTER by SCORE
abc_all <- filter(abc_all, Score >= 0.8)

## RENAME COLUMNS
#colnames(abc_all)
#abc_all <- rename(abc_all, Length = "Window")
abc_all <- rename(abc_all, Start = "Start.position")
abc_all <- rename(abc_all, ABCscore = "Score")
#abc_all <- select(abc_all, -X)

### FILTER BLANK SEQUENCES
#abc_all <- filter(abc_all, !is.na(Sequence))

## EXTRACT LENGTH
abc_all$Sequence <- as.character(abc_all$Sequence)
len <- c()
### IF EMPTY
if(dim(abc_all)[1] == 0){
  abc_all[1,] <- NA
  abc_all[,"Length"] <- NA
  abc_all <- abc_all[0,]
} else {
  for (x in 1:length(abc_all$Sequence)){
  x <- nchar(abc_all$Sequence[x])
  len <- c(len, x)
  }
abc_all$Length <- len
}

## EXTRACT END POSITION
End <- c()
if(dim(abc_all)[1] == 0){
  abc_all[1,] <- NA
  abc_all[,"End"] <- NA
  abc_all <- abc_all[0,]
} else {for (x in 1:length(abc_all$Start)){
  z <- abc_all$Start[x] + abc_all$Length[x] -1
  End <- c(End, z)
  }
abc_all$End <- End
}

## EXTRACT EPITOPE POSITIONS
Positions <- c()
if(dim(abc_all)[1] == 0){
  abc_all[1,] <- NA
  abc_all[,"Positions"] <- NA
  abc_all <- abc_all[0,]
} else {
abc_all$Positions = apply(abc_all,1, function(x) paste(x['Start']:x['End'],collapse=','))
}

## CREATE TOOL COLUMN
if(dim(abc_all)[1] == 0){
  abc_all[1,] <- NA
  abc_all[,"Tool"] <- NA
  abc_all <- abc_all[0,]
} else {
  abc_all$Tool <- "ABCpred"
}

## EXPORT RESULTS AS CSV
write.table(abc_all, file = paste0(argv$outdir, "/", argv$sample, ".csv"), quote = F, row.names = F, sep = ";")

## FINAL PRINT
print(paste("Find your output file at: ", argv$outdir, "/", argv$sample, ".csv", sep = ""))

### CONTAINS EPITOPES PRINT
if(dim(abc_all)[1] == 0){
  print("ABCpred did not find any epitope in your target sequence. You will get an empty dataframe but you can continue the pipeline with the other predictors")
} else {
  print("ABCpred found one or more epitopes in your target sequence. Go ahead!")
}

### EXPORT RDATA
# #p <- add_argument(p, "--save_rdata_dir", help = "Path to save rData image", type = "character", default = ".")
# dir.create(paste0(argv$save_rdata_dir, "/rdata", sep = ""))
# save.image(paste0(argv$save_rdata_dir, "/rdata/", argv$sample, ".rdata", sep = ""))
