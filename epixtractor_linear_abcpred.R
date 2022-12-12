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
p <- add_argument(p, "--input_10mers", help= "Path to ABCPRED 10-mers output dataset (.csv format)", type="character", default = ".")
p <- add_argument(p, "--input_12mers", help= "Path to ABCPRED 12-mers output dataset (.csv format)", type="character", default = ".")
p <- add_argument(p, "--input_14mers", help= "Path to ABCPRED 14-mers output dataset (.csv format)", type="character", default = ".")
p <- add_argument(p, "--input_16mers", help= "Path to ABCPRED 16-mers output dataset (.csv format)", type="character", default = ".")
p <- add_argument(p, "--input_18mers", help= "Path to ABCPRED 18-mers output dataset (.csv format)", type="character", default = ".")
p <- add_argument(p, "--input_20mers", help= "Path to ABCPRED 20-mers output dataset (.csv format)", type="character", default = ".")
p <- add_argument(p, "--sample", help = "Sample name to label output files", type = "character", default = "abcpred_results_extracted")
p <- add_argument(p, "--outdir", help = "Path to output files", type = "character", default = "brewpitopes/C_epixtractor")
#p <- add_argument(p, "--rdata", help = "Path to save rData image", type = "character", default = "brewpitopes/A_linear_predictions/abcpred")

# Parse the command line arguments
argv <- parse_args(p)

## IMPORT CSV - ABCPred RESULTS
abc_10 <- read.csv(argv$input_10mers)
abc_12 <- read.csv(argv$input_12mers)
abc_14 <- read.csv(argv$input_14mers)
abc_16 <- read.csv(argv$input_16mers)
abc_18 <- read.csv(argv$input_18mers)
abc_20 <- read.csv(argv$input_20mers)

## REMOVE COLUMN X
#abc_10 <- select(abc_10, -X)
#abc_12 <- select(abc_12, -X)
#abc_14 <- select(abc_14, -X)
#abc_16 <- select(abc_16, -X)
#abc_18 <- select(abc_18, -X)
#abc_20 <- select(abc_20, -X)

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
#abc_all$Sequence <- as.character(abc_all$Sequence)
#abc_all$Sequence
#Length <- c()
#for (x in 1:length(abc_all$Sequence)){
#  z <- nchar(abc_all$Sequence[x])
#  Length <- c(Length, z)
  #print(Length)
#}

#abc_all <- cbind(abc_all, Length)

## EXTRACT END POSITION
End <- c()
for (x in 1:length(abc_all$Start)){
  z <- abc_all$Start[x] + abc_all$Length[x] -1
  End <- c(End, z)
  #print(End)
}

abc_all <- cbind(abc_all, End)

## EXTRACT EPITOPE POSITIONS
abc_all$Positions = apply(abc_all,1, function(x) paste(x['Start']:x['End'],collapse=','))

## CREATE TOOL COLUMN
abc_all$Tool <- "ABCpred"

## EXPORT RESULTS AS CSV
write.table(abc_all, file = paste0(argv$outdir, "/", argv$sample, ".csv"), quote = F, row.names = F, sep = ";")

### EXPORT RDATA
# #p <- add_argument(p, "--save_rdata_dir", help = "Path to save rData image", type = "character", default = ".")
# dir.create(paste0(argv$save_rdata_dir, "/rdata", sep = ""))
# save.image(paste0(argv$save_rdata_dir, "/rdata/", argv$sample, ".rdata", sep = ""))
