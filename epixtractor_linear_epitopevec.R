## EPIXTRACTOR LINEAR EpitopeVec
# GOAL: extract epitopes from EpitopeVec results.
#       and merge them into a single table

## EPITOPE PREDICTION FROM EPITOPEVEC
# https://github.com/hzi-bifo/epitope-prediction
# Predict epitopes using a range from 10 to 20 kmers
# Save each of the .txt outputs in: A_linear_predictions/epitopevec/epitopevec_pred_10.txt
#                                   A_linear_predictions/epitopevec/epitopevec_pred_17.txt
# Proceed with the following script.

#library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(argparser, quietly = TRUE, warn.conflicts = FALSE)
library(purrr, quietly = TRUE, warn.conflicts = FALSE)

### SCRIPT ARGUMENTS
# Create a parser
p <- arg_parser("EPIXTRACTOR LINEAR EpitopeVec")

# Add command line arguments
#p <- add_argument(p, "--input_allmers", help= "Path to ABCPRED all-mers output dataset (.csv format)", type="character", default = "A_linear_predictions/abcpred/abcpred.tsv")
p <- add_argument(p, "--path", help= "Path to the brewpitope project", type="character")

# Parse the command line arguments
argv <- parse_args(p)

## IMPORT multiple CSV 
ipath <- paste0(argv$path,"/A_linear_predictions/epitopevec/")
filelist <- list.files(ipath, full.names = TRUE)
datalist <- lapply(filelist, FUN=read.table, header=TRUE)
datafr = do.call("rbind", datalist) 

## RENAME COLUMNS
datafr <- rename(datafr, Sequence = "peptide")
datafr <- rename(datafr, EpitopeVecscore = "score")
datafr <- rename(datafr, Start = "start")
datafr <- rename(datafr, End = "end")

## MAKE DUMMY POSITIONAL VECTOR
seqs = mapply(FUN = function(a, b) {
  paste(a:b, collapse=",")
}, a = datafr$Start, b = datafr$End)
seqs =  paste0('[',seqs,']')
datafr$Positions <- seqs

datafr$Length <- datafr$End - datafr$Start

## EXPORT RESULTS AS CSV
opath <- paste0(argv$path,"/C_epixtractor/epitopevec_results_extracted.csv")
write.table(datafr, opath, quote = F, row.names = F, sep = ";")

## FINAL PRINT
print(paste("Find your output file at: ", opath))
