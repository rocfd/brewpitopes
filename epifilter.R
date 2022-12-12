## EPIFILTER

## GOAL: IMPLEMENT THE FILTERS TO SELECT THE FINAL CANDIDATES
# TOPOLOGY = EXTRACELLULAR
# GLYCOSILATION = NON-GLYCOSILATED
# ACCESSIBILITY = ACCESSIBLE

library(dplyr, quietly = T, warn.conflicts = F)
library(argparser, quietly = T, warn.conflicts = F)

### SCRIPT ARGUMENTS
# Create a parser
p <- arg_parser("EPIFILTER")

# Add command line arguments
p <- add_argument(p, "--data", help= "Path to labeled epitope file output of episurf.py, by default 'access_extracted.csv' (.csv format)", type="character", default = "G_episurf/access_extracted.csv")
p <- add_argument(p, "--sample_df", help = "Sample name to label the output dataframe", type = "character", default = "brewpitopes_results_df")
p <- add_argument(p, "--sample_candidates", help = "Sample name to label the output list of candidates", type = "character", default = "brewpitopes_results_candidates")
p <- add_argument(p, "--outdir", help = "Path to output files", type = "character", default = "I_final_candidates")
#p <- add_argument(p, "--save_rdata_dir", help = "Path to save rData image", type = "character", default = "I_final_candidates")

# Parse the command line arguments
argv <- parse_args(p)

## IMPORT DATA
data <- read.csv(file = argv$data)

# REMOVE UNDESIRED COLUMNS
#data <- select(data, -X)
#data <- select(data, -Rank.1)

## FILTER BY SPECIFIED CONDITIONS
## TOPOLOGY
data_top <- filter(data, grepl('Extracellular', Extracellular))

## GLYCOSILATION
data_glyc <- filter(data_top, Glycosilation == "Non-glycosilated")

## ACCESSBILITY ICM
data_acc <- filter(data_glyc, accessibility_icm == "Accessible")

### FILTER BY LENGTH >= 5
data_len <- filter(data_acc, Length >= 5)

## FINAL CANDIDATES
data_final <- data_len
data_candidates <- data_len$Sequence

## REMOVE COLUMNS NAMED RANK
data_final <- select(data_final, -contains("Rank"))

## EXPORT DATA
write.table(data_final, sep = ";", row.names = F, quote = F, file = paste0(argv$outdir, "/", argv$sample_df, ".csv", sep = ""))
write.table(data_candidates, sep = ";", row.names = F, quote = F, file = paste0(argv$outdir, "/", argv$sample_candidates, ".csv", sep = ""))

### FINAL PRINT
print(paste("Find your output file at: ", argv$outdir, sep = ""))
