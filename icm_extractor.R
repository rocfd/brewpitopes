### ACCESSIBILITY ICM BROWSER

## GOAL: FILTER RESIDUES BY RSA >= 0.20

## LIBRARIES
library(dplyr, quietly = T, warn.conflicts = F)
library(tidyr, quietly = T, warn.conflicts = F)
library(stringr, quietly = T, warn.conflicts = F)
library(argparser, quietly = T, warn.conflicts = F)

### SCRIPT ARGUMENTS
# Create a parser
p <- arg_parser("ICM RSA PARSER")

# Add command line arguments
p <- add_argument(p, "--icm", help= "Path to icm output of ICM (.csv format)", type="character", default = "G_episurf/icm/rsa.csv")
p <- add_argument(p, "--sample", help = "Sample name to label output files (NO extension)", type = "character", default = "buried_positions")
p <- add_argument(p, "--outdir", help = "Path to output files", type = "character", default = "G_episurf")
p <- add_argument(p, "--chain", help = "Chain of interest of the target PDB", type = "character", default = "A")
#p <- add_argument(p, "--save_rdata_dir", help = "Path to save rData image", type = "character", default = "G_episurf")

# Parse the command line arguments
argv <- parse_args(p)

### LOAD DATA
icm <- read.csv(file = argv$icm)

## SPLIT RESIDUE COLUMN
#icm_res <- data.frame(do.call("rbind", strsplit(as.character(icm$residue), ".", fixed = TRUE)))

## BIND RESIDUE COLUMN
#icm <- cbind(icm, icm_res)

# ### RENAME AND REFORMAT COLUMNS
# icm <- rename(icm, chain = "X1")
# icm <- rename(icm, aa = "X2")
# icm <- rename(icm, position = "X3")
# icm$position <- as.numeric(as.character(icm$position))

### NEW PARSING

## MID FUNCTION
mid = function(text, start_num, num_char) {
  substr(text, start_num, start_num + num_char - 1)
}

### PARSE CHAIN
icm$chain <- mid(icm$residue, 1, 1)

### PARSE RESIDUE
icm$aa <- mid(icm$residue, 3, 1)

### PARSE POSITION
icm$position <- as.numeric(gsub("\\D", "", icm$residue))
#glimpse(icm)

### REORDER BY CHAIN
icm <- arrange(icm,chain)
icm <- filter(icm, chain == argv$chain)

### DEFINE TO FILTER NETSURF RESULTS
icm_position <- icm$position

### EXTRACT / FILTER BURIED POSITIONS
icm_buried <- filter(icm, rasa <= 0.2)

### EXPORT DATAFRAME
write.csv(icm_buried, row.names = F, quote = F, file = paste0(argv$outdir, "/", argv$sample, "_df.csv", sep = ""))

### EXPORT POSITIONS
icm_buried_positions <- icm_buried$position

## EXPORT DATA
write.table(icm_buried_positions, sep = ",", col.names = "buried", row.names = F, quote = F, file = paste0(argv$outdir, "/", argv$sample, "_list.csv", sep = ""))
print(paste("Find your output file with buried positions at: ", argv$outdir, "/", argv$sample, "_list.csv", sep = ""))

### EXPORT RDATA
#p <- add_argument(p, "--save_rdata_dir", help = "Path to save rData image", type = "character", default = ".")
#dir.create(paste0(argv$save_rdata_dir, "rdata", sep = ""))
#save.image(paste0(argv$save_rdata_dir, "rdata/", argv$sample, ".rdata", sep = ""))

