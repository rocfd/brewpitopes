## EPIGLYCAN EXTRACTOR
## GOAL: extract glycosilated positions from Net-N-Glyc and Net-O-glyc output data.

## GLYCAN PREDICTION
## N-GLYCOSILATIONS AT:
## https://services.healthtech.dtu.dk/service.php?NetNGlyc-1.0
## SAVE THE DATAFRAME HEADED SeqName	Position	Potential	Jury_agreement	NGlyc_result	Prediction
## AS CSV

## O-GLYCOSILATIONS AT:
## https://services.healthtech.dtu.dk/service.php?NetOGlyc-4.0
## SAVE THE DATAFRAME HEADED seqname	source	feature	start	end	score	strand	frame	comment
## AS CSV

## LIBRARIES
library(dplyr, quietly = T, warn.conflicts = F)
library(tidyr, quietly = T, warn.conflicts = F)
library(argparser, quietly = T, warn.conflicts = F)

### SCRIPT ARGUMENTS
# Create a parser
p <- arg_parser("EPIGLYCAN POSITION EXTRACTOR")

# Add command line arguments
p <- add_argument(p, "--oglyc", help= "file to output from NetOglyc 4.0 (.txt format)", type="character", default = "brewpitopes/F_epiglycan/netoglyc")
p <- add_argument(p, "--nglyc", help= "file to output from NetNglyc 1.0 (.csv format)", type="character", default = "brewpitopes/F_epiglycan/netoglyc")
p <- add_argument(p, "--sample", help = "Sample name to label output files", type = "character", default = "glycan_positions.csv")
p <- add_argument(p, "--outdir", help = "file to output files", type = "character", default = "brewpitopes/F_epiglycan")
#p <- add_argument(p, "--save_rdata_dir", help = "file to save rData image", type = "character", default = "F_epiglycan")

# Parse the command line arguments
argv <- parse_args(p)

## IMPORT O-GLYCOSILATION DATA
oglyc1 <- read.table(header = F, col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "comment"), file = argv$oglyc, sep = "\t")
oglyc1$comment <- as.factor(oglyc1$comment)

## FILTER GLYCOSILATED POSITIONS
oglyc1_filt <- filter(oglyc1, comment == "#POSITIVE")

## RENAME COLUMN
oglyc1_filt <- rename(oglyc1_filt, oglyc1_pos = "start")

## EXTRACT GLYCOSILATED POSITIONS VECTOR
oglyc1_pos <- oglyc1_filt$oglyc1_pos

## EXPORT AS TABLE
write.table(oglyc1_pos, file = paste0(argv$outdir, "/", "oglyc_positions.csv", sep = ""), row.names = T, quote = F, col.names = "Index   Glyc_positions")

### IMPORT N-GLYCOSILATION DATA
nglyc <- read.table(file = argv$nglyc, header = F, col.names = c("SeqName", "Position", "Potential", "Jury_agreement", "N-Glyc_result", "Prediction"), sep = "\t")

## FILTER GLYCOSILATED POSITIONS
nglyc_filt <- filter(nglyc, Prediction == "++" | Prediction == "+" | Prediction == "+++")
nglyc_filt <- rename(nglyc_filt, nglyc_pos = "Position")

## EXTRACT GLYCOSILATED POSITIONS VECTOR
nglyc_pos <- nglyc_filt$nglyc_pos

# EXPORT DATA
write.table(nglyc_pos, file = paste0(argv$outdir, "/", "nglyc_positions.csv"), row.names = T, quote = F, col.names = "Index   Glyc_positions")

### LABEL GLYCOSYLATION TYPE
oglyc1_filt$Type <- "O-glyc"

nglyc_filt$Type <- "N-glyc"

## SUBSET DATAFRAMES
oglyc1_df <- oglyc1_filt[,c("oglyc1_pos", "Type")]
oglyc1_df <- rename(oglyc1_df, glyc_pos = "oglyc1_pos")

nglyc_df <- nglyc_filt[,c("nglyc_pos", "Type")]
nglyc_df <- rename(nglyc_df, glyc_pos = "nglyc_pos")

## MERGE DATAFRAMES
glyc_df <- rbind(nglyc_df, oglyc1_df)

## EXPORT DATA
write.csv(glyc_df, file = paste0(argv$outdir, "/", argv$sample, sep = ""), row.names = F, quote = F)
print(paste("Find your output files at: ", argv$outdir, sep = ""))

### EXPORT RDATA
#p <- add_argument(p, "--save_rdata_dir", help = "file to save rData image", type = "character", default = ".")
#dir.create(paste0(argv$save_rdata_dir, "rdata", sep = ""))
#save.image(paste0(argv$save_rdata_dir, "rdata/", argv$sample, ".rdata", sep = ""))
