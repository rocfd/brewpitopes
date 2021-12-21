## EPIGLYCAN EXTRACTOR
## GOAL: extract glycosilated positions from Net-N-Glyc and Net-O-glyc output data.

## GLYCAN PREDICTION
## N-GLYCOSILATIONS AT:
## https://services.healthtech.dtu.dk/service.php?NetNGlyc-1.0
## SAVE THE DATAFRAME HEADED SeqName	Position	Potential	Jury_agreement	NGlyc_result	Prediction
## AS CSV

## O-GLYCOSILATIONS AT:
## https://services.healthtech.dtu.dk/service.php?NetNGlyc-1.0
## SAVE THE DATAFRAME HEADED seqname	source	feature	start	end	score	strand	frame	comment
## AS CSV


library(dplyr)
library(tidyr)

## IMPORT O-GLYCOSILATION DATA
oglyc1 <- read.csv(header = F, col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "comment"), sep = "", "path/to/oglyc_results.csv")
colnames(oglyc1)
glimpse(oglyc1)
oglyc1$comment <- as.factor(oglyc1$comment)
table(oglyc1$comment)

## FILTER GLYCOSILATED POSITIONS
oglyc1_filt <- filter(oglyc1, comment == "#POSITIVE")
oglyc1_filt

## RENAME COLUMN
oglyc1_filt <- rename(oglyc1_filt, oglyc1_pos = "start")
glimpse(oglyc1_filt)

## EXTRACT GLYCOSILATED POSITIONS VECTOR
oglyc1_pos <- oglyc1_filt$oglyc1_pos
oglyc1_pos

## EXPORT AS TABLE
write.table(oglyc1_pos, "path/to/F_epiglycan/oglyc_positions.csv", row.names = T, quote = F, col.names = "Index   Glyc_positions")

'''
## IF THE PROTEIN EXCEEDS THE LENGTH HOLDED BY THE SERVER UPLOAD IN PARTS
## IMPORT DATA
oglyc2 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210914_R1AB_brewpitopes/F_epiglycan/20210914_R1AB_oglyc_part2.csv")
colnames(oglyc2)
glimpse(golyc)

## FILTER GLYCOSILATED POSITIONS
oglyc2_filt <- filter(oglyc2, comment == "#POSITIVE")

## RENAME COLUMN
oglyc2_filt <- rename(oglyc2_filt, oglyc2_pos = "start")

## EXTRACT GLYCOSILATED POSITIONS VECTOR
oglyc2_pos <- oglyc2_filt$oglyc2_pos
oglyc2_pos

## EXPORT AS TABLE
write.table(oglyc2_pos, "20210920_r1ab_oglyc_part2_positions.csv", row.names = T, quote = F, col.names = "Index   Glyc_positions")

## IMPORT DATA
oglyc3 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210914_R1AB_brewpitopes/F_epiglycan/20210914_R1AB_oglyc_part3.csv")
colnames(oglyc3)
glimpse(golyc)

## FILTER GLYCOSILATED POSITIONS
oglyc3_filt <- filter(oglyc3, comment == "#POSITIVE")

## RENAME COLUMN
oglyc3_filt <- rename(oglyc3_filt, oglyc3_pos = "start")

## EXTRACT GLYCOSILATED POSITIONS VECTOR
oglyc3_pos <- oglyc3_filt$oglyc3_pos
oglyc3_pos

## EXPORT AS TABLE
write.table(oglyc3_pos, "20210920_r1ab_oglyc_part3_positions.csv", row.names = T, quote = F, col.names = "Index   Glyc_positions")
'''

### IMPORT N-GLYCOSILATION DATA
nglyc <- read.csv("path/to/F_epiglycan/nglyc_results.csv", header = F, col.names = c("SeqName", "Position", "Potential", "Jury_agreement", "N-Glyc_result", "Prediction"))
colnames(nglyc)            
dim(nglyc)
View(nglyc)

## FILTER GLYCOSILATED POSITIONS
nglyc_filt <- filter(nglyc, Prediction == "++" | Prediction == "+" | Prediction == "+++")
nglyc_filt <- rename(nglyc_filt, nglyc_pos = "Position")

## EXTRACT GLYCOSILATED POSITIONS VECTOR
nglyc_pos <- nglyc_filt$nglyc_pos
nglyc_pos

# EXPORT DATA
write.table(nglyc_pos, "path/to/F_epiglycan/nglyc_positions.csv", row.names = T, quote = F, col.names = "Index   Glyc_positions")

'''
# IMPORT N GLYC DAT
nglyc2 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210914_R1AB_brewpitopes/F_epiglycan/20210914_R1AB_nglyc_part2.csv", header = F, col.names = c("SeqName", "Position", "Potential", "Jury_agreement", "N-Glyc_result", "Prediction"))
colnames(nglyc2)            
dim(nglyc2)
View(nglyc2)

## FILTER GLYCOSILATED POSITIONS
nglyc2_filt <- filter(nglyc2, Prediction == "++" | Prediction == "+" | Prediction == "+++")
nglyc2_filt <- rename(nglyc2_filt, nglyc2_pos = "Position")

## EXTRACT GLYCOSILATED POSITIONS VECTOR
nglyc2_pos <- nglyc2_filt$nglyc2_pos
nglyc2_pos

# EXPORT DATA
write.table(nglyc2_pos, "20210920_r1ab_nglyc_part2_positions.csv", row.names = T, quote = F, col.names = "Index   Glyc_positions")

# IMPORT N GLYC DAT
nglyc3 <- read.csv("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210914_R1AB_brewpitopes/F_epiglycan/20210914_R1AB_nglyc_part3.csv", header = F, col.names = c("SeqName", "Position", "Potential", "Jury_agreement", "N-Glyc_result", "Prediction"))
colnames(nglyc3)            
dim(nglyc3)
View(nglyc3)

## FILTER GLYCOSILATED POSITIONS
nglyc3_filt <- filter(nglyc3, Prediction == "++" | Prediction == "+" | Prediction == "+++")
nglyc3_filt <- rename(nglyc3_filt, nglyc3_pos = "Position")

## EXTRACT GLYCOSILATED POSITIONS VECTOR
nglyc3_pos <- nglyc3_filt$nglyc3_pos
nglyc3_pos

# EXPORT DATA
write.table(nglyc3_pos, "20210920_r1ab_nglyc_part3_positions.csv", row.names = T, quote = F, col.names = "Index   Glyc_positions")
'''

## Create new columns for glycosilations types
oglyc1_filt$Type <- "O-glyc"
#oglyc2_filt$Type <- "O-glyc"
#oglyc3_filt$Type <- "O-glyc"
nglyc_filt$Type <- "N-glyc"
#nglyc2_filt$Type <- "N-glyc"
#nglyc3_filt$Type <- "N-glyc"

## SUBSET DATAFRAMES
oglyc1_df <- oglyc1_filt[,c("oglyc1_pos", "Type")]
oglyc1_df <- rename(oglyc1_df, glyc_pos = "oglyc1_pos")

'''
oglyc2_df <- oglyc2_filt[,c("oglyc2_pos", "Type")]
oglyc2_df <- rename(oglyc2_df, glyc_pos = "oglyc2_pos")

oglyc3_df <- oglyc3_filt[,c("oglyc3_pos", "Type")]
oglyc3_df <- rename(oglyc3_df, glyc_pos = "oglyc3_pos")
'''

nglyc_df <- nglyc_filt[,c("nglyc_pos", "Type")]
nglyc_df <- rename(nglyc_df, glyc_pos = "nglyc_pos")

'''
nglyc2_df <- nglyc2_filt[,c("nglyc2_pos", "Type")]
nglyc2_df <- rename(nglyc2_df, glyc_pos = "nglyc2_pos")

nglyc3_df <- nglyc3_filt[,c("nglyc3_pos", "Type")]
nglyc3_df <- rename(nglyc3_df, glyc_pos = "nglyc3_pos")
'''

## MERGE DATAFRAMES
glyc_df <- rbind(nglyc_df, oglyc1_df)#nglyc2_df, nglyc3_df, oglyc1_df, oglyc2_df, oglyc3_df)

## EXPORT DATA
write.csv(nglyc_df, "path/to/F_epiglycan/XXX_glycosilated_positions.csv", row.names = F, quote = F)
