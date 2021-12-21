## EPIFILTER

## GOAL: IMPLEMENT THE FILTERS TO SELECT THE FINAL CANDIDATES
# TOPOLOGY = EXTRACELLULAR
# GLYCOSILATION = NON-GLYCOSILATED
# ACCESSIBILITY = ACCESSIBLE

library(dplyr)

## IMPORT DATA
data <- read.csv("path/to/XXX_epitopology_glycans_surf.csv", sep = ";")
glimpse(data)

# REMOVE UNDESIRED COLUMNS
#data <- select(data, -X)
data <- select(data, -Rank.1)
glimpse(data)

## FILTER BY SPECIFIED CONDITIONS
## TOPOLOGY
data_top <- filter(data, extracellular == "Extracellular")

## GLYCOSILATION
data_glyc <- filter(data_top, Glycosilation == "Non-glycosilated")

## ACCESSBILITY ICM
data_acc <- filter(data_glyc, accessibility == "Accessible")

### FILTER BY LENGTH >= 5
data_len <- filter(data_acc, Length >= 5)
dim(data_len)

## FINAL CANDIDATES
data_final <- data_len
data_candidates <- data_len$Sequence

## EXPORT DATA
write.table(data_final, sep = ";", row.names = F, quote = F, file = "path/to/I_final_candidates/XXX_final_dataframe.csv")
write.table(data_candidates, sep = ";", row.names = F, quote = F, file = "path/to/I_final_candidates/XXX_final_candidates.csv")
