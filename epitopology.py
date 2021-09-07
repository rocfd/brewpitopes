## EPITOPOLOGY

## ## Locate epitopes on protein surface based on subcellular location information.

library(dplyr)
library(tidyr)
library(stringr)

## DOWNLOAD UNIPROT INFORMATION FOR SPIKE
## LOCATION
"/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210112_brewpitopes/C_epitopology/uniprot-reviewed_yes+AND+spike+sars2.tab"

## IMPORT DATAFRAME
surf <- read.table("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210112_brewpitopes/C_epitopology/uniprot-reviewed_yes+AND+spike+sars2_short.txt", sep = "\t", header = T)
glimpse(surf)


## PARSE
# deleting punctuations
surf$Topological.domain<-gsub("[[:punct:][:blank:]]+", " ", surf$Topological.domain)
# deleting trailing space
surf$Topological.domain<-gsub("\\n"," ", surf$Topological.domain)
surf$Topological.domain

for (i in 1:length(surf$Topological.domain)){    
  #comparing split with the state abbreviation   
  dom1 <- sub("note Extracellular.*", "", surf$Topological.domain)      
  dom2 <- sub(".*TOPO DOM ", "", dom1)
  dom3 <- gsub(" ", "-", dom2, fixed = TRUE)
  dom4 <- gsub('.{1}$', '', dom3)
  dom5 <- strsplit(dom4, "-")
  surf$extracellular.start <- dom5[[1]][1]
  surf$extracellular.end <- dom5[[1]][2]
  surf$extracellular.positions = apply(surf,1, function(x) paste(x['extracellular.start']:x['extracellular.end'],collapse=','))
  #surf$extracellular.end <- dom5[[2]]
  #adding states to the new column   
  #surf$location[i] <- state_split[1]  
  print(dom5)
}

#glimpse(surf)
#sub("note Extracellular.*", "", surf$Topological.domain)

## GLIMPSE MERGED
glimpse(merged)
#merged$extracellular = apply(merged,1, function(x) ifelse(merged$Start[x] >= surf$extracellular.start[x] & merged$End[x] <= surf$extracellular.end[x], "Extracellular", "Non-Extracellular"))

## CHARACTER TO NUMBER
surf$extracellular.start <- as.numeric(surf$extracellular.start)
surf$extracellular.end <- as.numeric(surf$extracellular.end)

## LABEL EXTRACELLULAR
extracellular <- c()
for (x in 1:length(merged$Sequence)){
  y <- ifelse(merged$Start[x] >= surf$extracellular.start[1] & merged$End[x] <= surf$extracellular.end[1], "Extracellular", "Non-extracellular")
  extracellular <- c(extracellular, y)
  print(extracellular)
}

## BIND LABEL TO DATAFRAME
merged <- cbind(merged, extracellular)
glimpse(merged)
merged <- rm(merged, "extracellular")

## EXPORT DATAFRAME
write.table(merged, sep = ";", "/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20210112_brewpitopes/C_epitopology/epitopology_spike.csv", quote = F, row.names = F)



