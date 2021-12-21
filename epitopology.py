## EPITOPOLOGY
## GOAL: label epitopes on the viral surface based on topology information.

library(dplyr)
library(tidyr)
library(stringr)

## UPLOAD TOPOLOGY DATA
## OPTION 1
## IMPORT DATAFRAME from CCTOP PREDICTIONS
#surf <- read.csv("path/to/cctop_prediction.csv", header = T)
glimpse(surf)

## MANUAL DOMAIN UPLOAD (info extracted from Uniprot)
## ADD ONLY EXTRAVIRAL DOMAINS
surf <- data.frame(extracellular.start = c(2247,2797,3121, 3608, 3656, 3751), extracellular.end = c(2317, 3044, 3127, 3608, 3673, 3778))
surf

## CHARACTER TO NUMBER
surf$extracellular.start <- as.numeric(surf$extracellular.start)
surf$extracellular.end <- as.numeric(surf$extracellular.end)

## LABEL EXTRACELLULAR DOMAINS (each extracellular correspond to one extracellular domain. Add or remove according to the extraviral domains of the protein.)
extracellular1 <- c()
for (x in 1:length(merged$Sequence)){
  y <- ifelse(merged$Start[x] >= surf$extracellular.start[1] & merged$End[x] <= surf$extracellular.end[1], "Extracellular", "Non-extracellular")
  extracellular1 <- c(extracellular1, y)
  print(extracellular1)
}

extracellular2 <- c()
for (x in 1:length(merged$Sequence)){
  y <- ifelse(merged$Start[x] >= surf$extracellular.start[2] & merged$End[x] <= surf$extracellular.end[2], "Extracellular", "Non-extracellular")
  extracellular2 <- c(extracellular2, y)
  print(extracellular2)
}

extracellular3 <- c()
for (x in 1:length(merged$Sequence)){
  y <- ifelse(merged$Start[x] >= surf$extracellular.start[2] & merged$End[x] <= surf$extracellular.end[3], "Extracellular", "Non-extracellular")
  extracellular3 <- c(extracellular3, y)
  print(extracellular3)
}

extracellular4 <- c()
for (x in 1:length(merged$Sequence)){
  y <- ifelse(merged$Start[x] >= surf$extracellular.start[4] & merged$End[x] <= surf$extracellular.end[4], "Extracellular", "Non-extracellular")
  extracellular4 <- c(extracellular4, y)
  print(extracellular4)
}

extracellular5 <- c()
for (x in 1:length(merged$Sequence)){
  y <- ifelse(merged$Start[x] >= surf$extracellular.start[5] & merged$End[x] <= surf$extracellular.end[5], "Extracellular", "Non-extracellular")
  extracellular5 <- c(extracellular5, y)
  print(extracellular5)
}

extracellular6 <- c()
for (x in 1:length(merged$Sequence)){
  y <- ifelse(merged$Start[x] >= surf$extracellular.start[6] & merged$End[x] <= surf$extracellular.end[6], "Extracellular", "Non-extracellular")
  extracellular6 <- c(extracellular6, y)
  print(extracellular6)
}

### UPLOAD MERGED INFORMATION (output from epimerger containing Bepipred, ABCpred and Discotope epitopes).
merged <- read.csv("path/to/epimerger_output.csv")

### LABEL EPITOPES AS EXTRACELLULAR
merged <- cbind(merged, extracellular1, extracellular2, extracellular3, extracellular4, extracellular5, extracellular6)
glimpse(merged)


## EXPORT RESULTS
write.table(merged, sep = ";", "/E_epitopology/XXX_epitopology.csv", quote = F, row.names = F)

'''
### UNIPROT PARSER (TO REFINE)
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
'''
