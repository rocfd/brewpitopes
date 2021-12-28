## FASTA MUTATOR

## IMPLEMENT FOR gamma VARIANT

## LIBRARIES
library(dplyr)
library(seqinr)

## IMPORT FASTA FROM REFERENCE SPIKE (UNIPROT)
fasta <- read.csv(header = F, col.names = "sequence", comment.char = ">", "path/to/spike_fasta.txt")
glimpse(fasta)

## MERGE FASTA SEQUENCE
y <- ""
for (x in fasta$sequence){
  y = paste(y,x, sep = "")
  print(y)
}
nchar(y)

## IMPORT MUTATIONS FROM THE CORRESPONDING VOC (IN THIS CASE GAMMA)
gamma <- read.csv("path/to/gamma_mutations.csv") ## MUTATIONS FILE FOR GAMMA : 20211203_spike_gamma_vocs.csv
glimpse(gamma)

## EXTRACT MUT POSITION
gamma$position <- regmatches(gamma$mut_code, gregexpr("[[:digit:]]+", gamma$mut_code))
gamma

### EXTRACT MUT AND WT
### DEFINE LEFT AND RIGHT FUNCTIONS FOR EXTRACTION
left <- function(text, n) {
  substr(text, 1, n)
}

right <- function(text, n) {
  substr(text, nchar(text) - (n - 1), nchar(text))
}

gamma$mut_code <- as.character(gamma$mut_code)
gamma$mut <- right(gamma$mut_code, 1)
gamma$wt <- left(gamma$mut_code, 1)
gamma

## REMOVE INSERTION (FOR THE MOMENT)
##gamma <- gamma[1:36,]

## MUTATE FASTA
### MUTATE FASTA LOOP
for (x in 1:length(gamma$position)){
  if (x == 1) {
    z <- paste(substring(y,1,as.numeric(gamma$position[x])-1), gamma$mut[x], substring(y, as.numeric(gamma$position[x])+1,), sep = "")  
  }
  else{
    z <- paste(substring(z,1,as.numeric(gamma$position[x])-1), gamma$mut[x], substring(z, as.numeric(gamma$position[x])+1,), sep = "")  
  }
}

fasta <- z

### CHECK MUTATED POSITIONS IN MUTATIONS DF
muts <- c()
for (x in 1:length(gamma$mut)){
  m <- gamma$mut[x]
  muts <- c(muts, m)
}
muts

### CHECK MUTATED POSITIONS IN MUTATED FASTA
muts_f <- c()
for (x in 1:length(gamma$position)){
  n <- substring(fasta, gamma$position[x], gamma$position[x])
  muts_f <- c(muts_f, n)
}

### CHECK IF EQUAL MUTATIONS AT MUTATED DF AND MUTATED FASTA
muts == muts_f
muts_f # YES


## SAVE AS TEXT
sink(file = "Z_fasta/gamma_fasta_seq.txt")
fasta
sink()

## ONCE SAVED REMOVE ""
