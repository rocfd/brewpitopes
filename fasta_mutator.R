## FASTA MUTATOR

## IMPLEMENT FOR mut VARIANT

## LIBRARIES
library(dplyr, quietly = T, warn.conflicts = F)
library(seqinr, quietly = T, warn.conflicts = F)
library(argparser, quietly = T, warn.conflicts = F)

### SCRIPT ARGUMENTS
# Create a parser
p <- arg_parser("FASTA MUTATOR")

# Add command line arguments
p <- add_argument(p, "--fasta", help= "Path to FASTA of the WT targeted protein (.fasta format)", type="character", default = ".")
p <- add_argument(p, "--mut", help= "Path to CSV containing the mutations (.csv format)", type="character", default = ".")
p <- add_argument(p, "--mut_header", help= "Header of the new (mutated) fasta. Do NOT include > !!", type="character", default = ".")
p <- add_argument(p, "--sample", help = "Sample name to label output files", type = "character", default = "fasta_mut.fasta")
p <- add_argument(p, "--outdir", help = "Path to output files", type = "character", default = "Z_fasta")
#p <- add_argument(p, "--save_rdata_dir", help = "Path to save rData image", type = "character", default = "Z_fasta")

# Parse the command line arguments
argv <- parse_args(p)

## IMPORT FASTA FROM REFERENCE SPIKE (UNIPROT)
fasta <- read.csv(header = F, col.names = "sequence", comment.char = ">", file = argv$fasta)

## MERGE FASTA SEQUENCE
y <- ""
for (x in fasta$sequence){
  y = paste(y,x, sep = "")
  print(y)
}
nchar(y)

## IMPORT MUTATIONS FROM THE CORRESPONDING VOC (IN THIS CASE mut)
mut <- read.csv(file = argv$mut) ## MUTATIONS FILE FOR mut : 20211203_spike_mut_vocs.csv

## EXTRACT MUT POSITION
mut$position <- regmatches(mut$mut_code, gregexpr("[[:digit:]]+", mut$mut_code))


### DEFINE LEFT AND RIGHT FUNCTIONS FOR EXTRACTION
left <- function(text, n) {
  substr(text, 1, n)
}

right <- function(text, n) {
  substr(text, nchar(text) - (n - 1), nchar(text))
}

### EXTRACT MUT AND WT
mut$mut_code <- as.character(mut$mut_code)
mut$mut <- right(mut$mut_code, 1)
mut$wt <- left(mut$mut_code, 1)

## MUTATE FASTA
### MUTATE FASTA LOOP
for (x in 1:length(mut$position)){
  if (x == 1) {
    z <- paste(substring(y,1,as.numeric(mut$position[x])-1), mut$mut[x], substring(y, as.numeric(mut$position[x])+1,), sep = "")  
  }
  else{
    z <- paste(substring(z,1,as.numeric(mut$position[x])-1), mut$mut[x], substring(z, as.numeric(mut$position[x])+1,), sep = "")  
  }
}

fasta_mut <- z

### CHECK MUTATED POSITIONS IN MUTATIONS DF
muts <- c()
for (x in 1:length(mut$mut)){
  m <- mut$mut[x]
  muts <- c(muts, m)
}

### CHECK MUTATED POSITIONS IN MUTATED FASTA
muts_f <- c()
for (x in 1:length(mut$position)){
  n <- substring(fasta, mut$position[x], mut$position[x])
  muts_f <- c(muts_f, n)
}


### CHECK IF EQUAL MUTATIONS AT MUTATED DF AND MUTATED FASTA
print("Mutation check:")
muts == muts_f

### DF TO FASTA FUNCTION

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

### FASTA HEADER
header <- argv$header

### FASTA DATAFRAME
fasta_df <- data.frame(name = argv$mut_header, seq = fasta_mut)

### DF to FASTA
writeFasta(fasta_df, paste0(argv$outdir, "/", argv$sample, sep = ""))
