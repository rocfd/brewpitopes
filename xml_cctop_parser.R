## CCTOP RESULTS PARSER
## GOAL: extract the extracellular domains from .xml CCTOP results.

## PREDICT TOPOLOGY USING CCTOP
## http://cctop.enzim.ttk.mta.hu/?_=/jobs/submit
## DOWNLOAD RESULTS AS .XML

## LIBRARIES
library(argparser, quietly = T, warn.conflicts = F)
library(dplyr, quietly = T, warn.conflicts = F)

## XML PARSER
library(XML, quietly = T, warn.conflicts = F)

### SCRIPT ARGUMENTS
# Create a parser
p <- arg_parser("XML CCTOP PARSER")

# Add command line arguments
p <- add_argument(p, "--xml", help= "Path to XML output of CCTOP (.xml format)", type="character", default = ".")
p <- add_argument(p, "--sample", help = "Sample name to label output files", type = "character", default = "cctop_domains.csv")
p <- add_argument(p, "--outdir", help = "Path to output files", type = "character", default = "E_epitopology/CCTOP")
#p <- add_argument(p, "--save_rdata_dir", help = "Path to save rData image", type = "character", default = "E_epitopology/CCTOP")

# Parse the command line arguments
argv <- parse_args(p)

## LOAD XML
data <- xmlParse("/Users/rocfarriolduran/Desktop/BSC/B_SARS-CoV-2/A_SARS2_DEVELOPMENT/20221212_dockerized_final_test/brewpitopes/Projects/marina/brewpitopes/E_epitopology/CCTOP/cctop.xml")

## TRANSFORM TO LIST
xml_data <- xmlToList(data)

## LIST TO DATAFRAME
xml_df <- xml_data$Topology
glimpse(xml_df)

### EXTRACT ALL THE TOPOLOGY DOMAINS
### ADD ONE DATAFRAME PER DOMAIN
## EXTRACT ALL REGIONS INTO SEPPARATE DATAFRAMES
## SKIP LAST DATAFRAME SINCE IT DOES NOT CONTAIN TOPOLOGY INFORMATION

### IF EMPTY
if(is.null(xml_df)){
  print("STOPPER!! Your target protein has no predicted extracellular domains. Hence, neutralizing antibodies will not recognize it. You should consider another protein from your target organism.")
} else {
  for(i in 1:(length(xml_df)-1)) {
    assign(paste0("topdf.", i), as.data.frame(t(xml_df[[i]])))
  }
  
  ### MERGE DOMAINS INTO A SINGLE DATAFRAME
  top_df <- mget(ls(pattern="^topdf\\.\\d+")) %>% bind_rows()
  
  ## EXPORT DATAFRAME
  write.csv(top_df, file = paste0(argv$outdir, "/", argv$sample, sep = ""), row.names = F, quote = F)

}

### EXPORT RDATA
#p <- add_argument(p, "--save_rdata_dir", help = "Path to save rData image", type = "character", default = ".")
#dir.create(paste0(argv$save_rdata_dir, "/rdata", sep = ""))
#save.image(paste0(argv$save_rdata_dir, "/rdata/", argv$sample, ".rdata", sep = ""))
