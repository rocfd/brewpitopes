## CCTOP RESULTS PARSER
## GOAL: extract the extracellular domains from .xml CCTOP results.

# Created by: Roc Farriol
#    XX XX 2022
# 
# Modified by: Victor Montal
#    7 Feb 2024

## PREDICT TOPOLOGY USING CCTOP
## http://cctop.enzim.ttk.mta.hu/?_=/jobs/submit
## DOWNLOAD RESULTS AS .XML


# --
## Libraries
# --
rm(list = ls())
library(argparser, quietly = T, warn.conflicts = F)
library(dplyr, quietly = T, warn.conflicts = F)
library(XML, quietly = T, warn.conflicts = F)

# --
## SCRIPT PARSER
# --
# Create a parser
p <- arg_parser("XML CCTOP PARSER")
p <- add_argument(p, "--path", help= "Path to the brewpitope project", type="character")
argv <- parse_args(p)

# --
## Define defaults 
# --
argv$path <- "/home/vmontalb/Desktop/todel/brewpitopes"
ixml <- paste0(argv$path,"/E_epitopology/CCTOP/cctop.xml")
imerged <- paste0(argv$path,"/D_epimerger/merged.csv")
ofile <- paste0(argv$path,"/E_epitopology/CCTOP/cctop_domains.csv")

# --
## Operate CCTOP xml file
# --
data <- xmlParse(ixml)        # parse xml
xml_data <- xmlToList(data)   # transform to table
xml_df <- xml_data$CCTOPItem$Topology   # transform to df

# Loop over all topology domains and merge
if(is.null(xml_df)){
  print("STOPPER!! Your target protein has no predicted extracellular domains. Hence, neutralizing antibodies will not recognize it. You should consider another protein from your target organism.")
} else {
  for(i in 1:(length(xml_df)-1)) {
    assign(paste0("topdf.", i), as.data.frame(t(xml_df[[i]])))
  }
  # merge
  top_df <- mget(ls(pattern="^topdf\\.\\d+")) %>% bind_rows()
  # export
  write.csv(top_df, file=ofile , row.names=F, quote=F)
}
