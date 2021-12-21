## CCTOP RESULTS PARSER
## GOAL: extract the extracellular domains from .xml CCTOP results.

## PREDICT TOPOLOGY USING CCTOP
## http://cctop.enzim.ttk.mta.hu/?_=/jobs/submit
## DOWNLOAD RESULTS AS .XML

## XML PARSER
require(XML)

## LOAD XML
data <- xmlParse("path/to/E_epitoplogy/results_cctop.xml")

## TRANSFORM TO LIST
xml_data <- xmlToList(data)
View(xml_data)

## LIST TO DATAFRAME
xml_df <- xml_data$Topology
xml_df

### EXTRACT ALL THE TOPOLOGY DOMAINS
### ADD ONE DATAFRAME PER DOMAIN
## EXTRACT REGION TO DATAFRAME
xml_df1 <- as.data.frame(xml_df[1])
xml_df1 <- t(xml_df1)
rownames(xml_df1) <- NULL
dim(xml_df1)

## EXTRACT REGION TO DATAFRAME
xml_df2 <- as.data.frame(xml_df[2])
xml_df2 <- t(xml_df2)
rownames(xml_df2) <- NULL
dim(xml_df2)

## EXTRACT REGION TO DATAFRAME
xml_df3 <- as.data.frame(xml_df[3])
xml_df3 <- t(xml_df3)
rownames(xml_df3) <- NULL
dim(xml_df3)

## EXTRACT REGION TO DATAFRAME
xml_df4 <- as.data.frame(xml_df[4])
xml_df4 <- t(xml_df4)
rownames(xml_df4) <- NULL
View(xml_df4)
dim(xml_df4)

## EXTRACT REGION TO DATAFRAME
'''xml_df5 <- as.data.frame(xml_df[5])
xml_df5 <- t(xml_df5)
rownames(xml_df5) <- NULL
View(xml_df5)
dim(xml_df5)'''

## MERGE INTO DATAFRAME
xml_df_all <- rbind(xml_df1, xml_df2, xml_df3, xml_df4)
View(xml_df_all)
class(xml_df_all)

## EXPORT DATAFRAME
write.csv(xml_df_all, "path/to/E_epitopology/CCTOP/cctop_results_extracted.csv")
