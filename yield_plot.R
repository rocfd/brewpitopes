### PLOT BREWPITOPE OVERVIEW
### GOAL: observe the yield of the brewpitopes pipeline.

### libraries
library(dplyr, quietly = T, warn.conflicts = F)
library(ggplot2, quietly = T, warn.conflicts = F)
library(ggthemr, quietly = T, warn.conflicts = F)
library(argparser, quietly = T, warn.conflicts = F)

### SCRIPT ARGUMENTS
# Create a parser
p <- arg_parser("YIELD PLOT")

# Add command line arguments
#p <- add_argument(p, "--bepi", help= "Path to Bepipred2.0 file output (.csv format)", type="character", default = "/brewpitopes/A_linear_predictions/bepipred/bepipred_results.csv")
#p <- add_argument(p, "--abc", help= "Path to ABCpred file output (.csv format)", type="character", default = "/brewpitopes/A_linear_predictions/abcpred/abcpred_results.csv")
#p <- add_argument(p, "--disc", help= "Path to Discotope2.0 file output (.csv format)", type="character", default = "/brewpitopes/B_structural_predictions/discotope/discotope_results.csv")
#p <- add_argument(p, "--merged", help= "Path to Epimerger.R file output (.csv format)", type="character", default = "/brewpitopes/D_epimerger/merged.csv")
#p <- add_argument(p, "--top", help= "Path to Epitopology.R file output (.csv format)", type="character", default = "/brewpitopes/E_epitopology/topology_extracted.csv")
#p <- add_argument(p, "--glyc", help= "Path to epiglycan.py file output (.csv format)", type="character", default = "/brewpitopes/F_epiglycan/glycan_extracted.csv")
#p <- add_argument(p, "--acc", help= "Path to episurf.py file output (.csv format)", type="character", default = "/brewpitopes/G_episurft/access_extracted.csv")
p <- add_argument(p, "--data", help = "Path to dataframe of unfiltered candidates (.csv format)", type = "character", default = "H_epifilter/brewpitopes_unfiltered_df.csv")
p <- add_argument(p, "--merged", help = "Path to dataframe output of epimerger.R (.csv format)", type = "character", default = "D_epimerger/merged.csv")
p <- add_argument(p, "--eregs", help = "Path to dataframe output of epimerger.R (.csv format)", type = "character", default = "K_epitope_regions/epitope_regions_extracted.csv")
p <- add_argument(p, "--outdir", help = "Path to output files", type = "character", default = "J_plots")
#p <- add_argument(p, "--save_rdata_dir", help = "Path to save rData image", type = "character", default = "I_final_candidates")

# Parse the command line arguments
argv <- parse_args(p)

### REVISE DATA
#bepi <- read.csv(file = argv$bepi)
#abc <- read.csv(file = argv$abc)
#disc <- read.csv(file = argv$disc)
#merged <- read.csv(file = argv$merged)
#topology <- read.csv(file = argv$top)
#glycan <- read.csv(file = argv$glyc)
#surface <- read.csv(file = argv$acc)

## UPLOAD FINAL CANDIDATES
data <- read.csv(file = argv$data, sep = ";")

## EXTRACT Nº CANDIDATES
## BEPIPRED 2.0
bepi <- filter(data, Tool == "Bepipred 2.0")
bepi_c <- length(bepi$Sequence)

## ABCpred
abc <- filter(data, Tool == "ABCpred")
abc_c <- length(abc$Sequence)

## DISCOTOPE
disc <- filter(data, Tool == "Discotope 2.0")
disc_c <- length(disc$Sequence)

## MERGED LINEAR + STRUCTURAL
merged <- read.csv(file = argv$merged, sep = ";")
merged_c <- length(merged$Sequence)

### TOPOLOGY
top <- filter(data, grepl('Extracellular', Extracellular))
top_c <- length(top$Sequence)

## GLYCAN
glyc <- filter(top, Glycosilation == "Non-glycosilated")
glyc_c <- length(glyc$Sequence)

## SURFACE
surf <- filter(glyc, accessibility_icm == "Accessible")
surf_c <- length(surf$Sequence)

## LENGTH
len <- filter(surf, Length >= 5)
len_c <- length(len$Sequence)

## LOAD epregs DATA
epregs <- read.csv(file = argv$eregs)

## epregs
epregs_c <- length(epregs$Sequence)

### PLOT
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## DATAFRAME
#bar_epregs <- data.frame(Dataset = c("Discotope 2.0", "Bepipred 2.0", "ABCpred", "Merged", "Topology", "Glycosilation", "Accessibility", "Length", "Epitope Regions"), Yield =  c(disc_c, bepi_c, abc_c, merged_c, top_c, glyc_c, surf_c, len_c, epregs_c))

## LOCK FACTOR LEVELS
#bar_epregs$Dataset <- factor(bar_epregs$Dataset, levels = c("Discotope 2.0", "Bepipred 2.0", "ABCpred", "Merged", "Topology", "Glycosilation", "Accessibility", "Length", "Epitope Regions"))

# PLOT
#ggthemr("dust")
#bar_plot_epregs <- ggplot(bar_epregs, aes(y=Yield, x = Dataset)) +
#  geom_bar(position="dodge", stat='identity') +#fill = "#0072B2", color = "#000000", alpha = 0.5) +
#  geom_text(aes(label=Yield),position=position_dodge(width = 1), vjust = -1) +
  #scale_fill_manual(values=c("#628395", "#C5D86D", "#DB2763"))
#  scale_colour_manual(values = cbPalette) +
#  ggtitle("Brewpitopes Yield") +
#  xlab("Pipeline steps") + ylab("Epitopes") +
#  scale_fill_hue(l=40)
#bar_plot_epregs

## EXPORT
#ggsave(filename= paste0(argv$outdir, "/", "brewpitopes_yield.pdf", sep = ""), plot=bar_plot_epregs, width=12, height=8, units="in", bg = "white")

### STACKED BARPLOT
## DATAFRAME UNFILTERED DATA
ep_pred <- data.frame(table(data$Tool))
ep_pred <- rename(ep_pred, Tool = "Var1")
ep_pred <- rename(ep_pred, Yield = "Freq")
ep_pred$Dataset <- "Epitope Prediction"

### DATAFRAME TOPOLOGY
topology <- data.frame(table(top$Tool))
topology <- rename(topology, Tool = "Var1")
topology <- rename(topology, Yield = "Freq")
topology$Dataset <- "Topology"

### DATAFRAME glycosylation
glycosylation <- data.frame(table(glyc$Tool))
glycosylation <- rename(glycosylation, Tool = "Var1")
glycosylation <- rename(glycosylation, Yield = "Freq")
glycosylation$Dataset <- "Glycosylation"

### DATAFRAME accessibility
accessibility <- data.frame(table(surf$Tool))
accessibility <- rename(accessibility, Tool = "Var1")
accessibility <- rename(accessibility, Yield = "Freq")
accessibility$Dataset <- "Accessibility"

### DATAFRAME long
long <- data.frame(table(len$Tool))
long <- rename(long, Tool = "Var1")
long <- rename(long, Yield = "Freq")
long$Dataset <- "Length"

### DATAFRAME EPITOPE REGIONS
#epitope_regions$Yield[1] <- if(epregs$Score_Bepipred_2_0 == 0 & epregs$Score_ABCpred != 0){
#                          length(epregs$Sequence)
#                          }

#glimpse(epregs)

epitope_regions <- data.frame(Yield = c(1,2,3), Tool = c("ABCpred", "Bepipred 2.0", "Discotope 2.0"), Dataset = "Epitope Regions")

## ABC COUNT IN EPITOPE REGIONS
abc_count <- 0
for(i in 1:length(epregs$Sequence)){
                            if(epregs$Score_Bebipred_2_0[i] == 0 & epregs$Score_ABCpred[i] != 0){
                            abc_count <- abc_count + 1
                            }
}

epitope_regions$Yield[1] <- abc_count

## BEPI COUNT IN EPITOPE REGIONS
bepi_count <- 0
for(i in 1:length(epregs$Sequence)){
                            if(epregs$Score_Bebipred_2_0[i] != 0){
                            bepi_count <- bepi_count + 1
                            }
}

epitope_regions$Yield[2] <- bepi_count

## DISC COUNT IN EPITOPE REGIONS
disc_count <- 0
for(i in 1:length(epregs$Sequence)){
                            if(epregs$Score_Bebipred_2_0[i] == 0 & epregs$Score_Discotope_2_0[i] != -10){
                            disc_count <- disc_count + 1
                            }
}

epitope_regions$Yield[3] <- disc_count

#epitope_regions$Yield[2] <- if(epregs$Sequence while epregs$Score_Bepipred_2_0 != 0){
#                          length(epregs$Sequence)
#                          }
#epitope_regions$Yield[3] <- if(epregs$Sequence while epregs$Score_Bepipred_2_0 == 0 & epregs$Score_Discotope_2_0 != -10){
#                          length(epregs$Sequence)
#                          }
epitope_regions$Tool[1] <- "ABCpred"
epitope_regions$Tool[2] <- "Bepipred 2.0"
epitope_regions$Tool[3] <- "Discotope 2.0"
epitope_regions$Dataset <- "Epitope Regions"

### RBIND DATAFRAME
bar_stacked <- rbind(ep_pred, topology, glycosylation, accessibility, long, epitope_regions)
#glimpse(bar_stacked)

## LOCK FACTOR LEVELS
bar_stacked$Dataset <- factor(bar_stacked$Dataset, levels = c("Epitope Prediction", "Topology", "Glycosylation", "Accessibility", "Length", "Epitope Regions"))
#levels(bar_stacked$Dataset)

### STACKED PLOT
ggthemr("dust")
bar_stacked_plot <- ggplot(bar_stacked, aes(x = Dataset, y = Yield, fill = Tool)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label=Yield),position = "stack", vjust = -0.3, hjust = 0.5) + # position_stack(vjust = .5)
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Pipeline steps") + ylab("Epitopes")

 ## EXPORT
ggsave(filename= paste0(argv$outdir, "/", "brewpitopes_yield.pdf", sep = ""), plot=bar_stacked_plot, width=12, height=8, units="in", bg = "white")

### FINAL PRINT
print(paste("Find your output plot at: ", argv$outdir, "/", "brewpitopes_yield.pdf", sep = ""))
