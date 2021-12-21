### PLOT BREWPITOPE OVERVIEW
### GOAL: observe the yield of the brewpitopes pipeline.

### libraries
library(dplyr)
library(ggplot2)
library(ggthemr)

### REVISE DATA
bebi
abc <- abc_all
linear
#structural
#lin_struct <- merged
topology <- data_top
glycan <- data_glyc
surface <- data_acc

## EXTRACT Nº CANDIDATES
## BEBI
glimpse(bebi)
bebi_c <- length(bebi$Sequence)
bebi_c # 22

## ABCpred
glimpse(abc)
#abc_filt <- filter(abcpred, ABCscore >= 0.9)
abc_c <- length(abc$Sequence)
abc_c #194

## DISCOTOPE
glimpse(disc_all)
disc_c <- length(disc_all$Sequence)
disc_c #14

## LINEAR STRUCTURAL
ls_c <- length(merged$Sequence)
ls_c #230

### TOPOLOGY
top_c <- length(data_top$Sequence)
top_c #222

## GLYCAN
glyc_c <- length(data_glyc$Sequence)
glyc_c # 164

## SURFACE
surf_c <- length(data_acc$Sequence)
surf_c # 13

## LENGTH
len_c <- length(data_len$Sequence)
len_c #11


## DATAFRAME
bar <- data.frame(Dataset = c("Discotope", "Bebipred", "ABCpred", "Linear-Structural", "Topology", "Glycosilation", "Accessibility", "Length"), Yield =  c(disc_c, bebi_c, abc_c, ls_c, top_c, glyc_c, surf_c, len_c))
glimpse(bar)
View(bar)

## LOCK FACTOR LEVELS
bar$Dataset <- factor(bar$Dataset, levels = c("Discotope", "Bebipred", "ABCpred", "Linear-Structural", "Topology", "Glycosilation", "Accessibility", "Length"))
levels(bar$Dataset)

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# PLOT
library(ggplot2)
library(ggthemr)
ggthemr("dust")
bar_plot <- ggplot(bar, aes(y=Yield, x = Dataset)) +
  geom_bar(position="dodge", stat='identity') +#fill = "#0072B2", color = "#000000", alpha = 0.5) +
  geom_text(aes(label=Yield),position=position_dodge(width = 1), vjust = -1) +
  #scale_fill_manual(values=c("#628395", "#C5D86D", "#DB2763")) +
  scale_colour_manual(values = cbPalette) +
  ggtitle("Spike Omicron Pipeline Barplot") +
  xlab("Pipeline steps") + ylab("Epitopes") +
  scale_fill_hue(l=40)
bar_plot


## EXPORT
ggsave(filename="J_plots/XXX_barplot.pdf", plot=bar_plot, width=12, height=8, units="in", bg = "white")

## LOAD CONTIG DATA
contig <- read.csv("path/to/K_contigs/XXX_contigs.csv")
glimpse(contig)

## CONTIG
contig_c <- length(contig$Sequence)
contig_c # 6


## DATAFRAME
bar_contig <- data.frame(Dataset = c("Discotope", "Bebipred", "ABCpred", "Linear-Structural", "Topology", "Glycosilation", "Accessibility", "Length", "Contig"), Yield =  c(disc_c, bebi_c, abc_c, ls_c, top_c, glyc_c, surf_c, len_c, contig_c))
glimpse(bar_contig)
View(bar_contig)

## LOCK FACTOR LEVELS
bar_contig$Dataset <- factor(bar_contig$Dataset, levels = c("Discotope", "Bebipred", "ABCpred", "Linear-Structural", "Topology", "Glycosilation", "Accessibility", "Length", "Contig"))
levels(bar_contig$Dataset)

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# PLOT
library(ggplot2)
library(ggthemr)
ggthemr("dust")
bar_plot_contig <- ggplot(bar_contig, aes(y=Yield, x = Dataset)) +
  geom_bar(position="dodge", stat='identity') +#fill = "#0072B2", color = "#000000", alpha = 0.5) +
  geom_text(aes(label=Yield),position=position_dodge(width = 1), vjust = -1) +
  #scale_fill_manual(values=c("#628395", "#C5D86D", "#DB2763")) +
  scale_colour_manual(values = cbPalette) +
  ggtitle("Spike Omicron Pipeline Barplot Contig") +
  xlab("Pipeline steps") + ylab("Epitopes") +
  scale_fill_hue(l=40)
bar_plot_contig

## EXPORT
ggsave(filename="J_plots/XXX_barplot_contig.pdf", plot=bar_plot_contig, width=12, height=8, units="in", bg = "white")
