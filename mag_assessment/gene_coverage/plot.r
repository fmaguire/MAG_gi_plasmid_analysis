library(tidyverse)
library(ggplot2)
library(stringr)

#function to correct the capitalization and spelling of tools
correctToolNames<- function(dataframe)
{
  #correct binners
  df <- dataframe %>% 
    mutate(Binning_Tool = str_replace(Binning_Tool, "concoct", "CONCOCT"))%>% 
    mutate(Binning_Tool = str_replace(Binning_Tool, "dastool", "DAS_Tool"))%>% 
    mutate(Binning_Tool = str_replace(Binning_Tool, "maxbin2", "MaxBin2"))%>% 
    mutate(Binning_Tool = str_replace(Binning_Tool, "metabat2", "MetaBAT2")) %>% 
    mutate(Binning_Tool = str_replace(Binning_Tool, "reference", "Reference"))
  
  df <- df %>% 
    mutate(Assembly_Tool = str_replace(Assembly_Tool, "idba_ud", "IDBA_UD"))%>% 
    mutate(Assembly_Tool = str_replace(Assembly_Tool, "megahit", "MEGAHIT"))%>% 
    mutate(Assembly_Tool = str_replace(Assembly_Tool, "metaspades", "metaSPAdes")) %>%
    mutate(Assembly_Tool = str_replace(Assembly_Tool, "reference", "Reference"))
  
  return(df)
}

#condensed function to draw bar plots to keep themes consistent
BuildBarPlot <- function(ggplotObject, fileName, colorSet = 1)
{
  if (colorSet == 1) {
    colors <- c("#3274a1", "#e1812c", "#3a923a")
  } else {
    colors <- c("#91bfdb", "#abdda4", "#fc8d59")
  }
  
  p <- ggplotObject +
    facet_grid(.~Assembly_Tool, labeller = label_both) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9))+
    ylim(0,100)+
    theme_light() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    scale_fill_manual(values=colors)
  ggsave(fileName, plot = p, device = "png", width = 9, height = 4)
}

#condensed function to draw box plots to keep themes consistent
BuildBoxPlot <- function(ggplotObject, fileName)
{
  p <- ggplotObject +
    geom_boxplot(outlier.colour = 'black', outlier.shape = 16, outlier.size=1, notch = FALSE) +
    facet_grid(.~Assembly_Tool, labeller = label_both, scales="free") +
    theme_light() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) 
  ggsave(fileName, plot = p, device = "png", width = 9, height = 4)
}

#condensed function to draw point plots to keep themes consistent
BuildPointPlot<- function(ggplotObject, fileName)
{
  p<- ggplotObject +
    geom_point(aes(shape = Assembly_Tool, color = Binning_Tool), size = 1.5, stroke=1) +
    coord_flip() +
    theme_light() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  ggsave(fileName, plot = p, device = "png", width = 9, height = 4)
}

BuildLinePlot<- function(ggplotObject, fileName)
{
  p <- ggplotObject +
    geom_line(aes(color=Group))+
    geom_point(aes(color=Group)) +
    facet_grid(Binning_Tool~Assembly_Tool) +
    ylim(0,100)+
    theme_light() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  ggsave(fileName, plot = p, device = "png", width = 9, height = 7)
}


#read in all the tsv files corresponding to each figure
topHitsPerBinData <- read_csv("./data/top_hits_per_bin.csv", col_names = TRUE) %>% select("species","binner","assembly","coverage") %>% 
  `colnames<-` (c("Species","Binning_Tool","Assembly_Tool","Coverage")) %>% correctToolNames() %>% group_by(Species) %>%
  mutate(Average = mean(Coverage)) %>% ungroup() %>% mutate(Species=fct_reorder(Species, Average)) #add an average column to reorder the species to highest coverage to lowest

binCoverageData <- read_csv("./data/top_hits_per_bin.csv", col_names = TRUE) %>% select("binner","assembly","coverage")%>% 
`colnames<-`(c("Binning_Tool","Assembly_Tool","Coverage")) %>% correctToolNames()

binPurityData <- read_csv("./data/purity.csv", col_names = TRUE) %>% select("binner","assembly","species")%>% 
  `colnames<-`(c("Binning_Tool","Assembly_Tool","Num_Species")) %>% correctToolNames()

plasmidRecoveryData <- read_csv("./data/plasmids.csv", col_names = TRUE) %>% select("binner","assembly","Category","% of Plasmids (>50% coverage)")%>% 
  `colnames<-`(c("Binning_Tool","Assembly_Tool","Category","Percent")) %>% correctToolNames()

giRecoveryData <- read_csv("./data/gi_data.csv", col_names = TRUE) %>% select("binner","assembly","Category","% of GIs (>50% Coverage)")%>% 
  `colnames<-`(c("Binning_Tool","Assembly_Tool","Category","Percent")) %>% correctToolNames()

predictedGenesData <- read_csv("./data/geneCounts.csv", col_names = c("bin", "count", "type","binner","assembler")) %>% 
  `colnames<-`(c("Bin", "Count", "Type", "Binning_Tool","Assembly_Tool")) %>% correctToolNames() %>%
  filter(!grepl('unbinned', Bin)) #remove unbinned rows

amrRecoveryData <- read_tsv("./data/amrRecoveryData.tsv", col_names = TRUE)%>% 
  `colnames<-`(c("Binning_Tool","Assembly_Tool", "Absolute", "Relative", "Category")) %>% correctToolNames()

amrLocalizationRecoveryData <- read_tsv("./data/amrLocalizationRecoveryData.tsv", col_names = TRUE)%>% 
  `colnames<-`(c("Binning_Tool","Assembly_Tool", "Absolute", "Relative", "Category")) %>% correctToolNames()

vfRecoveryData <- read_tsv("./data/vfRecovery.tsv", col_names = TRUE)%>% 
  `colnames<-`(c("Binning_Tool","Assembly_Tool", "Absolute", "Relative", "Category")) %>% correctToolNames()

vfLocalizationRecoveryData <- read_tsv("./data/vfLocalizationRecovery.tsv", col_names = TRUE) %>% 
  `colnames<-`(c("Binning_Tool","Assembly_Tool", "Absolute", "Relative", "Category")) %>% correctToolNames()

rateOfLossData <- read_tsv("./data/rateOfLoss.tsv", col_names = TRUE) %>% correctToolNames()
View(rateOfLossData)

#draw the plots from data
topHitsPerBinPlot <- (ggplot(topHitsPerBinData, aes(x=Species, y=Coverage)) + ylab("Percent Coverage of Top Hit Genome") + xlab("Source Genome")) %>% 
  BuildPointPlot("top_hits_per_bin.png")

binCoveragePlot <- (ggplot(binCoverageData, aes(x=Binning_Tool, y=Coverage, fill = Binning_Tool)) + ylab("Percent Coverage of Top Hit Genome") + xlab("Binner")) %>% 
  BuildBoxPlot("bin_coverage.png")

binPurityPlot <- (ggplot(binPurityData, aes(x=Binning_Tool, y=Num_Species, fill = Binning_Tool)) + ylab("Species Per Bin (>5% Coverage)") + xlab("Binner")) %>%
  BuildBoxPlot("bin_purity.png")

plasmidRecoveryPlot <- (ggplot(plasmidRecoveryData, aes(x=Binning_Tool, y=Percent, fill=Category)) + ylab("Percent of Plasmids Recovered (>50% Coverage)") + xlab("Binner")) %>% 
  BuildBarPlot("plasmid_recovery.png")

giRecoveryPlot <- (ggplot(giRecoveryData, aes(x=Binning_Tool, y=Percent, fill=Category)) + ylab("Percent of GIs Recovered (>50% Coverage)") + xlab("Binner")) %>% 
  BuildBarPlot("GI_recovery.png")

predictedGenesPlot <- (ggplot(predictedGenesData, aes(x=Binning_Tool, y=Count, fill=Binning_Tool)) + ylab("Number of predicted Genes Per Bin") + xlab("Binner")) %>% 
  BuildBoxPlot("number_of_predicted_genes.png")

amrRecoveryPlot <- (ggplot(amrRecoveryData, aes(x=Binning_Tool, y=Relative, fill=Category)) + ylab("Proportion of Reference AMR genes") + xlab("Binner")) %>% 
  BuildBarPlot("amr_recovery.png")

amrLocalizationRecoveryPlot <- (ggplot(amrLocalizationRecoveryData, aes(x=Binning_Tool, y=Relative, fill=Category)) + ylab("Proportion of Reference AMR genes") + xlab("Binner")) %>% 
  BuildBarPlot("amr_localization_recovery.png",2)

vfRecoveryPlot <- (ggplot(vfRecoveryData, aes(x=Binning_Tool, y=Relative, fill=Category)) + ylab("Proportion of Reference VF genes") + xlab("Binner")) %>% 
  BuildBarPlot("vf_recovery.png")

vfLocalizationRecoveryPlot <- (ggplot(vfLocalizationRecoveryData, aes(x=Binning_Tool, y=Relative, fill=Category)) + ylab("Proportion of Reference VF Genes") + xlab("Binner")) %>% 
  BuildBarPlot("vf_localization_recovery.png",2)

rateOfLostPlot <- (ggplot(rateOfLossData, aes(x=Category, y=Relative, group=Group)) + ylab("Percent Recovery") + xlab("Stage")) %>%
  BuildLinePlot("rate_of_loss.png")
