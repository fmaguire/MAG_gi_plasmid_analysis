library(tidyverse)
library(ggplot2)
library(stringr)

plasmidRecoveryData <- read_csv("./data/plasmids.csv", col_names = TRUE) %>% select("binner","assembly","Category","% of Plasmids (>50% coverage)")%>% 
  `colnames<-`(c("Binning_Tool","Assembly_Tool","Category","Relative")) %>% mutate(Group = "Plasmid")

giRecoveryData <- read_csv("./data/gi_data.csv", col_names = TRUE) %>% select("binner","assembly","Category","% of GIs (>50% Coverage)")%>% 
  `colnames<-`(c("Binning_Tool","Assembly_Tool","Category","Relative")) %>% mutate(Group = "Genomic Island")

amrRecoveryData <- read_tsv("./data/amrRecoveryData.tsv", col_names = TRUE)%>% 
  `colnames<-`(c("Binning_Tool","Assembly_Tool", "Absolute", "Relative", "Category")) %>% select(Binning_Tool, Assembly_Tool, Category, Relative) %>% mutate(Group = "AMR Genes")

vfRecoveryData <- read_tsv("./data/vfRecovery.tsv", col_names = TRUE)%>% 
  `colnames<-`(c("Binning_Tool","Assembly_Tool", "Absolute", "Relative", "Category")) %>% select(Binning_Tool, Assembly_Tool, Category, Relative) %>% mutate(Group = "VF Genes")

combinedData <- plasmidRecoveryData %>% rbind(giRecoveryData) %>% rbind(amrRecoveryData) %>% rbind(vfRecoveryData)
write.table(combinedData, "rateOfLoss.tsv", sep = "\t", row.names = FALSE)
assembled <- combinedData %>% filter(Category == "Assembled") %>% select(-Category) %>% rename(Assembled = Relative)
binned <- combinedData %>% filter(Category == "Binned")%>% select(-Category) %>% rename(Binned = Relative)
correctlyBinned <- combinedData %>% filter(Category == "Correctly Binned")%>% select(-Category) %>% rename(Correctly_Binned = Relative)

combined <- assembled %>% 
  full_join(binned, by = c("Assembly_Tool" = "Assembly_Tool","Binning_Tool" = "Binning_Tool", "Group" = "Group" )) %>%
  full_join(correctlyBinned, by = c("Assembly_Tool" = "Assembly_Tool","Binning_Tool" = "Binning_Tool", "Group" = "Group"))

View(combined)
colnames(combinedData)

p <- ggplot(combinedData, aes(x=Category, y=Relative, group=Group)) +
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group)) +
  facet_grid(Binning_Tool~Assembly_Tool) +
  ylim(0,100)+
  theme_light() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

plot(p)
