library(tidyverse)
library(ggplot2)
library(stringr)

#Overall:
#  - how many AMR genes are in each bin? how many are in each reference genome? 
#  - can we map the AMR genes in each bin back to the reference genome?
#  - represent this as AMR in bin normalized to AMR in reference genome?
#  - find all the AMR genes on plasmids/GI, where are they in the bins?

#step1: mapping files that maps:
# - reference chromosome and plasmid to bioaccession
# - bioaccession to species
# - bin to reference chromosome and plasmid
# - bins to tool and binner.
MagAssessment.AMR("H:/OneDrive/School_Work/_MAGAssessment/newest/RScripts/rgi_all.tsv")
MagAssessment("H:/OneDrive/School_Work/_MAGAssessment/newest/RScripts/vf_all.tsv")


MagAssessment.AMR <- function(filePath)
{
  taxaMetadata <- read_tsv("H:/OneDrive/School_Work/_MAGAssessment/newest/MAGAssessment/taxa_metadata.tsv", col_names=c("spp","bioAccession", "genomeAccession", "type"))
  
  #reference - accession map
  chromosomeRefBioAccMap <- taxaMetadata %>% filter(type=="chromosome") %>% select(genomeAccession, bioAccession)
  plasmideRefBioAccMap <- taxaMetadata %>% filter(type=="plasmid") %>% select(genomeAccession, bioAccession)
  
  #accession - spp map
  bioAccSpeciesMap <- taxaMetadata %>% select(bioAccession, spp) %>% unique()
  
  #binToReferenceMap <-
  binChromosomeRef <- read_csv("H:\\OneDrive\\School_Work\\_MAGAssessment\\newest\\binning\\binChrRefMap.txt", col_names=c("binner","assembler", "bin", "chrReference"))
  dastoolMegahitBinChromosomeMap <- binChromosomeRef %>% filter(binner=="dastool" & assembler=="megahit") %>% mutate(bin=str_remove(bin,".fa")) %>% 
    mutate(chrReference = str_remove(chrReference,".fasta")) %>% select(bin, chrReference) 
  
  binPlasmidRef <- read_csv("H:\\OneDrive\\School_Work\\_MAGAssessment\\newest\\binning\\binPlasmidRefMap.txt", col_names=c("bin", "plasReference"))
  dastoolMegahitBinPlasmidMap <- binPlasmidRef %>% mutate(bin=str_remove(bin,".fa")) %>% mutate(plasReference = str_remove(plasReference,".fasta"))
  
  
  #first question, how many AMR genes are in each bin? how many are in each reference genome?
  #rgiDF <- read_tsv("H:/OneDrive/School_Work/_MAGAssessment/newest/RScripts/rgi_all.tsv", col_names=TRUE)
  rgiDF <- read_tsv(filePath, col_names=TRUE)
  
  #filter rgiDF to have only dastool + megahit results
  dastoolMegahitBinAMRResults <- rgiDF %>% filter(BinningTool == "dastool" & AssemblyTool == "megahit") %>% mutate(binID.t = str_remove(binID,".fa.txt"))
  plasmidReferenceAMRResults <- rgiDF %>% filter(BinningTool == "plasmid" & AssemblyTool == "plasmid") %>% mutate(binID.t = str_remove(binID,".fasta.txt"))
  chromosomeReferenceAMRResults <- rgiDF %>% filter(BinningTool == "chromosome" & AssemblyTool == "chromosome") %>% mutate(binID.t = str_remove(binID,".fasta.txt"))
  
  #count total AMR per bin and per reference replicon
  dastoolMegahitBinAMRSummary <- dastoolMegahitBinAMRResults %>% group_by(binID.t) %>% summarise(total = n())
  plasmidReferenceAMRSummary <- plasmidReferenceAMRResults %>% group_by(binID.t) %>% summarise(total = n())
  chromosomeReferenceAMRSummary <- chromosomeReferenceAMRResults %>% group_by(binID.t) %>% summarise(total = n())
  
  #count totalAMR thats not nudged
  
  dastoolMegahitBinAMRSummary.noNudge <- dastoolMegahitBinAMRResults %>% filter(is.na(Nudged)) %>% group_by(binID.t) %>% summarise(total = n()) %>% mutate (type="MAG")
  plasmidReferenceAMRSummary.noNudge <- plasmidReferenceAMRResults %>% filter(is.na(Nudged)) %>% group_by(binID.t) %>% summarise(total = n()) %>% mutate (type="Ref_Plasmid")
  chromosomeReferenceAMRSummary.noNudge <- chromosomeReferenceAMRResults %>% filter(is.na(Nudged)) %>% group_by(binID.t) %>% summarise(total = n()) %>% mutate (type="Ref_Chromosome")
  
  colnames(dastoolMegahitBinAMRSummary.noNudge)
  colnames(plasmidReferenceAMRSummary.noNudge)
  colnames(chromosomeReferenceAMRSummary.noNudge)
  
  
  totalAMRCount.noNudge <- dastoolMegahitBinAMRSummary.noNudge %>% rbind(plasmidReferenceAMRSummary.noNudge) %>% rbind(chromosomeReferenceAMRSummary.noNudge)
  
  #graph the counts
  p <- ggplot(totalAMRCount.noNudge, aes(x=type, y=total)) +
    ggtitle("Total Number Of AMR Genes per type") +
    ylab("Number Of AMR Genes") + 
    geom_boxplot(outlier.colour = 'black', outlier.shape = 16, outlier.size=1, notch = FALSE) 
  plot(p)
  ggsave("reboot_total_amr_noNudge.svg", plot = last_plot(), device = "svg", width = 5, height = 8)
  
  
  
  #part 2: how good are the AMR genes recovered when we map each bin to their reference chromosome
  #normalize count to reference chromosome and plasmid
  dastoolMegahitBinAMRSummary.noNudge.mapped <- dastoolMegahitBinAMRSummary.noNudge %>% select(-type) %>%
    full_join(dastoolMegahitBinChromosomeMap, by = c("binID.t"="bin")) %>% rename(totalMAG=total) %>%
    full_join(chromosomeReferenceAMRSummary.noNudge, by = c("chrReference"="binID.t")) %>% select(-type) %>%
    full_join(chromosomeRefBioAccMap, by = c("chrReference"="genomeAccession")) %>% rename(totalChr=total) %>%
    full_join(plasmideRefBioAccMap, by = c("bioAccession"="bioAccession")) %>% rename(plasReference=genomeAccession) %>%
    full_join(plasmidReferenceAMRSummary.noNudge, by = c("plasReference"="binID.t")) %>% select(-type) %>% rename(totalPlas=total) %>% 
    mutate(totalMAG=replace(totalMAG, is.na(totalMAG),0)) %>%
    mutate(totalChr=replace(totalChr, is.na(totalChr),0)) %>%
    mutate(totalPlas=replace(totalPlas, is.na(totalPlas),0))
  
  dastoolMegahitBinAMRSummary.noNudge.mapped.chromosome <- dastoolMegahitBinAMRSummary.noNudge.mapped %>% 
    select (bioAccession, binID.t, totalMAG, chrReference, totalChr) %>% unique() %>%
    filter (!is.na(binID.t) | totalChr > 0)
  dastoolMegahitBinAMRSummary.noNudge.mapped.chromosome.normalized <- dastoolMegahitBinAMRSummary.noNudge.mapped.chromosome %>%
    mutate(percentage = totalMAG / totalChr) %>% 
    mutate(percentage=replace(percentage, percentage=="NaN",1)) %>%
    mutate(type="chromosome")
  
  
  dastoolMegahitBinAMRSummary.noNudge.mapped.plasmid<- dastoolMegahitBinAMRSummary.noNudge.mapped %>% 
    select (bioAccession, binID.t, totalMAG, plasReference, totalPlas) %>% unique()%>%
    filter (totalPlas > 0)
  dastoolMegahitBinAMRSummary.noNudge.mapped.plasmid.normalized <- dastoolMegahitBinAMRSummary.noNudge.mapped.plasmid %>%
    mutate(percentage = totalMAG / totalPlas)%>%
    mutate(type="plasmid")
  
  dastoolMegahitBinAMRSummary.noNudge.normalizedCount <- dastoolMegahitBinAMRSummary.noNudge.mapped.chromosome.normalized %>% select (bioAccession,type,percentage) %>%
    rbind((dastoolMegahitBinAMRSummary.noNudge.mapped.plasmid.normalized %>% select(bioAccession,type,percentage))) %>% filter(percentage!="Inf")
  
  #full_join(dastoolMegahitBinPlasmidMap, by = c("binID.t"="bin")) %>%
  #full_join(plasmidReferenceAMRSummary.noNudge, by = c("plasReference"="binID.t")) %>% select(-type) 
  #dastoolMegahitBinAMRSummary.noNudge.mapped[is.na(dastoolMegahitBinAMRSummary.noNudge.mapped)] <- 0 
  
  p <- ggplot(dastoolMegahitBinAMRSummary.noNudge.normalizedCount, aes(x=type, y=percentage)) +
    ggtitle("AMR Gene Percent Recovery Of Each Bin to Reference Genome") +
    ylab("% Recovered") + 
    geom_boxplot(outlier.colour = 'black', outlier.shape = 16, outlier.size=1, notch = FALSE) +
    geom_jitter(position=position_jitter(width=.1, height=0))
  plot(p)
  ggsave("reboot_normalized_amr_noNudge.svg", plot = last_plot(), device = "svg", width = 5, height = 8)
  
  
  #part 3: are genes on plasmids recovered in the bins?
  
  plasmidAROs <- plasmidReferenceAMRResults %>% filter(is.na(Nudged)) %>% select(Best_Hit_ARO) %>% unique() %>% mutate(refLocation="plasmid")
  MagAROs <- dastoolMegahitBinAMRResults %>% select(Best_Hit_ARO, binID.t) %>% rename(MAGLocation=binID.t)
  plasmidAROLocation <- left_join(plasmidAROs, MagAROs, by=c("Best_Hit_ARO"="Best_Hit_ARO")) %>% unique()  %>%
    mutate (found=case_when(MAGLocation == "NA" ~ "NOT Found", MAGLocation == "unbinned" ~"Found in Unbinned", TRUE ~ "Found in bin")) %>% 
    group_by(found) %>% summarize(foundCount=n()) %>%
    rename (Location=found) %>% rename(count = foundCount) %>%
    mutate(Source="plasmid") %>% mutate(total=sum(count)) %>%
    mutate(relative =  100 * count / total)
  
  chromosomeAROs <- chromosomeReferenceAMRResults %>% filter(is.na(Nudged)) %>% select(Best_Hit_ARO) %>% unique() %>% mutate(refLocation="chromosome") 
  MagAROs <- dastoolMegahitBinAMRResults %>% select(Best_Hit_ARO, binID.t) %>% rename(MAGLocation=binID.t)
  chromosomeAROLocation <- left_join(chromosomeAROs, MagAROs, by=c("Best_Hit_ARO"="Best_Hit_ARO")) %>% unique() %>% 
    group_by(Best_Hit_ARO, refLocation) %>% summarise(MAGLocation=paste(MAGLocation, collapse=";")) %>%
    mutate (found=case_when(MAGLocation == "NA" ~ "NOT Found", MAGLocation == "unbinned" ~"Found in Unbinned", TRUE ~ "Found in bin")) %>% 
    group_by(found) %>% summarize(foundCount=n())%>%
    rename (Location=found) %>% rename(count = foundCount) %>%
    mutate(Source="chromosome")%>% mutate(total=sum(count))%>%
    mutate(relative = 100 * count / total)
  
  combinedARO <- rbind(plasmidAROLocation, chromosomeAROLocation) %>% select(Source, Location, relative)
  p <- ggplot(combinedARO, aes(Source,relative)) +  
    ggtitle("Where did the Reference Genome AMR Genes go?") +
    ylab("% of total AMR genes") + 
    geom_bar(aes(fill=Location), stat = "identity") + theme(legend.position = "bottom")
  plot(p)
  ggsave("reboot_location_amr_noNudge.svg", plot = last_plot(), device = "svg", width = 5, height = 8)
  
}
