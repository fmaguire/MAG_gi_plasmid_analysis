library(tidyverse)
library(ggplot2)
library(stringr)

#Overall:
#  - how many VF genes are in each bin? how many are in each reference genome? 
#  - can we map the VF genes in each bin back to the reference genome?
#  - represent this as VF in bin normalized to VF in reference genome?
#  - find all the VF genes on plasmids/GI, where are they in the bins?

#step1: mapping files that maps:
# - reference chromosome and plasmid to bioaccession
# - bioaccession to species
# - bin to reference chromosome and plasmid
# - bins to tool and binner.

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
vfDF <- read_tsv("H:\\OneDrive\\School_Work\\_MAGAssessment\\newest\\RScripts\\vf_all0.001.tsv", col_names=TRUE)
#View(vfDF)

#filter rgiDF to have only dastool + megahit results
dastoolMegahitBinVFesults <- vfDF %>% filter(BinningTool == "dastool" & AssemblyTool == "megahit") %>% mutate(binID.t = str_remove(binID,".fa.vf"))
plasmidReferenceVFResults <- vfDF %>% filter(BinningTool == "plasmid" & AssemblyTool == "plasmid") %>% mutate(binID.t = str_remove(binID,".fasta.vf"))
chromosomeReferenceVFResults <- vfDF %>% filter(BinningTool == "chromosome" & AssemblyTool == "chromosome") %>% mutate(binID.t = str_remove(binID,".fasta.vf"))


#count total AMR per bin and per reference replicon
dastoolMegahitBinVFSummary <- dastoolMegahitBinVFesults %>% group_by(binID.t) %>% summarise(total = n()) %>% mutate (type="MAG")
plasmidReferenceVFSummary <- plasmidReferenceVFResults %>% group_by(binID.t) %>% summarise(total = n()) %>% mutate (type="Ref_Plasmid")
chromosomeReferenceVFSummary <- chromosomeReferenceVFResults %>% group_by(binID.t) %>% summarise(total = n()) %>% mutate (type="Ref_Chromosome")

plasmidReferenceVFSummary <- chromosomeReferenceVFSummary %>% mutate(total=0) %>% mutate(type="Ref_Plasmid")

#View(dastoolMegahitBinVFSummary)
#View(plasmidReferenceVFSummary)
#View(chromosomeReferenceVFSummary)


totalVFCount <- dastoolMegahitBinVFSummary %>% rbind(plasmidReferenceVFSummary) %>% rbind(chromosomeReferenceVFSummary) %>% filter(binID.t!="unbinned")
#View(totalVFCount)

#graph the counts
p <- ggplot(totalVFCount, aes(x=type, y=total)) +
  ggtitle("Total Number Of Virulence Factor Genes per type") +
  ylab("Number Of VF Genes") + 
  geom_boxplot(outlier.colour = 'black', outlier.shape = 16, outlier.size=1, notch = FALSE) 
plot(p)
ggsave("reboot_total_VF.svg", plot = last_plot(), device = "svg", width = 5, height = 8)

#part 2: how good are the AMR genes recovered when we map each bin to their reference chromosome
#normalize count to reference chromosome and plasmid
dastoolMegahitBinVFSummary.mapped <- dastoolMegahitBinVFSummary %>% select(-type) %>%
  full_join(dastoolMegahitBinChromosomeMap, by = c("binID.t"="bin")) %>% rename(totalMAG=total) %>%
  full_join(chromosomeReferenceVFSummary, by = c("chrReference"="binID.t")) %>% select(-type) %>%
  full_join(chromosomeRefBioAccMap, by = c("chrReference"="genomeAccession")) %>% rename(totalChr=total) %>%
  full_join(plasmideRefBioAccMap, by = c("bioAccession"="bioAccession")) %>% rename(plasReference=genomeAccession) %>%
  full_join(plasmidReferenceVFSummary, by = c("plasReference"="binID.t")) %>% select(-type) %>% rename(totalPlas=total) %>% 
  mutate(totalMAG=replace(totalMAG, is.na(totalMAG),0)) %>%
  mutate(totalChr=replace(totalChr, is.na(totalChr),0)) %>%
  mutate(totalPlas=replace(totalPlas, is.na(totalPlas),0))

#View(dastoolMegahitBinVFSummary.mapped)

dastoolMegahitBinVFSummary.mapped.chromosome <- dastoolMegahitBinVFSummary.mapped %>% 
  select (bioAccession, binID.t, totalMAG, chrReference, totalChr) %>% unique() %>%
  filter (!is.na(binID.t) | totalChr > 0)
dastoolMegahitBinVFSummary.mapped.chromosome.normalized <- dastoolMegahitBinVFSummary.mapped.chromosome %>%
  mutate(percentage = totalMAG / totalChr) %>% 
  mutate(percentage=replace(percentage, percentage=="NaN",1)) %>%
  mutate(type="chromosome")
  

dastoolMegahitBinVFSummary.mapped.plasmid<- dastoolMegahitBinVFSummary.mapped %>% 
  select (bioAccession, binID.t, totalMAG, plasReference, totalPlas) %>% unique()%>%
  filter (totalPlas > 0)
dastoolMegahitBinVFSummary.mapped.plasmid.normalized <- dastoolMegahitBinVFSummary.mapped.plasmid %>%
  mutate(percentage = totalMAG / totalPlas)%>%
  mutate(type="plasmid")

dastoolMegahitBinVFSummary.normalizedCount <- dastoolMegahitBinVFSummary.mapped.chromosome.normalized %>% select (bioAccession,type,percentage) %>%
  rbind((dastoolMegahitBinVFSummary.mapped.plasmid.normalized %>% select(bioAccession,type,percentage))) %>% filter(percentage!="Inf")

  #full_join(dastoolMegahitBinPlasmidMap, by = c("binID.t"="bin")) %>%
  #full_join(plasmidReferenceAMRSummary.noNudge, by = c("plasReference"="binID.t")) %>% select(-type) 
  #dastoolMegahitBinAMRSummary.noNudge.mapped[is.na(dastoolMegahitBinAMRSummary.noNudge.mapped)] <- 0 
#View(dastoolMegahitBinVFSummary.normalizedCount)

p <- ggplot(dastoolMegahitBinVFSummary.normalizedCount, aes(x=type, y=percentage)) +
  ggtitle("Virulence Factor Percent Recovery Of Each Bin to Reference Genome") +
  ylab("% Recovered") + 
#  geom_violin() +
  geom_boxplot(outlier.colour = 'black', outlier.shape = 16, outlier.size=1, notch = FALSE) +
  geom_jitter(position=position_jitter(width=.1, height=0))
plot(p)
ggsave("reboot_normalized_vf.svg", plot = last_plot(), device = "svg", width = 5, height = 8)


#part 3: are genes on plasmids recovered in the bins?

#plasmidAROs <- plasmidReferenceAMRResults %>% filter(is.na(Nudged)) %>% select(Best_Hit_ARO) %>% unique() %>% mutate(refLocation="plasmid")
#MagAROs <- dastoolMegahitBinAMRResults %>% select(Best_Hit_ARO, binID.t) %>% rename(MAGLocation=binID.t)
#plasmidAROLocation <- left_join(plasmidAROs, MagAROs, by=c("Best_Hit_ARO"="Best_Hit_ARO")) %>% unique()  %>%
#  mutate (found=case_when(MAGLocation == "NA" ~ "NOT Found", MAGLocation == "unbinned" ~"Found in Unbinned", TRUE ~ "Found in bin")) %>% 
#  group_by(found) %>% summarize(foundCount=n()) %>%
#  rename (Location=found) %>% rename(count = foundCount) %>%
#  mutate(Source="plasmid") %>% mutate(total=sum(count)) %>%
#  mutate(relative =  100 * count / total)

colnames(chromosomeReferenceVFResults)

#View(chromosomeReferenceVFResults %>% select(sseqid)%>% unique() )

chromosomeVFHits <- chromosomeReferenceVFResults %>% select(sseqid) %>% unique() %>% mutate(refLocation="chromosome") 
#View(chromosomeVFHits)
MagVFHits <- dastoolMegahitBinVFesults %>% select(sseqid, binID.t) %>% rename(MAGLocation=binID.t)
#View(MagVFHits%>% select(sseqid)%>% unique())
chromosomeVFLocation <- left_join(chromosomeVFHits, MagVFHits, by=c("sseqid"="sseqid")) %>% unique() %>% 
  group_by(sseqid, refLocation) %>% summarise(MAGLocation=paste(MAGLocation, collapse=";")) %>%
  mutate (found=case_when(MAGLocation == "NA" ~ "NOT Found", MAGLocation == "unbinned" ~"Found in Unbinned", TRUE ~ "Found in bin")) %>% 
  group_by(found) %>% summarize(foundCount=n())%>%
  rename (Location=found) %>% rename(count = foundCount) %>%
  mutate(Source="chromosome")%>% mutate(total=sum(count))%>%
  mutate(relative = 100 * count / total)
#View(chromosomeVFLocation)

#combinedARO <- rbind(plasmidAROLocation, chromosomeAROLocation) %>% select(Source, Location, relative)
p <- ggplot(chromosomeVFLocation, aes(Source,relative)) +  
  ggtitle("Where did the Reference Genome Virulence Factors go?") +
  ylab("% of total VF genes") + 
  geom_bar(aes(fill=Location), stat = "identity") + theme(legend.position = "bottom")
plot(p)
ggsave("reboot_location_vf.svg", plot = last_plot(), device = "svg", width = 5, height = 8)

