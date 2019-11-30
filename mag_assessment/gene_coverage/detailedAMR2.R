library(tidyverse)
library(ggplot2)
library(stringr)

#1. AMR and Virulence factor recovery:
#  a. facet = assembly
#  b. x-axis = binner
#  c. y-axis = %recovery
#  d. category = assembled (category_assembly), binned, correctly binned. 


taxaMetadata <- read_tsv("H:/OneDrive/School_Work/_MAGAssessment/newest/MAGAssessment/taxa_metadata.tsv", col_names=c("spp","bioAccession", "genomeAccession", "type"))

#reference - accession map
chromosomeRefBioAccMap <- taxaMetadata %>% filter(type=="chromosome") %>% select(genomeAccession, bioAccession)
plasmideRefBioAccMap <- taxaMetadata %>% filter(type=="plasmid") %>% select(genomeAccession, bioAccession)
plasmideRefBioAccMap <- taxaMetadata %>% filter(type=="plasmid") %>% select(genomeAccession, bioAccession)


#accession - spp map
bioAccSpeciesMap <- taxaMetadata %>% select(bioAccession, spp) %>% unique()

#binToReferenceMap <-
binChromosomeRef <- read_csv("H:\\OneDrive\\School_Work\\_MAGAssessment\\newest\\binning\\last_chrRefMap.txt", col_names=c("binner","assembler", "bin", "chrReference"))%>% 
  mutate(bin=str_remove(bin,".tabr")) %>% 
  mutate(bin=str_remove(bin,".fasta")) %>% 
  mutate(bin=str_remove(bin,".fa")) %>% 
  mutate(chrReference = str_remove(chrReference,".fasta"))

#dastoolMegahitBinChromosomeMap <- binChromosomeRef %>% filter(binner=="dastool" & assembler=="megahit") %>% mutate(bin=str_remove(bin,".fa")) %>% 
#  mutate(chrReference = str_remove(chrReference,".fasta")) %>% select(bin, chrReference) 

binPlasmidRef <- read_csv("H:\\OneDrive\\School_Work\\_MAGAssessment\\newest\\binning\\last_plasRefMap.txt", col_names=c("binner","assembler", "bin", "plasReference")) %>% 
  mutate(bin=str_remove(bin,".tabp")) %>% 
  mutate(bin=str_remove(bin,".fasta")) %>% 
  mutate(bin=str_remove(bin,".fa")) %>% 
  mutate(plasReference = str_remove(plasReference,".fasta"))

binGIRef <- read_csv("H:\\OneDrive\\School_Work\\_MAGAssessment\\newest\\binning\\last_giRefMap.txt", col_names=c("binner","assembler", "bin", "plasReference")) %>% 
  mutate(bin=str_remove(bin,".tabg")) %>% 
  mutate(bin=str_remove(bin,".fasta")) %>% 
  mutate(bin=str_remove(bin,".fa")) %>% 
  mutate(plasReference = str_remove(plasReference,".fasta"))
#dastoolMegahitBinPlasmidMap <- binPlasmidRef %>% mutate(bin=str_remove(bin,".fa")) %>% mutate(plasReference = str_remove(plasReference,".fasta"))


#first question, how many AMR genes are in each bin? how many are in each reference genome?
rgiDF <- read_tsv("H:/OneDrive/School_Work/_MAGAssessment/newest/RScripts/results/rgi.tsv", col_names=TRUE)
#rgiDF <- read_tsv(filePath, col_names=TRUE)

#filter rgiDF to have only dastool + megahit results
magBinAMRResults <- rgiDF %>% filter(Type=="mag") %>% mutate(binID.t = str_remove(binID,".txt")) %>% mutate(binID.t = str_remove(binID.t,".fasta"))%>% mutate(binID.t = str_remove(binID.t,".fa"))
plasmidReferenceAMRResults <- rgiDF %>% filter(BinningTool == "plasmid" & AssemblyTool == "plasmid") %>% mutate(binID.t = str_remove(binID,".fasta.txt"))
chromosomeReferenceAMRResults <- rgiDF %>% filter(BinningTool == "chromosome" & AssemblyTool == "chromosome") %>% mutate(binID.t = str_remove(binID,".fasta.txt"))
giReferenceAMRResults <- rgiDF %>% filter(BinningTool == "gi" & AssemblyTool == "gi") %>% mutate(binID.t = str_remove(binID,".txt"))
assemblyAMRResults <- rgiDF %>% filter (Type =="assembly") %>% mutate(binID.t = str_remove(binID,".txt")) %>% mutate(binID.t = str_remove(binID.t,".fasta"))%>% mutate(binID.t = str_remove(binID.t,".fa"))
View(assemblyAMRResults)

#count total AMR per bin and per reference replicon
magBinAMRSummary <- magBinAMRResults %>% group_by(BinningTool,AssemblyTool, binID.t) %>% summarise(total = n()) %>% mutate (type="MAG")
plasmidReferenceAMRSummary <- plasmidReferenceAMRResults %>% group_by(binID.t) %>% summarise(total = n()) %>% mutate (type="Ref_Plasmid")
chromosomeReferenceAMRSummary <- chromosomeReferenceAMRResults %>% group_by(binID.t) %>% summarise(total = n()) %>% mutate (type="Ref_Chromosome")
giReferenceAMRSummary <- giReferenceAMRResults %>% group_by(binID.t) %>% summarise(total = n()) %>% mutate (type="Ref_GI")
assemblyAMRSummary <- assemblyAMRResults %>% group_by(binID.t) %>% summarise(total = n()) %>% mutate (type="Assembly")

#count totalAMR thats not nudged
magBinAMRSummary.noNudge <- magBinAMRResults %>% filter(is.na(Nudged)) %>% group_by(BinningTool,AssemblyTool, binID.t) %>% summarise(total = n()) %>% mutate (type="MAG")
plasmidReferenceAMRSummary.noNudge <- plasmidReferenceAMRResults %>% filter(is.na(Nudged)) %>% group_by(BinningTool,AssemblyTool,binID.t) %>% summarise(total = n()) %>% mutate (type="Ref_Plasmid")
chromosomeReferenceAMRSummary.noNudge <- chromosomeReferenceAMRResults %>% filter(is.na(Nudged)) %>% group_by(BinningTool,AssemblyTool,binID.t) %>% summarise(total = n()) %>% mutate (type="Ref_Chromosome")
giReferenceAMRSummary.noNudge <- giReferenceAMRResults %>% filter(is.na(Nudged)) %>% group_by(BinningTool,AssemblyTool,binID.t) %>% summarise(total = n()) %>% mutate (type="Ref_GI")
assemblyAMRSummary.noNudge <-  assemblyAMRResults %>% filter(is.na(Nudged)) %>% group_by(BinningTool,AssemblyTool,binID.t) %>% summarise(total = n()) %>% mutate (type="Assembly")

totalAMRCount.noNudge <- magBinAMRSummary.noNudge %>% rbind(plasmidReferenceAMRSummary.noNudge) %>% rbind(chromosomeReferenceAMRSummary.noNudge)%>% rbind(giReferenceAMRSummary.noNudge) %>% rbind(assemblyAMRSummary.noNudge)

#part 2: how good are the AMR genes recovered when we map each bin to their reference chromosome
#normalize count to reference chromosome and plasmid
 
#total assembled
totalAMRCount.noNudge.all <- totalAMRCount.noNudge %>% group_by(type, AssemblyTool, BinningTool) %>% mutate(assembled.n=sum(total)) %>% select(-binID.t, -total) %>% unique() 
refAMR.noNudge <- sum((totalAMRCount.noNudge.all %>% filter(type=="Ref_Plasmid" || type=="Ref_Chromosome"))$assembled.n)

totalAMRCount.noNudge.assembled <- totalAMRCount.noNudge %>% filter(type=="Assembly") %>% group_by(type, AssemblyTool, BinningTool) %>% mutate(assembled.n=sum(total)) %>% select(-binID.t, -total) %>% unique() 
category_assembly <- totalAMRCount.noNudge.assembled %>% mutate(assembled.p=100*(assembled.n)/refAMR.noNudge) %>% ungroup() %>% select(-type) %>% mutate(Category="Assembled")

#View(refAMR.noNudge)

#total binned
totalAMRCount.noNudge.binned <- totalAMRCount.noNudge %>% filter(!grepl('unbinned', binID.t))%>% group_by(type, AssemblyTool, BinningTool)%>% mutate(binned.n=sum(total)-4) %>% select(-binID.t, -total) %>% unique() 
totalAMRCount.noNudge.unbinned <- totalAMRCount.noNudge %>% filter(grepl('unbinned', binID.t)) %>% group_by(type, AssemblyTool, BinningTool)%>% mutate(binned.n=sum(total)) %>% select(-binID.t, -total) %>% unique()
category_binned <- totalAMRCount.noNudge.binned %>% filter(type=="MAG") %>% mutate(binned.p=100*binned.n/refAMR.noNudge) %>% ungroup() %>% select(-type) %>% mutate(Category="Binned")
#View(totalAMRCount.noNudge.binned)

magBinAMRResults.short <- magBinAMRResults %>% filter(!grepl('unbinned', binID.t)) %>% 
  filter(is.na(Nudged)) %>% select(ARO, BinningTool,AssemblyTool,binID.t)
plasmidReferenceAMRResults.short <- plasmidReferenceAMRResults %>% 
  filter(is.na(Nudged)) %>% select(ARO,binID.t) %>% rename(ARO.Plas=ARO)
chromosomeReferenceAMRResults.short <- chromosomeReferenceAMRResults %>% 
  filter(is.na(Nudged))%>% select(ARO, binID.t)  %>% rename(ARO.Chr=ARO)
giReferenceAMRResults.short <- giReferenceAMRResults %>% 
  filter(is.na(Nudged))%>% select(ARO, binID.t)  %>% rename(ARO.GI=ARO)

magBinAMRResults.short.mapped <- magBinAMRResults.short %>% 
  full_join(binChromosomeRef, by = c("BinningTool"="binner","AssemblyTool"="assembler","binID.t"="bin"))

ChrPlasMap <- chromosomeRefBioAccMap %>% full_join(plasmideRefBioAccMap, by = c("bioAccession"="bioAccession")) %>% select(-bioAccession) %>% rename(chr=genomeAccession.x)%>% rename(plas=genomeAccession.y)
plasmidReferenceAMRResults.short.mapped <- plasmidReferenceAMRResults.short %>% left_join(ChrPlasMap, by = c("binID.t" = "plas")) %>% rename(plas = binID.t)

#1. MAG ARO inc reference genome- magBinAMRResults.short.mapped
#2. Chr ARO - chromosomeReferenceAMRResults.short
#3. Plas ARO inc reference chromosome - plasmidReferenceAMRResults.short.mapped
#4. GI ARO inc reference chromosme - giReferenceAMRResults.short

mag <- magBinAMRResults.short.mapped
chr <- chromosomeReferenceAMRResults.short
plas <- plasmidReferenceAMRResults.short.mapped
gi <- giReferenceAMRResults.short

magxchr <- mag %>% left_join(chr, by = c("chrReference" = "binID.t")) %>% filter(!grepl('unbinned', binID.t))
magxplas <- mag %>% left_join(plas, by = c("chrReference" = "chr")) %>% filter(!grepl('unbinned', binID.t))
magxgi <- mag %>% left_join(gi, by = c("chrReference" = "binID.t")) %>% filter(!grepl('unbinned', binID.t))
View(magxgi.t)
magxchr.t <- magxchr %>% filter(ARO==ARO.Chr) %>% unique() %>% group_by(BinningTool, AssemblyTool) %>% summarise(correctlyBinned.n = n()) %>% mutate(correctlyBinned.p=100*correctlyBinned.n/refAMR.noNudge) %>% mutate(Category="Correctly Binned")
magxplas.t <- magxplas %>% filter(ARO==ARO.Plas) %>% unique() %>% group_by(BinningTool, AssemblyTool) %>% summarise(correctlyBinned.n = n()) %>% mutate(correctlyBinned.p=100*correctlyBinned.n/refAMR.noNudge) %>% mutate(Category="Correctly Binned")
magxgi.t <- magxgi %>% filter(ARO==ARO.GI) %>% unique() %>% group_by(BinningTool, AssemblyTool) %>% summarise(correctlyBinned.n = n()) %>% mutate(correctlyBinned.p=100*correctlyBinned.n/refAMR.noNudge) %>% mutate(Category="Correctly Binned")

category_correctlyBinned <- magxchr.t
colnames(category_assembly) <-colnames(category_binned) <-colnames(category_correctlyBinned) <- c("BinningTool", "AssemblyTool", "absolute","relative","Category")
firstFigureData <- category_assembly %>% rbind(as.data.frame(category_binned)) %>% rbind(as.data.frame(category_correctlyBinned))
View(secondFigureData)
p <- ggplot(firstFigureData, aes(x=BinningTool, y=relative, fill=Category)) +
  facet_grid(.~AssemblyTool, labeller = label_both) +
  ylab("Proportion of Reference AMR genes") + 
  xlab("Binner") +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))+
  ylim(0,100)+
  theme_light() +
  scale_fill_manual(values=c("#3274a1", "#e1812c", "#3a923a"))
plot(p)
ggsave("amr_recovery.png", plot = last_plot(), device = "png", width = 8.5, height = 3)



#2. all the amr gene thats been correctly binned. are they in the right context (unstacked), reflecting all 12 pairs.

totPlas <- (totalAMRCount.noNudge.all %>% filter(AssemblyTool=="plasmid"))$assembled.n
totGI <- (totalAMRCount.noNudge.all %>% filter(AssemblyTool=="gi"))$assembled.n
totChr <- (totalAMRCount.noNudge.all %>% filter(AssemblyTool=="chromosome"))$assembled.n - totGI

magxchr.2 <- magxchr.t %>% mutate(Category="Chromosome") %>% mutate(correctlyBinned.p=correctlyBinned.n * 100 / totChr)
magxplas.2 <- magxchr.t %>% mutate(correctlyBinned.n=0) %>% mutate(correctlyBinned.p=0) %>% mutate(Category="Plasmid")#empty frame hax #magxplas.t %>% mutate(Category="Plasmid") %>% mutate(correctlyBinned.p=correctlyBinned.n * 100 / totPlas)
tempDF<-data.frame("metabat2","idba_ud",as.numeric(0),as.numeric(0),"Genomic Island")
colnames(tempDF)<-colnames(magxplas.2)
magxgi.2 <- magxgi.t %>% mutate(Category="Genomic Island") %>% mutate(correctlyBinned.p=correctlyBinned.n * 100 / totGI) %>% as.data.frame() %>% rbind(as.data.frame(tempDF))

secondFigureData <- as.data.frame(magxchr.2) %>% rbind(as.data.frame(magxplas.2)) %>% rbind(as.data.frame(magxgi.2))
colnames(secondFigureData) <- c("BinningTool","AssemblyTool","absolute","relative","Category")
p <- ggplot(secondFigureData, aes(x=BinningTool, y=relative, fill=Category)) +
  facet_grid(.~AssemblyTool, labeller = label_both) +
  ylab("Proportion of Reference AMR genes in Each Category") + 
  xlab("Binner") +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))+
  ylim(0,100)+
  theme_light() +
  scale_fill_manual(values=c("#3274a1", "#e1812c", "#3a923a"))
plot(p)
ggsave("amr_localization_recovery.png", plot = last_plot(), device = "png", width = 8.5, height = 4)
