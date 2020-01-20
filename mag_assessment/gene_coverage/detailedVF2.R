library(tidyverse)
library(ggplot2)
library(stringr)

#1. VF and Virulence factor recovery:
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


#first question, how many VF genes are in each bin? how many are in each reference genome?
vfDF <- read_tsv("H:/OneDrive/School_Work/_MAGAssessment/newest/RScripts/results/vf_all0.001.tsv", col_names=TRUE)
#vfDF <- read_tsv(filePath, col_names=TRUE)

#filter vfDF to have only dastool + megahit results
magBinVFResults <- vfDF %>% filter(Type=="mag") %>% filter(evalue<=1e-30)  %>% filter(pident > 90) %>% mutate(binID.t = str_remove(binID,".vf")) %>% mutate(binID.t = str_remove(binID.t,".fasta"))%>% mutate(binID.t = str_remove(binID.t,".fa"))
plasmidReferenceVFResults <- vfDF %>% filter(BinningTool == "plasmid" & AssemblyTool == "plasmid") %>% mutate(binID.t = str_remove(binID,".fasta.vf"))
chromosomeReferenceVFResults <- vfDF %>% filter(BinningTool == "chromosome" & AssemblyTool == "chromosome") %>% mutate(binID.t = str_remove(binID,".fasta.vf"))
giReferenceVFResults <- vfDF %>% filter(BinningTool == "gi" & AssemblyTool == "gi") %>% mutate(binID.t = str_remove(binID,".vf"))
assemblyVFResults <- vfDF %>% filter (Type =="assembly") %>% filter(evalue<=1e-30)  %>% filter(pident > 90) %>% mutate(binID.t = str_remove(binID,".txt")) %>% mutate(binID.t = str_remove(binID.t,".fasta"))%>% mutate(binID.t = str_remove(binID.t,".fa"))
#View(assemblyVFResults)

magBinVFSummary <- magBinVFResults %>% select(sseqid, Type,	BinningTool,	AssemblyTool,	binID, binID.t) %>% unique() %>% group_by(BinningTool,AssemblyTool, binID.t) %>% summarise(total = n()) %>% mutate (type="MAG")
plasmidReferenceVFSummary <- plasmidReferenceVFResults  %>% select(sseqid,Type,	BinningTool,	AssemblyTool,	binID, binID.t) %>% unique()%>% group_by(BinningTool,AssemblyTool,binID.t) %>% summarise(total = n()) %>% mutate (type="Ref_Plasmid")
chromosomeReferenceVFSummary <- chromosomeReferenceVFResults  %>% select(sseqid,Type,	BinningTool,	AssemblyTool,	binID, binID.t) %>% unique()%>% group_by(BinningTool,AssemblyTool,binID.t) %>% summarise(total = n()) %>% mutate (type="Ref_Chromosome")
giReferenceVFSummary <- giReferenceVFResults %>% select(sseqid,Type,	BinningTool,	AssemblyTool,	binID, binID.t) %>% unique() %>% group_by(BinningTool,AssemblyTool,binID.t) %>% summarise(total = n()) %>% mutate (type="Ref_GI")
assemblyVFSummary <- assemblyVFResults %>% select(sseqid,Type,	BinningTool,	AssemblyTool,	binID, binID.t) %>% unique() %>% group_by(BinningTool,AssemblyTool,binID.t) %>% summarise(total = n()) %>% mutate (type="Assembly")

#View(totalVFCount)
totalVFCount <- magBinVFSummary %>% rbind(plasmidReferenceVFSummary) %>% rbind(chromosomeReferenceVFSummary)%>% rbind(giReferenceVFSummary) %>% rbind(assemblyVFSummary)

#part 2: how good are the VF genes recovered when we map each bin to their reference chromosome
#normalize count to reference chromosome and plasmid

#total assembled
totalVFCount.all <- totalVFCount %>% group_by(type, AssemblyTool, BinningTool) %>% mutate(assembled.n=sum(total)) %>% select(-binID.t, -total) %>% unique() 


totalVFCount.assembled <- totalVFCount %>% filter(type=="Assembly") %>% group_by(type, AssemblyTool, BinningTool) %>% mutate(assembled.n=sum(total)) %>% select(-binID.t, -total) %>% unique() 
refVF <- sum((totalVFCount.all %>% filter(type=="Ref_Plasmid" || type=="Ref_Chromosome"))$assembled.n)
category_assembly <- totalVFCount.assembled %>% mutate(assembled.p=100*(assembled.n)/refVF) %>% ungroup() %>% select(-type) %>% mutate(Category="Assembled")

#View(refVF)
#total binned
totalVFCount.binned <- totalVFCount %>% filter(!grepl('unbinned', binID.t))%>% group_by(type, AssemblyTool, BinningTool)%>% mutate(binned.n=sum(total)-4) %>% select(-binID.t, -total) %>% unique() 
totalVFCount.unbinned <- totalVFCount %>% filter(grepl('unbinned', binID.t)) %>% group_by(type, AssemblyTool, BinningTool)%>% mutate(binned.n=sum(total)) %>% select(-binID.t, -total) %>% unique()
category_binned <- totalVFCount.binned %>% filter(type=="MAG") %>% mutate(binned.p=100*binned.n/refVF) %>% ungroup() %>% select(-type) %>% mutate(Category="Binned")

magBinVFResults.short <- magBinVFResults %>% filter(!grepl('unbinned', binID.t)) %>% select(sseqid, BinningTool,AssemblyTool,binID.t)
plasmidReferenceVFResults.short <- plasmidReferenceVFResults %>% select(sseqid,binID.t) %>% rename(sseqid.Plas=sseqid)
chromosomeReferenceVFResults.short <- chromosomeReferenceVFResults %>% select(sseqid, binID.t)  %>% rename(sseqid.Chr=sseqid)
giReferenceVFResults.short <- giReferenceVFResults %>% select(sseqid, binID.t)  %>% rename(sseqid.GI=sseqid)

magBinVFResults.short.mapped <- magBinVFResults.short %>% 
  full_join(binChromosomeRef, by = c("BinningTool"="binner","AssemblyTool"="assembler","binID.t"="bin"))
ChrPlasMap <- chromosomeRefBioAccMap %>% full_join(plasmideRefBioAccMap, by = c("bioAccession"="bioAccession")) %>% select(-bioAccession) %>% rename(chr=genomeAccession.x)%>% rename(plas=genomeAccession.y)
plasmidReferenceVFResults.short.mapped <- plasmidReferenceVFResults.short %>% left_join(ChrPlasMap, by = c("binID.t" = "plas")) %>% rename(plas = binID.t)

#1. MAG sseqid inc reference genome- magBinVFResults.short.mapped
#2. Chr sseqid - chromosomeReferenceVFResults.short
#3. Plas sseqid inc reference chromosome - plasmidReferenceVFResults.short.mapped
#4. GI sseqid inc reference chromosme - giReferenceVFResults.short

mag <- magBinVFResults.short.mapped
chr <- chromosomeReferenceVFResults.short
plas <- plasmidReferenceVFResults.short.mapped
gi <- giReferenceVFResults.short


magxchr <- mag %>% left_join(chr, by = c("chrReference" = "binID.t")) %>% filter(!grepl('unbinned', binID.t))
magxplas <- mag %>% left_join(plas, by = c("chrReference" = "chr")) %>% filter(!grepl('unbinned', binID.t))
magxgi <- mag %>% left_join(gi, by = c("chrReference" = "binID.t")) %>% filter(!grepl('unbinned', binID.t))

magxchr.t <- magxchr %>% filter(sseqid==sseqid.Chr) %>% unique() %>% group_by(BinningTool, AssemblyTool) %>% summarise(correctlyBinned.n = n()) %>% mutate(correctlyBinned.p=100*correctlyBinned.n/refVF) %>% mutate(Category="Correctly Binned")
magxplas.t <- magxplas %>% filter(sseqid==sseqid.Plas) %>% unique() %>% group_by(BinningTool, AssemblyTool) %>% summarise(correctlyBinned.n = n()) %>% mutate(correctlyBinned.p=100*correctlyBinned.n/refVF) %>% mutate(Category="Correctly Binned")
magxgi.t <- magxgi %>% filter(sseqid==sseqid.GI) %>% unique() %>% group_by(BinningTool, AssemblyTool) %>% summarise(correctlyBinned.n = n()) %>% mutate(correctlyBinned.p=100*correctlyBinned.n/refVF) %>% mutate(Category="Correctly Binned")


category_correctlyBinned <- magxchr.t
colnames(category_assembly) <-colnames(category_binned) <-colnames(category_correctlyBinned) <- c("BinningTool", "AssemblyTool", "absolute","relative","Category")
firstFigureData <- category_assembly %>% rbind(as.data.frame(category_binned)) %>% rbind(as.data.frame(category_correctlyBinned))
View(secondFigureData)
p <- ggplot(firstFigureData, aes(x=BinningTool, y=relative, fill=Category)) +
  facet_grid(.~AssemblyTool, labeller = label_both) +
  ylab("Proportion of Reference VF genes") + 
  xlab("Binner") +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))+
  ylim(0,100)+
  theme_light() +
  scale_fill_manual(values=c("#3274a1", "#e1812c", "#3a923a"))
plot(p)
ggsave("VF_recovery.png", plot = last_plot(), device = "png", width = 8.5, height = 3)


 
#2. all the VF gene thats been correctly binned. are they in the right context (unstacked), reflecting all 12 pairs.

totPlas <- (totalVFCount.all %>% filter(AssemblyTool=="plasmid"))$assembled.n
totGI <- (totalVFCount.all %>% filter(AssemblyTool=="gi"))$assembled.n
totChr <- (totalVFCount.all %>% filter(AssemblyTool=="chromosome"))$assembled.n - totGI

magxchr.2 <- magxchr.t %>% mutate(Category="Chromosome") %>% mutate(correctlyBinned.p=correctlyBinned.n * 100 / totChr)
magxplas.2 <- magxchr.t %>% mutate(correctlyBinned.n=0) %>% mutate(correctlyBinned.p=0) %>% mutate(Category="Plasmid")#empty frame hax #magxplas.t %>% mutate(Category="Plasmid") %>% mutate(correctlyBinned.p=correctlyBinned.n * 100 / totPlas)
tempDF<-data.frame("metabat2","idba_ud",as.numeric(0),as.numeric(0),"Genomic Island")
colnames(tempDF)<-colnames(magxplas.2)
magxgi.2 <- magxgi.t %>% mutate(Category="Genomic Island") %>% mutate(correctlyBinned.p=correctlyBinned.n * 100 / totGI) %>% as.data.frame() %>% rbind(as.data.frame(tempDF))

secondFigureData <- as.data.frame(magxchr.2) %>% rbind(as.data.frame(magxplas.2)) %>% rbind(as.data.frame(magxgi.2))
colnames(secondFigureData) <- c("BinningTool","AssemblyTool","absolute","relative","Category")
#View(secondFigureData)
p <- ggplot(secondFigureData, aes(x=BinningTool, y=relative, fill=Category)) +
  facet_grid(.~AssemblyTool, labeller = label_both) +
  ylab("Proportion of Reference VF genes in Each Category") + 
  xlab("Binner") +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))+
  ylim(0,100)+
  theme_light() +
  scale_fill_manual(values=c("#3274a1", "#e1812c", "#3a923a"))
plot(p)
ggsave("VF_localization_recovery.png", plot = last_plot(), device = "png", width = 8.5, height = 4)

