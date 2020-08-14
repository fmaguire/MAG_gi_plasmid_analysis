library(tidyverse)
library(ggplot2)
library(stringr)

#script to generate boxplots of the total number of recovered AMR/VF genes from each binner/assembler pair. 

refType <-c("chromosome","plasmid")
binningTools <- c("concoct","dastool","maxbin2","metabat2")
assemblyTools <- c("idba_ud","megahit","metaspades")

#txt are RGI, vf are VFDB
#lets make 1 GIANT RGI/VFDB file
rgiDF = data.frame()

#parse the individual VF/AMR result files in results/$type/$binner/$assembler
t = "mag"
for (b in binningTools)
{
  for (a in assemblyTools)
  {
    #lets make 1 GIANT RGI/VFDB file
    for (value in list.files(path=paste0("results/",t,"/",b,"_results/",a,"_results/"), pattern="*.txt", full.names=TRUE, recursive=FALSE) %>% as.vector())
    {
      print(value)
      rgi <- read_tsv(value, col_names = dfColNames) %>% 
        filter(Cut_Off=="Strict" | Cut_Off =="Perfect") %>% #used for amr genes
        mutate(Type = t) %>% 
        mutate (BinningTool = b) %>% 
        mutate (AssemblyTool = a) %>% 
        mutate(binID=basename(value))
      rgiDF <- rbind(rgiDF, rgi)
    }
  }
}
t="references"
for (r in refType)
{
    #lets make 1 GIANT RGI/VFDB file
    for (value in list.files(path=paste0("results/references/",r,"_results/"), pattern="*.vf", full.names=TRUE, recursive=FALSE) %>% as.vector())
    {
      print(value)
      rgi <- read_tsv(value, col_names = dfColNames) %>% 
        filter(Cut_Off=="Strict" | Cut_Off =="Perfect") %>% #used for amr genes
        mutate(Type = t) %>% 
        mutate (BinningTool = r) %>% 
        mutate (AssemblyTool = r) %>% 
        mutate(binID=basename(value))
      rgiDF <- rbind(rgiDF, rgi)
    }
}

#save as combined datatable so it's easier for downstream analysis to avoid parsing. 
write.table(rgiDF, "rgi.tsv", sep = "\t", col.names=TRUE, row.names =FALSE)

#transform data into difference views
allHits <- rgiDF %>% group_by(Type, AssemblyTool, BinningTool) %>% count(binID) #get a total count of genes
allHits_excludingNudge <- rgiDF %>% filter(is.na(Nudged)) %>% group_by(Type, AssemblyTool, BinningTool) %>% count(binID) #total count of genes without the nudge flag
allHits_bins <- allHits %>% filter(!stringr::str_detect(binID,'unbinned')) #bins only 
allHits_bins_ref <- allHits_bins %>% filter(Type=="references") #references only
allHits_Unbinned <- allHits %>% filter(stringr::str_detect(binID,'unbinned')) #unbinned only
allHits_excludingNudge_bins <- allHits_excludingNudge %>% filter(!stringr::str_detect(binID,'unbinned')) #binned only, excluding nudge
allHits_excludingNudge_Unbinned <- allHits_excludingNudge %>% filter(stringr::str_detect(binID,'unbinned')) #unbinned only, excluding nudge
allHit_bins_total <- allHits_bins %>% group_by(Type,BinningTool, AssemblyTool) %>% summarize(avg = sum(n)) #summarize countss in binned only
allHits_Unbinned_total <- allHits_Unbinned  %>% group_by(Type,BinningTool, AssemblyTool) %>% summarize(avg = sum(n)) #summarize counts in unbinned only
allHit_excludingNudge_bins_total <- allHits_excludingNudge_bins %>% group_by(Type,BinningTool, AssemblyTool) %>% summarize(avg = sum(n)) #summarize counts in binned only, without nudge
allHits_excludingNudge_Unbinned_total <- allHits_excludingNudge_Unbinned  %>% group_by(Type,BinningTool, AssemblyTool) %>% summarize(avg = sum(n))#summarize counts in unnbined only, without nudge

#plot the above summaries
p <- ggplot(allHits_bins, aes(x=BinningTool, y=n, color=AssemblyTool)) +
  ggtitle("MAG Bins") +
  ylab("Number Of AMR Genes per bin") + 
  geom_boxplot(outlier.colour = 'black', outlier.shape = 16, outlier.size=1, notch = FALSE)
png(paste0("amr_bins_all.png"), width = 1000, height = 500)
plot(p)
dev.off()

p <- ggplot(allHits_excludingNudge_bins, aes(x=BinningTool, y=n, color=AssemblyTool)) +
  ggtitle("MAG Bins, excluding Nudge") +
  ylab("Number Of AMR Genes per bin") + 
  geom_boxplot(outlier.colour = 'black', outlier.shape = 16, outlier.size=1, notch = FALSE)
png(paste0("amr_bins_noNudge.png"), width = 1000, height = 500)
plot(p)
dev.off()

p <- ggplot(allHits_Unbinned, aes(x=BinningTool, y=n, fill=AssemblyTool)) +
  ggtitle("MAG unbinned reads") +
  ylab("Number Of AMR Genes per bin") + 
  geom_boxplot(outlier.colour = 'black', outlier.shape = 16, outlier.size=1, notch = FALSE)
png(paste0("amr_unbins.png"), width = 1000, height = 500)
plot(p)
dev.off()

p <- ggplot(allHits_excludingNudge_Unbinned, aes(x=BinningTool, y=n, color=AssemblyTool)) +
  ggtitle("MAG unbinned reads, excluding Nudge") +
  ylab("Number Of AMR Genes per bin") + 
  geom_boxplot(outlier.colour = 'black', outlier.shape = 16, outlier.size=1, notch = FALSE)
png(paste0("amr_unbins_noNudge.png"), width = 1000, height = 500)
plot(p)
dev.off()

p <- ggplot(allHit_bins_total, aes(x=BinningTool, y=avg, fill=AssemblyTool)) +
  ggtitle("MAG bins") +
  ylab("total number of AMR in bins") + 
  geom_bar(stat = "identity", position = "dodge", color = "grey20")
png(paste0("amr_bins_total.png"), width = 1000, height = 500)
plot(p)
dev.off()

p <- ggplot(allHits_Unbinned_total, aes(x=BinningTool, y=avg, fill=AssemblyTool)) +
  ggtitle("MAG bins") +
  ylab("total number of AMR in bins") + 
  geom_bar(stat = "identity", position = "dodge", color = "grey20")
png(paste0("amr_unbins_total.png"), width = 1000, height = 500)
plot(p)
dev.off()

p <- ggplot(allHit_excludingNudge_bins_total, aes(x=BinningTool, y=avg, fill=AssemblyTool)) +
  ggtitle("MAG bins, no nudge") +
  ylab("total number of AMR in bins") + 
  geom_bar(stat = "identity", position = "dodge", color = "grey20")
png(paste0("amr_bins_noNudge_total.png"), width = 1000, height = 500)
plot(p)
dev.off()

p <- ggplot(allHits_excludingNudge_Unbinned_total, aes(x=BinningTool, y=avg, fill=AssemblyTool)) +
  ggtitle("MAG bins, no nudge") +
  ylab("total number of AMR in bins") + 
  geom_bar(stat = "identity", position = "dodge", color = "grey20")
png(paste0("amr_unbins_noNudge_total.png"), width = 1000, height = 500)
plot(p)
dev.off()


#now copy/paste the exact code process for virulence factors

refType <-c("chromosome","plasmid")
binningTools <- c("concoct","dastool","maxbin2","metabat2")
assemblyTools <- c("idba_ud","megahit","metaspades")

#txt are RGI, vf are VFDB
#lets make 1 GIANT RGI/VFDB file
vfDF = data.frame()
dfColNames =  c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

t = "mag"
for (b in binningTools)
{
  for (a in assemblyTools)
  {
    #lets make 1 GIANT RGI/VFDB file
    for (value in list.files(path=paste0("results/",t,"/",b,"_results/",a,"_results/"), pattern="*.vf", full.names=TRUE, recursive=FALSE) %>% as.vector())
    {
      if (!file.size(value)==0)
      {
        print(value)
        #value="results/mag/concoct_results/idba_ud_results/10.fasta.vf"
        vf <- read_tsv(value, col_names = dfColNames) %>% 
          filter(evalue <= 0.001) %>% 
          #group_by(qseqid, bitscore) %>% top_n(n=1) %>%
          mutate(Type = t) %>% 
          mutate (BinningTool = b) %>% 
          mutate (AssemblyTool = a) %>% 
          mutate(binID=basename(value)) %>% 
          group_by(qseqid,bitscore) %>% top_n(n=1) %>% as.data.frame()

        vfDF <- rbind(vfDF, vf)
      }
    }
  }
}
t="references"
for (r in refType)
{
  #lets make 1 GIANT RGI/VFDB file
  for (value in list.files(path=paste0("results/references/",r,"_results/"), pattern="*.vf", full.names=TRUE, recursive=FALSE) %>% as.vector())
  {
    if (!file.size(value)==0)
    {
      print(value)
      vf <- read_tsv(value, col_names = dfColNames) %>% 
        filter(evalue <= 0.001) %>%  #have to filter by evalue to remove some of the hits.
        #group_by(qseqid, bitscore) %>% top_n(n=1) %>%
        mutate(Type = t) %>% 
        mutate (BinningTool = r) %>% 
        mutate (AssemblyTool = r) %>% 
        mutate(binID=basename(value)) %>% 
        group_by(qseqid,bitscore) %>% top_n(n=1) %>% as.data.frame()
      
        vfDF <- rbind(vfDF, vf)
    }
  }
}

write.table(vfDF, "vf_all0.001.tsv", sep = "\t", col.names=TRUE, row.names =FALSE)
View(allHits_bins)
allHits <- vfDF  %>% group_by(Type, AssemblyTool, BinningTool) %>% count(binID) 

allHits_bins <- allHits %>% filter(!stringr::str_detect(binID,'unbinned'))
allHits_Unbinned <- allHits %>% filter(stringr::str_detect(binID,'unbinned')) 

allHit_bins_total <- allHits_bins %>% group_by(Type,BinningTool, AssemblyTool) %>% summarize(avg = sum(n))
allHits_Unbinned_total <- allHits_Unbinned  %>% group_by(Type,BinningTool, AssemblyTool) %>% summarize(avg = sum(n)) 

p <- ggplot(allHits_bins, aes(x=BinningTool, y=n, color=AssemblyTool)) +
  ggtitle("MAG Bins") +
  ylab("Number Of VF Genes per bin") + 
  geom_boxplot(outlier.colour = 'black', outlier.shape = 16, outlier.size=1, notch = FALSE)
png(paste0("vf_bins_all.png"), width = 1000, height = 500)
plot(p)
dev.off()

p <- ggplot(allHits_Unbinned, aes(x=BinningTool, y=n, fill=AssemblyTool)) +
  ggtitle("MAG unbinned reads") +
  ylab("Number Of VF Genes per bin") + 
  geom_boxplot(outlier.colour = 'black', outlier.shape = 16, outlier.size=1, notch = FALSE)
png(paste0("vf_unbins.png"), width = 1000, height = 500)
plot(p)
dev.off()

p <- ggplot(allHit_bins_total, aes(x=BinningTool, y=avg, fill=AssemblyTool)) +
  ggtitle("MAG bins") +
  ylab("total number of VF in bins") + 
  geom_bar(stat = "identity", position = "dodge", color = "grey20")
png(paste0("vf_bins_total.png"), width = 1000, height = 500)
plot(p)
dev.off()

p <- ggplot(allHits_Unbinned_total, aes(x=BinningTool, y=avg, fill=AssemblyTool)) +
  ggtitle("MAG bins") +
  ylab("total number of VF in bins") + 
  geom_bar(stat = "identity", position = "dodge", color = "grey20")
png(paste0("vf_unbins_total.png"), width = 1000, height = 500)
plot(p)
dev.off()

