#Process read coverages and depths of each accssio and generate figure.
#input: 
#1. Accession - BioSample ID map
#2. List of Chromosomal accessions
#3. Coverage and depth generated using `Samtools coverage` and `samtools depth`
#output:
#1. Average coverage across chromosome, plasmid and GI of all accessions
#2. Depth and coverage of each species(biosample id)
#3. Per base depth and coverage of each species (biosample id)

library(tidyverse)
library(ggpubr)
library(png)
options(scipen=10000)


#global parameters
basePath = "H:\\OneDrive\\School_Work\\_MAGAssessment\\newest\\alignment\\genomecov\\"
chrPrefix="chromosome\\"
plasmidPrefix="plasmid\\"
depthSuffix=".fasta.sam.sorted.bam.depth"
coverageSuffix=".fasta.sam.sorted.bam.coverage"
taxaMapPath="H:/OneDrive/School_Work/_MAGAssessment/newest/MAGAssessment/taxa_metadata.tsv"
colors <- c("#91bfdb", "#abdda4", "#fc8d59")
dir.create("coveragePlots\\summary", showWarnings = FALSE)
dir.create("coveragePlots\\perbase", showWarnings = FALSE)

#chromosome accession list
chromosomeAccs <- scan(paste0(basePath,"chromosome.txt"), what="", sep="\n")
#mapping files
taxaMetadata <- read_tsv(taxaMapPath, col_names=c("spp","bioAccession", "genomeAccession", "type"))
chromosomeRefBioAccMap <- taxaMetadata %>% filter(type=="chromosome") %>% select(genomeAccession, bioAccession, spp)
plasmideRefBioAccMap <- taxaMetadata %>% filter(type=="plasmid") %>% select(genomeAccession, bioAccession, spp)


#function to append a few extra columns (region, binary coverage) into the depth df. 
FillMissingBases <- function(depthDF, id, start, end, gis, cov, type=1){
  filledDepthDF <- depthDF

  if (nrow(gis) > 0){
    for (i in 1:nrow(gis)){
      #theree is a start and a stop.
      coverageLabelDF <- filledDepthDF %>% mutate (region = ifelse(X2 %in% c(gis$Start[i]:gis$End[i]), "GI", "Chromosome"))
    }
  } else if (type == 2){
    coverageLabelDF <- filledDepthDF %>% mutate (region = "Plasmid")
  }
  else {
    coverageLabelDF <- filledDepthDF %>% mutate (region = "Chromosome")
  }
  
  coverageLabelDF <- coverageLabelDF %>% mutate(coverage = cov) %>% mutate (covered = ifelse(X3 == 0, 0, 1))
  
  colnames(coverageLabelDF) <- c("acc", "base", "depth", "region", "coverage", "covered" )
  
  coverageLabelDF <- coverageLabelDF %>% mutate(acc = id)
  
  return(coverageLabelDF)
}

#function to take 1 spp's depth df and generate per base depth and average depth across region plots.
PlotCoverageByBase <- function(depthDF, gis, windowSize, acc, coverage, yAxisLimit=NULL, type=1, color="#91bfdb"){
  #sliding window average coverage
  n <- windowSize;
  avgDF<- aggregate(depthDF,list(rep(1:(nrow(depthDF)%/%n+1),each=n,len=nrow(depthDF))),mean)[-1]  %>% select(base, depth)

  #plot the coverage
  print("Plotting")
  
  if (is.null(yAxisLimit)){
    yAxisLimit = max(avgDF$depth)
    }
  
  p<- ggplot(data=avgDF, aes(x=base, y=depth)) +
    geom_line(size = 0.75) +
    # geom_point(size = 0.75) +
    xlim(0, max(avgDF$base)) +
    #ylim(0, yAxisLimit) +
    ylab("Depth") +
    ggtitle(paste0(acc, "\nCoverage: ", coverage, "%" )) +
    theme(plot.title = element_text(hjust = 0.5, size = 9)) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(panel.background = element_rect(fill = "transparent", colour = color, size = 2)) + 
    scale_y_log10(limits = c(1, 10000))

  
  if (type != 2){
    p<- p+ xlab("Genomic position")
  } else {
    p<- p+ xlab("Plasmid position")
    
  }
  
  if (nrow(gis)>1 && type != 2){
    p <- p + geom_rect(data=gis, inherit.aes=FALSE, aes(xmin=Start, xmax=End, ymin=1,
                                               ymax=10000), color="transparent", fill="#abdda4", alpha=0.3) 
  }

  #ggsave(paste0("coveragePlots\\coverage_", acc, ".png"), plot = p, device = "png", width = 8, height = 4)
  
  return (p)
}


#function for main logic for taking in a chromosomal accession and generating figures
MakeCoveragePlots <- function(accssion){
  #accssion <- "NC_004337"
  depth.chr.path <- paste0(basePath,chrPrefix,accssion,depthSuffix)
  coverage.chr.path <- paste0(basePath,chrPrefix,accssion,coverageSuffix)
  #find the plasmids associated with this chr accession. 
  bioAcc <- chromosomeRefBioAccMap %>% filter(genomeAccession==accssion) %>% select(bioAccession) %>% toString()
  spp <- chromosomeRefBioAccMap %>% filter(genomeAccession==accssion) %>% select(spp) %>% toString()
  
  plasmidAccessions <- plasmideRefBioAccMap %>% filter(bioAccession==bioAcc) %>% select(genomeAccession) %>% c()
  
  print("Reading files")
  depthDF.chr <- read_tsv(depth.chr.path, col_names = FALSE)#"H:\\OneDrive\\School_Work\\_MAGAssessment\\newest\\alignment\\genomecov\\chromosome\\NC_004337.fasta.sam.sorted.bam.depth", col_names = F)
  coverageDF.chr <- read_tsv(coverage.chr.path, col_names = TRUE)#"H:\\OneDrive\\School_Work\\_MAGAssessment\\newest\\alignment\\genomecov\\chromosome\\NC_004337.fasta.sam.sorted.bam.coverage",  col_names = T)

  name = (coverageDF.chr$`#rname` %>% strsplit(".", fixed = T))[[1]][1]
  start <- coverageDF.chr$startpos
  end <- coverageDF.chr$endpos
  coverage <- coverageDF.chr$coverage

  #find GI regions
  gis <- read_tsv("H:\\OneDrive\\School_Work\\_MAGAssessment\\data\\all_gis_islandpath_dimob_iv4.txt", col_names = T) %>% 
    filter(Accession_number==name) %>% select(-Accession_number, -Prediction_method)

  print("combining")
  #populate a df with depth data and labels
  combinedDepthDF <- FillMissingBases(depthDF.chr, name, start, end, gis, coverage)
  #View(combinedDepthDF)
  
  print("combining plasmid")
  maxDepth <- 1000
  for (plasmid in plasmidAccessions[[1]]){
 #   plasmid<-"NC_004632"
    print(plasmid)
    depth.plasmid.path <- paste0(basePath,plasmidPrefix,plasmid,depthSuffix)
    coverage.plasmid.path <- paste0(basePath,plasmidPrefix,plasmid,coverageSuffix)
    
    depthDF.plasmid<- read_tsv(depth.plasmid.path, col_names = FALSE)#"H:\\OneDrive\\School_Work\\_MAGAssessment\\newest\\alignment\\genomecov\\chromosome\\NC_004337.fasta.sam.sorted.bam.depth", col_names = F)
    coverageDF.plasmid <- read_tsv(coverage.plasmid.path, col_names = TRUE)#
    
    name = (coverageDF.plasmid$`#rname` %>% strsplit(".", fixed = T))[[1]][1]
    start <- coverageDF.plasmid$startpos
    end <- coverageDF.plasmid$endpos
    coverage <- coverageDF.plasmid$coverage
    plasmidDF<- FillMissingBases(depthDF.plasmid, name, start, end, data.frame(), coverage, 2)
    if (max(plasmidDF$depth) >= maxDepth){
      maxDepth = max(plasmidDF$depth)
    }
    combinedDepthDF <- rbind(combinedDepthDF, plasmidDF)
  }
  
  p.summary <- ggplot(combinedDepthDF, aes(x=region, y=depth, color=region)) + 
    geom_boxplot(outlier.colour = 'black', outlier.shape = NA, outlier.size=0.25, notch = FALSE) +
    theme_light() +
    xlab("Region") +
    ylab("Depth") +
    scale_color_manual(values=colors)+
    stat_summary(fun=mean, geom="point", shape=23, size=1) +
    ggtitle(paste0(bioAcc, "\n", spp)) +
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    #theme(legend.title = element_blank()) + 
    scale_y_log10(limits = c(1, 10000))

  ggsave(paste0("coveragePlots\\summary\\summary_coverage_",bioAcc, ".png"), plot = p.summary, device = "png", width = 5, height = 4)
  
  maxDepth =10000
  
  
  print("plotting chr")
  combinedDepthDF.chr <- combinedDepthDF %>% filter(region!="Plasmid")
  plotList <- list()
  p <- PlotCoverageByBase(combinedDepthDF.chr, gis, 10000, combinedDepthDF.chr$acc[1], combinedDepthDF.chr$coverage[1], maxDepth)
  plotCount=1
  plotList[[plotCount]] <- p

  print("plotting plasmid")
  for (plasmid in plasmidAccessions[[1]]){
   # plasmid<-"NC_004632"
    plotCount <- plotCount + 1
    combinedDepthDF.plas <- combinedDepthDF %>% filter(acc == plasmid)
    p1 <- PlotCoverageByBase(combinedDepthDF.plas, gis, 100, combinedDepthDF.plas$acc[1], combinedDepthDF.plas$coverage[1], maxDepth, 2, "#fc8d59")
    plotList[[plotCount]] <- p1
  }
  
  width <- c(3)
  for (i in 2:plotCount){
    width <- c(width, 1)
  }
  
  print("plotting all")
  figure <- ggarrange(plotlist = plotList, ncol=(length(plasmidAccessions[[1]]) + 1), nrow=1, widths = width) %>%
    annotate_figure(top = text_grob(paste0(bioAcc," - ", spp)))
  
  ggsave(paste0("coveragePlots\\perbase\\per_base_coverage_",bioAcc, ".png"), plot = figure, device = "png", width = 14, height = 4)
  
  
  df <- data.frame()
 
  combinedDepthDF.chr <- combinedDepthDF %>% filter(region=="Chromosome" )
  combinedDepthDF.gi <- combinedDepthDF %>% filter(region=="GI")
  combinedDepthDF.plas <- combinedDepthDF %>% filter(region=="Plasmid")
  coverage.chr <- ((combinedDepthDF.chr  %>% filter(covered == 1)%>% count()) / nrow(combinedDepthDF.chr)) * 100
  coverage.gi <- ((combinedDepthDF.gi  %>% filter(covered == 1)%>% count()) / nrow(combinedDepthDF.gi)) * 100
  coverage.plas <- ((combinedDepthDF.plas  %>% filter(covered == 1) %>% count()) / nrow(combinedDepthDF.plas)) * 100

  
  df <- rbind(df, data.frame("Region" = "Chromosome", "Coverage" = coverage.chr))
  if (nrow(combinedDepthDF.gi) > 0)
    df <- rbind(df, data.frame("Region" = "GI", "Coverage" = coverage.gi))
  if (nrow(combinedDepthDF.plas) > 0)
    df <- rbind(df, data.frame("Region" = "Plasmid", "Coverage" = coverage.plas))
  
  
  return(df)
}

#Region Main()
coveragesAll<-data.frame() #df to store coverage data across all spp

for (acc in (chromosomeAccs)){
  
  acc =  (acc %>% strsplit(".", fixed = T))[[1]][1]
  #acc <- "NC_004337"
  print(acc)
  temp <- MakeCoveragePlots(acc)
  coveragesAll <- rbind(coveragesAll, temp)
}

colnames(coveragesAll) <- c("Region", "Coverage")

#plot the average coverage figure
p.coverage <- ggplot(coveragesAll, aes(x=Region, y=Coverage, color=Region)) + 
  geom_boxplot(outlier.colour = 'black', outlier.shape = 16, outlier.size=1, notch = FALSE) +
  theme_light() +
  ylim(75, 100) +
  scale_color_manual(values=colors)+
  stat_summary(fun=mean, geom="point", shape=23, size=1) +
  ggtitle("Average Coverage By Genome Region") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) 

ggsave(paste0("coverage.png"), plot = p.coverage, device = "png", width = 5, height = 4)

#condense the individual summary and perbase plots into 1 figure.
pngList.summary <- list.files(path = "coveragePlots\\summary", pattern = "*.png", full.names = TRUE)
pngList.perbase <- list.files(path = "coveragePlots\\perbase", pattern = "*.png", full.names = TRUE)

#function to combine all the indiviudal perbase/average depth into 1 figure
PlotFinalFigures <- function(FilesList, title, fileName, numPerPage=6, numCol=2, numRow=3){
  count = 1
  figures <- list()
  for (files in FilesList){
    print(files)
    img <- readPNG(files)
    plot <- ggplot() + 
      background_image(img) + 
      theme(plot.margin = margin(t=0.5, l=0.5, r=0.5, b=0.5, unit = "cm"))
    figures[[count]] <- plot
    count <- count + 1
  }
  
  splitLists<- split(figures, rep(1:ceiling((length(figures)/numPerPage)), each=numPerPage)[1:(length(figures))])
  count <- 1
  for (list in splitLists) {
    fig <- ggarrange(plotlist = list, ncol=numCol, nrow = numRow) %>%
      annotate_figure(top = text_grob(paste0(title, " (", LETTERS[count], ")"), face = "bold", size = 25))

    ggsave(paste0(fileName, "_", LETTERS[count], ".png"), plot = fig, device = "png", width = 8.5, height = 11.75)
    count <- count + 1
  }
}

summaryFigList <- PlotFinalFigures(pngList.summary, "Depth of Simulated Reads By Species", "DepthByspp", 15, 3, 5)
perbaseFigList <- PlotFinalFigures(pngList.perbase, "Per Base Depth of Simulated Reads By Species", "PerBaseDepthBySpp", 16, 2, 8)

