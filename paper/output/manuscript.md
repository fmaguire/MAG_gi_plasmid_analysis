---
author-meta:
- Finlay Maguire*
- Baofeng Jia*
- Kristen Gray
- Venus Lau
- Robert G. Beiko
- Fiona S.L. Brinkman
date-meta: '2019-10-28'
keywords:
- markdown
- publishing
- manubot
lang: en-CA
title: Metagenome-Assembled Genome Binning Methods Disproportionately Fail for Plasmids
  and Genomic Islands
...






<small><em>
This manuscript
was automatically generated
on October 28, 2019.
</em></small>

## Authors
Please note the current author order is chronological and does not reflect the final order.



+ **Finlay Maguire***<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [0000-0002-1203-9514](https://orcid.org/0000-0002-1203-9514)
    · ![GitHub icon](images/github.svg){.inline_icon}
    [fmaguire](https://github.com/fmaguire)
    · ![Twitter icon](images/twitter.svg){.inline_icon}
    [fmaguire](https://twitter.com/fmaguire)<br>
  <small>
     Faculty of Computer Science, Dalhousie University
     · Funded by ['Genome Canada', 'Donald Hill Family Fellowship']
  </small>

+ **Baofeng Jia***<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [XXXX-XXXX-XXXX-XXXX](https://orcid.org/XXXX-XXXX-XXXX-XXXX)
    · ![GitHub icon](images/github.svg){.inline_icon}
    [imasianxd](https://github.com/imasianxd)<br>
  <small>
     Department of Biochemistry and Molecular Biology, Simon Fraser University
  </small>

+ **Kristen Gray**<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [XXXX-XXXX-XXXX-XXXX](https://orcid.org/XXXX-XXXX-XXXX-XXXX)<br>
  <small>
     Department of Biochemistry and Molecular Biology, Simon Fraser University
  </small>

+ **Venus Lau**<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [XXXX-XXXX-XXXX-XXXX](https://orcid.org/XXXX-XXXX-XXXX-XXXX)<br>
  <small>
     Department of Biochemistry and Molecular Biology, Simon Fraser University
  </small>

+ **Robert G. Beiko**<br><br>
  <small>
     Faculty of Computer Science, Dalhousie University
  </small>

+ **Fiona S.L. Brinkman**<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [XXXX-XXXX-XXXX-XXXX](https://orcid.org/XXXX-XXXX-XXXX-XXXX)<br>
  <small>
     Department of Biochemistry and Molecular Biology, Simon Fraser University
  </small>




## Abstract {#abstract}

**final numbers to add**

Metagenomic methods, in which all the DNA in sample is simultanously sequenced, is an increasingly popular method in the life sciences.
They have a major advantage over genomic or phenotypic methods as they do not require time-intensive and bias-inducing culturing steps.
This means a much greater diversity can be profiled with minimal _a priori_ assumptions.
Due to this strength, metagenomics is emerging as a key tool in public health microbiology for surveillance of virulence and antimicrobial resistance (AMR) genes.
The most important sequences for surveillance purposes are those associated with mobile genetic elements such as plasmids and genomic islands (GIs).
Unfortunately, metagenomic data, even when assembled, results in complex mixed set of DNA fragments rather than nicely resolved individual genomes.
Recently, methods have been developed that attempt to group these fragments into bins likely to have been derived from the same underlying genome.
These bins are commonly known as metagenome-assembled genomes (MAGs). 
MAG based approaches have been used to great effect in revealing huge amounts of previously uncharacterised microbial diversity.
These methods perform this grouping using aspects of the sequence composition and the relative abundance of that sequence in the dataset.
Unfortunately, plasmids are often represented at different copy numbers than the corresponding chromosomes. 
Additionally, both plasmids and genomic islands often feature significantly different sequence composition than the rest of the source genome as a whole.
Due to this we hypothesise these types of sequences will be highly under-represented in MAG based approaches.

To evaluate this we generated a simulated metagenomic dataset comprised of genomes with large numbers of plasmids and considerable proportion of chrosomomal DNA consisting of GIs at varying relative abundances.
MAGs were then recovered from this data using a variety of different established MAG pipelines and parameterisations and correct binning of plasmid and GI sequences calculated relative to the genomes as a whole.
We show that regardless of the MAG approach used, plasmid and GI dominated sequences will systematically be left unbinned or incorrectly binned.
This indicates the importance of read based approaches for thorough evalaution of resistome complements in metagenomic data.


## Introduction {#intro}

Metagenomics, the untargeted sequencing of all DNA within a sample, has become the dominant approach for characterising viral and microbial communities over the last 17 years [@8PLOeAH6; @7RV1Ygsv]. By targeting all genomic contents, these methods allow researchers to profile the functional potential and the taxonomic composition of a sample. This is in contrast to barcoding based approaches such as 16S or 18S rRNA sequencing which only provide taxonomic information [@Trbw5ZnC] (although you can attempt to predict functional potential from taxonomic data [@93ojQnSg; @L09Dmaq3]). One of many areas where metagenomics has been very useful is in the analysis of antimicrobial resistance (AMR). Using these approaches has been instrumental in developing our understanding of the distribution and evolutionary history of AMR genes [@rwhLEYRY; @QSe5BqFk; @2xaXclNM]. It has also formed a very useful tool for pathogen tracking in public health outbreak analyses [@khJQfjDf]

While 3rd generation long-read technology has begun to be adopted in metagenomics analyses [@U4vhNZoB; @1759XyDVi] the majority of analyses still involve high-throughput 2nd generation sequencing. These 2nd generation platforms such as Illumina's MiSeq provide high numbers (10s-100s of millions) of relatively short reads (150-250bp) randomly sampled from the underlying DNA in the sample. This sampling is therefore in proportion to the relative abundance of different organisms (i.e. more abundant organisms will be more represented in the reads). There are 2 main approaches to the analysis of 2nd generation metagenomic data: read homology and metagenome assembly. Read-based approaches involve using reference databases and BLAST based sequence similarity search tools (e.g. DIAMOND [@4R96QRcV]), read mapping (e.g. Bowtie 2 [@PiS0h6Mu]), Hidden Markov Models (e.g. HMMER3 [@77xWEk9S]) or k-mer hashing (e.g. CLARK [@OoKZ0WcH]). These read-based approaches allow analysis of all reads with detectable similarity to the genes you are interested even if the organism is relatively under-represented in the dataset. However, read-based methods are reliant on quality of the reference database (i.e. you don't detect things you don't already know about) and does not provide any information about the genomic organisation of the genes.

In order to get more data about the relative genomic context and organisation of your genes of interest it is possible (although computationally demanding) to assemble the short reads into longer fragments of DNA (contigs). This approach has been used successfully in even very early metagenomic analysis papers [@F7RexqdF]. There are a variety of specialised de Bruijn graph assemblers developed to handle the particular challenges of this type of assembly (such as metaSPAdes [@KP5SjPXN] , IDBA-UD [@a4mT7fuU], and megahit [@1EUV0Ejkr]) each with a range of different strengths and weaknesses [@f0X2CPKM]. While metagenomic assembly does provide longer stretches of DNA incorporating information about multiple genes without further analysis it still leaves you with a large collection of DNA fragments with no obvious groupings.

An increasingly common way to deal with this is to attempt to group these assembled contigs into bins all derived from the same underlying genome in the sample. These resulting bins are known as metagenome assembled genomes (MAGs). This binning is typically performed by grouping all the contigs with similar abundance and similar sequence composition into the same bin. A range of tools have been released to perform this binning including CONCOCT [@WRoCf6pg], MetaBAT 2 [@b2WO18xh], and MaxBin 2 [@sG4CX8Uj]. There is also the metabinning tool DAS Tool [@DfIRBmdF] which combines predictions from multiple binning tools together. These MAG approaches have been used to great effect in unveiling huge amounts of previously uncharacterised genomic diversity [@4rsFboY4; @wrBRBdFb; @Rk2NATlI].

Unfortunately, only a relatively small proportion of reads are successfully assembled and binned in large complex metagenome datasets e.g. 24.2-36.4% of reads from permafrost [@buqrbdBh] and soil metagenomes [@d5Hh0941]. Additionally, a large number of detected genomes are not reconstructed at all with only ~23% of all detected genomes recovered in the soil metagenome [@d5Hh0941]. There have been attempts to benchmark and compare these tools such as the Critical Assessment of Metagenome Interpretation (CAMI) challenge's (https://data.cami-challenge.org/) Assessment of Metagenome BinnERs (AMBER) [@Y8sHlHi] however these only investigate the completeness and purity of recovered MAGs relative to true genomes in the sample. They don't attempt to assess whether there are specific components of the underlying genomes that are disproportionately lost. Two such genomic elements of great importance to the study and lateral gene transfer of AMR are genomic islands and plasmid sequences.

Genomic islands (GIs) are clusters of genes known or predicted to have been acquired through lateral gene transfer (LGT) events. These include integrons, transposons, integrative and conjugative elements (ICEs) and prophages (integrated phages) [@DET3tBYj; @1Af4oXwEX]. They have been shown to disproportionately encode virulence factors [@LxGqo7iq] and are a major mechanism of LGT of AMR genes [@x7HhCKyS; @17U91060Y]. However, these GIs often have different nucleotide composition compared to the rest of the genome [@DET3tBYj], which is exploited by tools such as SIGI-HMM [@19UeQywMr] and IslandPath-DIMOB [@M1pdcdMy] to detect GIs. Additionally, GIs may exist as multiple copies within a genome [@5g9Xc4ot] leading to potential difficulties in correctly assembling these regions in metagenome assemblies as well as likely biases in the calculation of coverage statistics. Similarly, plasmids are a major source of the dissemination and transistion of AMR genes throughout microbial ecosystems [@X9j9vETu; @x7HhCKyS]. They also exist at variable copy number [@qtpTcNWp; @Z1irb7eF] and with markedly different sequence composition to the genome they are associated with [@QK9dmRUA; @ps1aOiRU]

As MAG binning is performed on the basis of sequence composition and coverage this suggests that these sequences are liable to being incorrectly binned or lost in the process of recovering MAGs. Due to the importance of these genomic components in the function and spread of pathogenic traits such as AMR and virulence it is vital that we assess the impact of assembly and binning on the representation of these elements. This is particularly important with the increasing popularity of MAG approaches within microbial and public health research. Therefore, to address this issue we performed an analysis of GI and plasmid recovery accuracy across a range of assembly and binning approaches using a simulated medium complexity metagenome comprised of GI and plasmid rich taxa.


## Materials and Methods {#methods}

All analyses presented in this paper can be reproduced and inspected using the Jupyter (version XX) [citation] notebook (.ipynb) within the associated github repository (https://github.com/fmaguire/MAG\_gi\_plasmid\_analysis).
The specific code version used for this paper is also archived under DOI:XXXYYY.

### Metagenome Simulation

All genomes were selected from the set of completed RefSeq genomes as of April 2019.  
Genomic islands for these genomes were previously predicted using IslandPath-DIMOB [@M1pdcdMy] and collated into the IslandViewer database (http://www.pathogenomics.sfu.ca/islandviewer) [@4eEyIkDg].
Plasmid sequences and numbers were recovered for each genome using the linked genbank Project IDs.

30 genomes were arbitrarily chosen to exemplify the following criteria:
    - 10 with high numbers of plasmids
    - 10 with a very high proportion (\>10%) of chromosomes corresponding to composition detected GIs.
    - 10 with a very low proportion (\<1%) of chromosomes corresponding to composition detected GIs.

The data used to select the taxa is listed in Supplemental Table 1 and the details of the selected subset taxa are listed in Supplemental Table 2 with their NCBI accessions.
The sequences themselves are in `data/sequences/sequences.tar.bz2`

In accordance to the recommendation in the CAMI challenge [@lsbnKJf8] the genomes were randomly assigned a relative abundance following a log-normal distribituion (\mu = 1, \sigma = 2).
Plasmid copy number estimates could not be accurately found for all organisms, therefore, plasmids were randomly assigned a copy number regime: low (1-20), medium (20-100), or high (500-1000) at a 2:1:1 rate.
Within each regime the exact copy number was selected using an appropriately scaled gamma distribution (\alpha = 4, \beta = 1) or the minimum edge of the regime.
Finally, the effective plasmid relative abundance was determined by multiplying the plasmid copy number with the genome relative abundance.
The full set of randomly assigned relative abundances and copy numbers can be found in Supplemental Table 3.

Sequences were then concatenated into a single fasta file with the appropriate relative abudance. 
MiSeq v3 250bp paired-end reads with a mean fragment length of 1000bp (standard deviation of 50bp) were then simulated using art\_illumina (v2016.06.05) [@znONJtTo] at a fold coverage of 2.9 resulting in a simulate metagenome of 31,174,411 read pairs.

The selection of relative abundance and metagenome simulation itself was performed using the `data_simluation/simulate_metagenome.py` script.

### Metagenome Assembled Genome Recovery

Reads were trimmed using sickle (v1.33) [@1CBlSILo4] resulting in 25,682,644 surviving read pairs.
The trimmed reads were then assembled using 3 different metagenomic assemblers: metaSPAdes (v3.13.0)[@KP5SjPXN] , IDBA-UD (v1.1.3) [@a4mT7fuU], and megahit (v1.1.3) [@1EUV0Ejkr]).
The resulting assemblies were summarised using metaQUAST (v5.0.2) [@TeRvtMCl].
The assemblies were indexed and reads mapped back using Bowtie 2 (v2.3.4.3) [@PiS0h6Mu].

Samtools (v1.9) were then used to sort the read mappings and the read coverage calculated using the MetaBAT 2 accessory script (jgi\_summarize\_bam\_contig\_depths. 
The 3 metagenome assemblies were then separately binned using CONCOCT (v0.4.2) [@WRoCf6pg], MetaBAT 2 (v2.13) [@b2WO18xh], and MaxBin 2 (v2.2.6) [@sG4CX8Uj]. 
As per the specified manual instructions CONCOCT used a slightly different approach to estimate read coverage.
The supplied accessory scripts were also used to cut contigs into 10 kilobase fragments (cut\_up\_fasta.py) and read coverage calculated for the fragments (concoct\_coverage\_table.py).
These fragment coverages were then used to bin the 10kb fragments before the clustered fragments were merged (merge\_cutup\_clustering.py) to create the final CONCOCT MAG bins (extra\_fasta\_bins.py)
Finally, for each metagenome assembly the predicted bins from these 3 binners were combined using DAS Tool (v1.1.1) [@DfIRBmdF].
This resulted in 12 separate sets of MAGs (one set for each assembler and binner pair).

### MAG assessment

#### Chromosomal Coverage

The MAG assessment for chromosomal coverage was performed by creating a BLASTN 2.9.0+ [@nEsJGUWa] database 
consisting of all the chromosomes of the input reference genomes.
Each MAG contig was then used as a query against this database and the coverage of the underlying chromosomes tallied by merging the overlapping aligning regions and summing the total length of aligned MAG contigs.
The most represented genome in each MAG was assigned as the "identity" of that MAG for further analyses.
Coverages less than 5% were filtered out and the number of different genomes that a MAG contain contigs aligning to were tallied.
Finally, the overall proportion of chromosomes that were not present in any MAG were tallied for each binner and assembler.

#### Plasmid and GI Coverage

Plasmid and GI coverage were assessed in the same way as one another.
Firstly, each a BLASTN database was generated for each set of MAG contigs.
Then each MAG database was searched for plasmid and GI sequences.
Any plasmid or GI with greater than 50% coverage in a MAG was retained.
All plasmids or GIs which could be found in the unbinned contigs or the MAGs was recorded as having been successfully assembled.
The subset of these which were found in the binned MAGs was then seperately tallied.
Finally, we evaluated the proportion of plasmids or GIs which were binned correctly in the bin which was maximally composed of chromosomes from the same source gneome.
This was determined using the bin "IDs" from the chromosomal coverage analysis.

### Antimicrobial Resistance and Virulence Factors Assessment
#### Detection of AMR/VF Genes
For each of the 12 MAGs, and the reference chromosome and plasmids, AMR genes were predicted using Resistance Gene Identifier (RGI v5.0.0; default parameters) and the Comprehensive Antibiotic Resistance Database (CARD v3.0.2) [@1ByMfX8Y1]. Virulence factors were predicted using BLASTX against the Virulence Factors Database (VFDB; obtained on Aug 26, 2019) with an e-value cut-off of 0.001. [@pYB1SP5]. Each MAG was then assigned to a reference chromosome and plasmid using the above mentioned mapping criteria for downstream analysis.

#### AMR/VF Gene Recovery
The ability for MAGs to properly recover AMR and VF genes was assessed using only the megahit-DasTool assembler-binner combination as it was the best performing pair. For each bin, we counted the total number of AMR/VF genes recovered then compared this to the number predicted in their assigned reference chromosome and plasmids to determine MAG’s gene recovery ability. We then mapped the location of reference replicon’s predicted genes to the bins to determined the location of those genes in MAGs. 

#### Protein subcellular localization predictions
The MAG bins from megahit-DasTool assembler-binner combination was inputted into prodigal [@lX665mdh] to predict open reading frames (ORFs) using the default parameter. The list of predicted proteins is inputted into PSORTb v3.0 with default parameters [@19bHWbO47]. 


## Results {#results}

The overall ability of MAG methods to recapitulate the original chromosomal source genome results varied widely.
We consider the identity of a given MAG bin to be that of the genome that composes the largest proportion of sequence within that bin.
In other words if a bin is identifiably 70% species A and 30% species B we consider that to be a bin of species A.
Ideally, we wish to generate a single bin for each source genome comprised of the entire genome and no contigs from other genomes.
Some genomes are cleanly and accurately binned regardless of the assembler and binning method used (See Fig. (@fig:supspeciescov)).
Specifically, greater than 90% of _Streptomyces parvulus_ (minimum 91.8%) and _Clostridium baratii_ (minimum 96.4%) chromosomes are represented in individual bins across all methods.
However, no other genomes were consistently recovered by all methods for more than a 3rd of the chromosomes.
The 3 _Streptococcus_ genomes were particularly problematic, likely due to their similarity, with the best recovery for each ranging from 1.7% to 47.49%.

![Top Species Coverage](images/s1_species_top_coverage.png){#fig:supspeciescov}

In terms of assembler, megahit resulted in the highest median chromosomal coverage across all binners (81.9%) with metaSPAdes performing worst (76.8%) (see Fig. (@fig:1topcoverage)).
In terms of binning tool, CONCOCT performed very poorly with a median 26% coverage for top hit per bin, followed by maxbin2 (83.1%), and metabat2 (88.5%).
It is perhaps unsuprising that the best performing binner in terms of bin top hit coverage was the metabinner DAS-TOOL that combines predictions from the other 3 binners (94.3% median top hit chromosome coverage per bin, (@fig:1topcoverage)).

![Chromosomal coverages of most prevalent genome in each bin across binners and metagenome assemblies](images/1a_top_coverage.png){#fig:1topcoverage}

Bin purity, i.e. the number of genomes present in a bin at >5% coverage, was largely equivalent across assemblers (see Fig. (@fig:1purity)), with a very marginally higher purity for IDBA.
In terms of binning tool, however, maxbin2 proved an outlier with nearly twice as many bins containing multiple species as the next binner.
The remaining binning tools were largely equivalent, producing chimeric bins at approximately the same rates.

![Distribution of bin purities. Showing the number of genomes with >5% chromosomal coverage across the bins and methods](images/1b_purity.png){#fig:1purity}

Regardless of method, a very small proportion of plasmids were correctly binned in the bin that mostly contained chromosomal contigs from the same source genome.
Specifically, between 1.5% (IDBA-UD assembly with DAS Tool bins) and 29.2% (metaSPAdes with CONCOCT bins) were correctly binned at over 50% coverage.
In terms of metagenome assembly MetaSPAdes was far and away the most successful assembler at assembling plasmids with 66.2% of plasmids identifiable at greater than 50% coverage.
IDBA-UD performed worst with 17.1% of plasmids recovered, and megahit recovered 36.9%.
If the plasmid was successfully assembled it was placed in a bin by maxbin2 and CONCOCT, although a much smaller fraction correctly binned (typically less than 1/3rd).
Interestingly, metabat2 and DAS tool binners were a lot more conservative in assigning plasmid contigs to bins, however, of those assigned to bins nearly all were correctly binned (see Fig. (@fig:2plasmids)).

![Plasmid Coverage](images/2_plasmid_coverage.png){#fig:2plasmids}

GIs displayed a similar pattern of assembly and correct binning performance as plasmids (see Fig. (@fig:3gis)).
These sequences were assembled uniformly badly (37.8-44.1%) with metaSPAdes outperforming the other two assembly approaches.
For CONCOCT and maxbin2 binning tools all GIs that were assembled were assigned to a bin although the proportion of binned GIs that were correctly binned was lower than for DAS Tool and metabat2.
DAS Tool, metabat2 and CONCOCT didn't display the same precipitious drop-off between those assembled and those correctly binned as was observed for plasmids.
In terms of overall correct binning with the chromosomes from the same genome the metaSPAdes assembly with CONCOCT (44.1%) and maxbin2 (43.3%) binners performed best.

![Genomic Island Coverage](images/3_gi_coverage.png){#fig:3gis}

With respect to AMR genes, in total, MAGs were only able to recover between 40-53% of the AMR genes predicted in our reference genomes across all assembler-binner pairs (Fig. (@fig:4MAGBinTotal)). We then took the best assembler-binner pair (MegaHit-DasTools) and examined the AMR genes recovered in detail. We noticed that, for majority of the bins (85%), MAGs were able to correctly recover either 100% or 0% of the AMR genes (Median value 100%) that are contained in the reference chromosome assigned to that bin. However, MAGs were not able to correctly recover any of the AMR genes that were present on plasmids (Fig. (@fig:5MAGDMBinTotal), Fig. (@fig:6AMRGenePercentRecovery)). Lastly, we asked the question of where reference replicon AMR genes went in the MAGs. For chromosome, majority (81%) of the AMR genes was found in a bin of the MAG. A small portion (12%) was left unbinned and 7% were not found in MAGs at all. On the other hand, for plasmid born AMR genes, all of the recovered genes (n=20) were identified in the unbinned fraction of our MAG (Fig. (@fig:7LocationOfReferenceGenomeAMR)).  

![Total AMR Genes Detected Across Tools](images/4MAGBinTotal.png){#fig:4MAGBinTotal}
![Total AMR Genes Detected In Best Pair](images/5MAGBinTotal.png){#fig:5MAGDMBinTotal}
![AMR Gene Percent Recovery](images/6AMRGenePercentRecovery.png){#fig:6AMRGenePercentRecovery}
![Location of Reference AMR Genes](images/7LocationOfReferenceGenomeAMR.png){#fig:7LocationOfReferenceGenomeAMR}

Aside from AMR genes, we also examined virulence factors in our dataset. The number of VF varied dramatically across assembler-binner combinations. Maxbin2 and metabat2 binned MAGs, regardless of assembly method, predicted over 20,000 virulence genes, roughly 50x more than our reference genome. The DasTool-megahit MAG predicted just under 250 VFs, which represent 61% of the number of VFs in our reference replicons (Fig. (@fig:8vfBinTotal)). Furthermore, each bin’s percent recovery of predicted VF compared to the number in their assigned reference chromosome varied much more compared to AMR genes. The %recovered value ranges from 0% to 112% (Median 0). (Fig. (@fig:9vfDMBinTotal), Fig. (@fig:10VFGenePercentRecovery)). Finally, we examined the location of reference chromosomal VF genes in our MAG bins. 63% of VF genes were found in a MAG bin, 19% were found in the unbinned portion and 18% were not found al all. There were no VF predicted on reference plasmids (Fig. (@fig:11LocationOfReferenceGenomeVF)).  

![Total VF Genes Detected Across Tools](images/8vfBinTotal.png){#fig:8vfBinTotal}
![Total VF Genes Detected In Best Pair](images/9vfDMBinTotal.png){#fig:9vfDMBinTotal}
![VF Gene Percent Recovery](images/10VFGenePercentRecovery.png){#fig:10VFGenePercentRecovery}
![Location of Reference VF Genes](images/11LocationOfReferenceGenomeVF.png){#fig:11LocationOfReferenceGenomeVF}

Lastly, we looked at the ability for MAGs to predict subcellular localization of proteins using PSORTb. Overall, the localization distribution of predicted proteins were very similar in MAGs compared to the reference genome (Fig. (@fig:12subcellularLocalization)).

![Distribution of Predicted Protein Subcellular Localization](images/12subcellularLocalization.png){#fig:12subcellularLocalization}




## Discussion {#discussion}

In this paper, we evaluated the ability and accuracy of metagenome-assembled genomes (MAGs) to correctly recover mobile genetic elements (i.e. genomic islands and plasmids) from metagenomic samples across different tools used to assemble and bin MAGs. 

Overall, the best assembler-binner pair was megahit-DASTOOL in term of both chromosomal coverage (94.3%) and bin purity (1). Looking at genomes with the lowest coverage, the 3 Streptococcus genomes were particularly problematic, likely due to their similarity, with the best recovery for each ranging from 1.7% to 47.49%. This suggest that MAGs might not be able to distinguish between closely related species (COMMENT: Point 1, MAGs cannot distinguish closely related species). While CONCOCT performed significantly worse compared to the other binners, we did notice that CONCOCT seems to display a trend of generating lots of small partial bins. Perhaps CONCOCT bins might be able to distinguish between closely related species to a higher resolution (COMMENT: Small partial bins… what does it mean overall. is this assumption correct? Would it be able to distinguish closely related species?)

While the overall recovery of chromosomes was okay, we were interested in MAG’s ability to correctly bin mobile genetic elements due to their importance in the functions and spread of pathogenic traits such as AMR and virulence. In term of plasmids, a very small proportion of plasmids were correctly binned regardless of the method (<33% at best). Similarly, the same trend exists for genomic islands (<43.3%). This poor result is not unexpected as genomic islands and plasmids have divergent composition features relative to the chromosomes. Furthermore, the difference between the percentages suggest that binning plasmids are harder than GIs. This difference might be due to the problem of plasmid assembly. Therefore, the binning efficiency might improve if we use an assembler targeted at assembling plasmids [@12zFifp5x].

Due to the importance of mobile genetic elements to disseminate clinically relevant antimicrobial resistance genes and virulence factors, we explored whether or not MAGs can be used to provide useful lateral gene transfer insights. 

With respect to AMR genes, MAGs were able to recover roughly half of all AMR genes present in our reference genome. The correct bins were assigned for majority of the chromosomally located AMR genes (81%). The accuracy of chromosomal AMR genes were as expected given the accuracy of MAGs to recover chromosomes as discussed previously. However, while MAGs were able to detect all of plasmid-born AMR genes, none of these were placed in any of the bins. We specifically included a few high threat AMR genes in our dataset: namely KPC and OXA, which are plasmid borne carbapenemases of increasing prevalence in the clinics that are rendering our last resort antibiotics useless. These genes were successfully detected from the metagenomics assembly, but they were not assigned to a bin. This could mean a limited ability for MAGs to be used in the public health research to pinpoint the lateral transfer of AMR genes and to conduct epidemiological analysis (COMMENT: does this make sense?). 

Virulence factors had shown a similar trend as AMR genes, recovering ~60% of virulence factors present in the reference genome. Interestingly, while the detection of virulence factors is better than AMR genes, the binning accuracy was worse, with more being present in the unbinned fraction. Previous studies has found that VFs are disproportionally present on GIs[@LxGqo7iq], which might be the reason to why the binning accuracy was worse compared to AMR genes. 

Lastly, previous works have shown that AMR genes that are on mobile genetic elements disproportionally encode secrete proteins. Given that the recovery of plasmid-borne genes were not great, we asked if MAGs would affect the ability to predict the subcellular localization of proteins. We found that the proportion of predicted localizations were very similar between MAGs and our reference genomes, suggesting that there is not a significant penalty to use MAGs as input for protein localization predictions. 


## Conclusions {#conclusions}

Using a simulated medium complexity metagenome, this study had shown that MAGs provides a great tool to study a bacterial species’ chromosomal elements but presented difficulties in the recovery of mobile genetic elements from metagenomic samples. These mobile genetic elements are liable to being incorrectly binned or lost in this process. Due to the importance of these mobile genomic components in the function and spread of pathogenic traits such as AMR and virulence, it is vital that we utilize a combination of MAGs and other methods (e.g. read-based methods) in public health metagenomic researches. This would allow both the detection of the sample microbial diversity and the thorough evaluation of resistome in metagenomic data to provide meaningful epidemiological information.

## Supplementals {#supplementals}


![Top Species Coverage](images/s1_species_top_coverage.png){#fig:supspeciescov}

