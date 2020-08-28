
[![DOI](https://zenodo.org/badge/179372262.svg)](https://zenodo.org/badge/latestdoi/179372262)


# MAG Simulation

## Dataset Generation 
### Species Selection

Arbitrarily selected from the following sets:

- 10 with high numbers of plasmids (linked via project ID).
- 10 with high proportion of chromosomes corresponding to composition detected genomic islands (detected using IslandPath-DIMOB).
- 10 with low proportion of chromosomes corresponding to composition detected genomic islands (detected using IslandPath-DIMOB).

The data used to select the taxa is listed in `data/gi_plasmid_analysis.tsv` and the details of the select taxa are listed in `data/taxa_metadata.tsv`
The sequences themselves are in `data/sequences/sequences.tar.bz2`

### Plasmid Copy Number

Due to a dearth of database resources with known copy numbers each plasmid was randomly selected to be low copy number (1-20), medium copy number (20-100), and high copy number (500-1000).
The exact copy number was randomly selected using a gamma distribution parameterised towards the lower bound for each regime.


### Organism Relative Abundance

Relative abundance of organisms was selected according to the log-normal distribution.

### Read Simulation

For chromosome the underlying fasta was copied the corresponding number of times before read simulation as per the relative abundance estimate.
For the plasmids the abundance of the host taxa was multipled by the copy number and the underlying sequence copied the corresponding number of times.

`art_illumina` in MSv3 mode with 250bp PE reads at an overall coverage of 1x was used to simulated the actual read data.

This and the previous steps were performed using the `data_simluation/simulate_metagenome.py` script

## Metagenome Assembled Genome Pipeline

Use approach specified in Laura Hug's MAG textbook chapter

- Sickle quality control

- fq2fa in idba\_ud sampler to interleave paired reads (drop singletons) 

- subsample 5%, 10%, 25%, 33%, 66%, and 100% of read

- assemble each with idba\_ud

- map reads to contigs for each subsample with with `Bowtie2` and assess coverage of each contig with `BEDtools`

- merge all assemblies with Minimus2 (increasing overlap to 200 and increase minimum percent identity for overlap to 97%)

- Metabat and CONCOCT to do binning and compare
 (Kang et al. provide a
reasonable comparison of these methods)

## MAG assessment

- Annotate recovered genomes and original genomes with prokka on Vanilla

- predict GI from recovered genome with islandpath-DIMOB

- create list of genes in GI and non-GI of each genome

- look at proportion of missed genes in GI and non-GI of each genome


