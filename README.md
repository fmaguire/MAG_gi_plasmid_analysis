# MAG Simulation

## Dataset Generation 
### Species Selection

Arbitrarily selected from the following sets:

- 10 with high numbers of plasmids (linked via project ID).
- 10 with high proportion of chromosomes corresponding to composition detected genomic islands (detected using IslandPath-DIMOB).
- 10 with low proportion of chromosomes corresponding to composition detected genomic islands (detected using IslandPath-DIMOB).

### Plasmid Copy Number

Due to a dearth of database resources with known copy numbers each plasmid was randomly selected to be low copy number (1-20), medium copy number (20-100), and high copy number (500-1000).
The exact copy number was randomly selected using a poisson distribution.

### Organism Relative Abundance

- Relative copy number of 

### Read Simulation

`art_illumina` in MSv3 mode with 250bp PE reads at an overall coverage of 1x.

## Metagenome Assembled Genome Pipeline

Use approach specified in Laura Hug's MAG textbook chapter

- Sickle quality control

- fq2fa in idba\_ud sampler to interleave paired reads (drop singletons) 

- subsample 5%, 10%, 25%, 33%, 66%, and 100% of read

- assemble each with idba\_ud

- map reads to contigs for each subsample with with `Bowtie2` and assess coverage of each contig with `BEDtools`

- merge all assemblies with Minimus2 (increasing overlap to 200 and increase minimum percent identity for overlap to 97%)

- Metabat and CONCOCT to do binning and compare

## MAG assessment

...
