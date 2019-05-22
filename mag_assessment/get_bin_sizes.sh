#!/bin/bash

# count nt in each bin
echo "bin,nts" > bin_sizes.csv
for bin in ../mag_recovery/metabat_bins/bins/*.fa;
do
    echo $(basename $bin),$(grep -v ">" $bin | wc | awk '{print $3-$1}') >> bin_sizes.csv
done

# get ref chromosome genome sizes
echo "genome,nts" > genome_sizes.csv
for genome in ../data/sequences/*/*/chromosome/*.fasta
do
    genome_name=$(echo $genome | cut -d '/' -f 4)
    echo $genome_name,$(grep -v ">" $genome | wc | awk '{print $3-$1}') >> genome_sizes.csv
done
