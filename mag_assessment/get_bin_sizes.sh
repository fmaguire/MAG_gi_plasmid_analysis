#!/bin/bash

# count nt in each bin
echo "assembler,binner,bin,nts" > bin_sizes.csv


for binning_tool in maxbin2 metabat2 concoct dastool; do
    for assembly in idba_ud metaspades megahit; do 
            
        if [ "$binning_tool" == "maxbin2" ]; then
            bins=(../mag_recovery/binning/${binning_tool}/${assembly}/*.fasta)
        elif [ "$binning_tool" == "metabat2" ]; then
            bins=(../mag_recovery/binning/${binning_tool}/${assembly}/*.fa)
        elif [ "$binning_tool" == "concoct" ]; then
            bins=(../mag_recovery/binning/${binning_tool}/${assembly}/*.fa)
        elif [ "$binning_tool" == "dastool" ]; then
            bins=(../mag_recovery/binning/${binning_tool}/${assembly}/${assembly}_DASTool_bins/*.fa)
        else
            echo "Error"
            exit 1
        fi
        
        echo $binning_tool
        echo $assembly
      
        for bin in "${bins[@]}"; do
            echo $assembly,$binning_tool,$(basename $bin),$(grep -v ">" $bin | wc | awk '{print $3-$1}') >> bin_sizes.csv
        done

    done
done 

# get ref chromosome genome sizes
echo "genome,nts" > genome_sizes.csv
for genome in ../data/sequences/*/*/chromosome/*.fasta
do
    genome_name=$(echo $genome | cut -d '/' -f 4)
    echo $genome_name,$(grep -v ">" $genome | wc | awk '{print $3-$1}') >> genome_sizes.csv
done
