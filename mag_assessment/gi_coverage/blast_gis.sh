#!/bin/bash

for gis in ../../data/sequences/*/gis/*.fasta
do 
    genome=$(echo $gis | cut -d '/' -f 5 | sed 's/ /_/g')
    sed "s/^>/>$genome:/" $gis >> all_gis.fna
done


mkdir -p {maxbin2,metabat2,concoct,dastool}/{idba_ud,metaspades,megahit}
for binning_tool in maxbin2 metabat2 concoct dastool; do
    for assembly in idba_ud metaspades megahit; do 
        bins=(../bin_blast_dbs/${binning_tool}/${assembly}/*.nhr)
        
        echo $binning_tool $assembly
        echo
      
        for bin in "${bins[@]}"; do
            bin_name=$(echo $bin | sed 's/\.nhr//')
            blastn -query all_gis.fna -db $bin_name -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" -out ${binning_tool}/${assembly}/$(basename $bin_name)_gis.out6
        done

        for bin_hits in ${binning_tool}/${assembly}/*_gis.out6
        do 
            binname=$(basename $bin_hits)
            sed "s/$/\t$binning_tool\t$assembly\t$binname/" $bin_hits >> gis_blast_all.out6
        done 
    done
done 

