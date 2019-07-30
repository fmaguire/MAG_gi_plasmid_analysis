#!/bin/bash

for chrom in ../../data/sequences/*/*/chromosome/*.fasta
do 
    genome=$( echo $chrom | cut -d '/' -f 5 | sed 's/ /_/g')
    sed "s/^>/>$genome:/" $chrom>> all_chroms.fna
done

makeblastdb -in all_chroms.fna -out all_chroms -dbtype nucl


mkdir -p {maxbin2,metabat2,concoct,dastool}/{idba_ud,metaspades,megahit}
for binning_tool in maxbin2 metabat2 concoct dastool; do
    for assembly in idba_ud metaspades megahit; do 
            
        if [ "$binning_tool" == "maxbin2" ]; then
            bins=(../../mag_recovery/binning/${binning_tool}/${assembly}/*.fasta)
        elif [ "$binning_tool" == "metabat2" ]; then
            bins=(../../mag_recovery/binning/${binning_tool}/${assembly}/*.fa)
        elif [ "$binning_tool" == "concoct" ]; then
            bins=(../../mag_recovery/binning/${binning_tool}/${assembly}/*.fa)
        elif [ "$binning_tool" == "dastool" ]; then
            bins=(../../mag_recovery/binning/${binning_tool}/${assembly}/${assembly}_DASTool_bins/*.fa)
        else
            echo "Error"
            exit 1
        fi
        
        echo $binning_tool $assembly
        echo
      
        for bin in "${bins[@]}"; do
            blastn -query $bin -num_hits 1 -db all_chroms -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" -out ${binning_tool}/${assembly}/$(basename $bin)_chroms.out6
        done

        for bin_hits in ${binning_tool}/${assembly}/*_chroms.out6
        do 
            binname=$(basename $bin_hits)

            sed "s/$/\t$binning_tool\t$assembly\t$binname/" $bin_hits >> chroms_blast_all.out6
        done 
    done
done 



