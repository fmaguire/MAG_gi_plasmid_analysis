#!/bin/bash

for chrom in ../../data/sequences/*/*/chromosome/*.fasta
do 
    genome=$( echo $chrom | cut -d '/' -f 5 | sed 's/ /_/g')
    sed "s/^>/>$genome:/" $chrom>> all_chroms.fna
done

makeblastdb -in all_chroms.fna -out all_chroms -dbtype nucl

for bin in ../../mag_recovery/metabat_bins/bins/*.fa 
do
    blastn -query $bin -db all_chroms -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" -out $(basename $bin)_chroms.out6
done

for i in *_chroms.out6
do 
    b=$(echo $i | cut -d '_' -f 1)
    sed "s/$/\t$b/" $i >> chroms_blast_all.out6
done 
