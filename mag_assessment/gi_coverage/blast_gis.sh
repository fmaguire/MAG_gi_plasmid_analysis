#!/bin/bash

for gis in ../../data/sequences/*/gis/*.fasta
do 
    genome=$(echo $gis | cut -d '/' -f 5 | sed 's/ /_/g')
    sed "s/^>/>$genome:/" $gis >> all_gis.fna
done


for bin in ../bin_blast_dbs/bin.*.nhr
do
    bin_name=$(echo $bin | sed 's/\.nhr//')
    blastn -query all_gis.fna -db $bin_name -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" -out $(basename $bin_name)_gis.out6
done

for i in *_gis.out6; 
do 
    b=$(echo $i | cut -d '_' -f 1)
    sed "s/$/\t$b/" $i >> gs_blast_all.out6
done

