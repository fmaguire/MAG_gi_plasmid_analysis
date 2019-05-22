#!/bin/bash

for plasmid in ../../data/sequences/*/*/plasmid/*.fasta
do 
    genome=$( echo $plasmid | cut -d '/' -f 5 | sed 's/ /_/g')
    sed "s/^>/>$genome:/" $plasmid >> all_plasmids.fna
done


for bin in ../bin_blast_dbs/bin.*.nhr
do
    bin_name=$(echo $bin | sed 's/\.nhr//')
    blastn -query all_plasmids.fna -out $(basename $bin_name)_plasmids.out6 -db $bin_name -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" 
done

for i in *_plasmids.out6
do 
    b=$(echo $i | cut -d '_' -f 1)
    sed "s/$/\t$b/" $i >> plasmid_blast_all.out6
done
