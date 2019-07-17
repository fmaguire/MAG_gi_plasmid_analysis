#!/bin/bash
set -e 

# just do the 100% assemblies
for assembly in ../assembly/idba_ud/perc100/contig.fa ../assembly/megahit/perc100/final.contigs.fa ../assembly/metaspades/perc100/contigs.fasta
    do
        assembler=$(echo $assembly | cut -d '/' -f 3)
        subsample=$(echo $assembly | cut -d '/' -f 4)
        depth="../read_mapping/${assembler}_${subsample}_depth.txt"
        mkdir -p {concoct,metabat2,maxbin2}/${assembler}
        metabat2 -i $assembly -o metabat2/${assembler}/${assembler}_all --unbinned --sed-arg 42 -a $depth 
        concoct --coverage_file $depth -s 42 --composition_file $assembly -b concoct/${assembler}/${assembler}_all
        run_MaxBin.pl -abund $depth -contig $assembly -out maxbin2/${assembler}/${assembler}_all
    done
