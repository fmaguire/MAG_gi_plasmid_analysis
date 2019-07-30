#!/bin/bash

#extract the unbinned for each binner and assembly

for assembly in ../assembly/idba_ud/perc100/contig.fa ../assembly/megahit/perc100/final.contigs.fa ../assembly/metaspades/perc100/contigs.fasta
    do
    for binning_tool in maxbin2 metabat2 concoct dastool; do

        assembler=$(echo $assembly | cut -d '/' -f 3)
        if [ "$binning_tool" == "dastool" ]; then
            cut -f 1 ${binning_tool}/${assembler}/*scaffolds2bin.txt | sort > binned_contigs
        else
            cut -f 1 ${binning_tool}/${assembler}/*scaffolds2bin.tsv | sort > binned_contigs
        fi

        grep "^>" ${assembly} | sed 's/^>//' | sort > all_contigs
        comm -23 all_contigs binned_contigs | sed 's/^/>/' > unbinned_contigs
        fasta_formatter -i ${assembly} > clean_assembly

        if [ "$binning_tool" == "maxbin2" ]; then
            grep -A1 -f unbinned_contigs clean_assembly > ${binning_tool}/${assembler}/unbinned.fasta
        elif [ "$binning_tool" == "metabat2" ]; then
            grep -A1 -f unbinned_contigs clean_assembly > ${binning_tool}/${assembler}/unbinned.fa
        elif [ "$binning_tool" == "concoct" ]; then
            grep -A1 -f unbinned_contigs clean_assembly > ${binning_tool}/${assembler}/unbinned.fa
        elif [ "$binning_tool" == "dastool" ]; then
            grep -A1 -f unbinned_contigs clean_assembly > ${binning_tool}/${assembler}/${assembler}_DASTool_bins/unbinned.fa
        else
            echo "Error"
            exit 1
        fi
        rm all_contigs binned_contigs unbinned_contigs clean_assembly
    done
done
