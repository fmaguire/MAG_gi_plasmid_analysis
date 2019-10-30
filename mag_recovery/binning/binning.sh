#!/bin/bash
set -e 

# just do the 100% assemblies
for assembly in ../assembly/idba_ud/perc100/contig.fa ../assembly/megahit/perc100/final.contigs.fa ../assembly/metaspades/perc100/contigs.fasta
    do
        assembler=$(echo $assembly | cut -d '/' -f 3)
        subsample=$(echo $assembly | cut -d '/' -f 4)
        depth="../read_mapping/${assembler}_${subsample}_depth.txt"
        mkdir -p {dastool,concoct,metabat2,maxbin2}/${assembler}
        # metabat2
        metabat2 -i $assembly -o metabat2/${assembler}/${assembler}_all --unbinned --sed-arg 42 -a $depth 

        # maxbin2
        run_MaxBin.pl -abund $depth -contig $assembly -out maxbin2/${assembler}/${assembler}_all

        #concoct
        mkdir -p concoct/${assembler}/run_files
        python CONCOCT_binary/scripts/cut_up_fasta.py -c 10000 -o 0 --merge_last -b ${assembly}_contigs10k.bed ${assembly} > ${assembly}_contigs10k.fa
        samtools index ../read_mapping/${assembler}_${subsample}_sorted.bam  
        python CONCOCT_binary/scripts/concoct_coverage_table.py ${assembly}_contigs10k.bed ../read_mapping/${assembler}_${subsample}_sorted.bam > ../read_mapping/${assembler}_${subsample}_coverage.tsv
        concoct --coverage_file ../read_mapping/${assembler}_${subsample}_coverage.tsv -s 42 --composition_file ${assembly}_contigs10k.fa -b concoct/${assembler}/run_files/${assembler}_all
        python CONCOCT_binary/scripts/merge_cutup_clustering.py concoct/${assembler}/run_files/${assembler}_all_clustering_gt1000.csv > concoct/${assembler}/run_files/${assembler}_all_clustering_merged.csv
        python CONCOCT_binary/scripts/extract_fasta_bins.py --output_path concoct/${assembler} ${assembly} concoct/${assembler}/run_files/${assembler}_all_clustering_merged.csv 

        # das_tool to combine predictions from previous tools
        mkdir -p dastool/${assembler}/${assembler}
        # to prevent unbinned being included as a bin
        mv metabat2/${assembler}/${assembler}_all.unbinned.fa metabat2/${assembler}/${assembler}_all.unbinned.fasta
        Fasta_to_Scaffolds2Bin.sh -i metabat2/${assembler} -e 'fa' > metabat2/${assembler}/metabat2_scaffolds2bin.tsv
        Fasta_to_Scaffolds2Bin.sh -i maxbin2/${assembler} -e 'fasta' > maxbin2/${assembler}/maxbin2_scaffolds2bin.tsv
        Fasta_to_Scaffolds2Bin.sh -i concoct/${assembler} -e 'fa' > concoct/${assembler}/concoct_scaffolds2bin.tsv

        DAS_Tool -i metabat2/${assembler}/metabat2_scaffolds2bin.tsv,maxbin2/${assembler}/maxbin2_scaffolds2bin.tsv,concoct/${assembler}/concoct_scaffolds2bin.tsv -c ${assembly} --write_bins --write_unbinned -l metabat2,maxbin2,concoct --search_engine diamond -t 4 -o dastool/${assembler}/${assembler}

        mv metabat2/${assembler}/${assembler}_all.unbinned.fasta metabat2/${assembler}/${assembler}_all.unbinned.fa
    done
