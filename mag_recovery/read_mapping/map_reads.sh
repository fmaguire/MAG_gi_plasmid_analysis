for assembly in $(find ../assembly -name 'contig.fa' -o -name 'contigs.fasta' -o -name 'final.contigs.fa')
do
    echo
    echo Building Index
    assembler=$(echo $assembly | cut -d '/' -f 3)
    subsample=$(echo $assembly | cut -d '/' -f 4)
    reads=../subsample/${subsample}.fq
    bowtie2-build $assembly $assembly
    echo Mapping
    bowtie2 -p 4 -x $assembly --interleaved $reads | samtools view -bS - > ${assembler}_${subsample}.bam
    echo Sorting
    samtools sort -O BAM ${assembler}_${subsample}.bam > ${assembler}_${subsample}_sorted.bam
    echo Summarising
    jgi_summarize_bam_contig_depths --outputDepth ${assembler}_${subsample}_depth.txt ${assembler}_${subsample}_sorted.bam 
done


