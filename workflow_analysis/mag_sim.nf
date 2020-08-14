#!/usr/bin/env nextflow

Channel
    .fromFilePairs( params.reads_pe )
    .ifEmpty{ error "Cannot find any reads matching: ${params.reads_pe}" }
    .set{ metagenome_fastq_pair }

Channel
    .fromPath(params.run_params_csv)
    .splitCsv(header:true)
    .map{ row -> tuple(row.tool, row.label, row.run_params) }
    .branch {
        METASPADES_params: it[0] == "metaspades"
        MEGAHIT_params: it[0] == "megahit"
        IDBA_UD_params: it[0] == "idba_ud"
        CONCOCT_params: it[0] == "concoct"
        MAXBIN2_params: it[0] == "maxbin2"
        METABAT2_params: it[0] == "metabat2"
        DASTOOL_params: it[0] == "dastool"
    }
    .set{ runs_ch }


process trim_reads {
    publishDir "results/0.trimmed_reads", pattern: '*.fq', saveAs: { "${file(it).getSimpleName()}.${file(it).getExtension()}"}
    conda "$baseDir/conda_envs/trim.yml"
    input:
        set val(label), path(reads) from metagenome_fastq_pair
    output:
        path("interleaved_paired_only_trimmed_shuffled_temp.fq") into (MEGAHIT_input, METASPADES_input, IDBA_UD_input, read_mapping_input)
    script:
        """
        sickle pe -t sanger -f ${reads[0]} -r ${reads[1]} -o r1_trimmed.fq -p r2_trimmed.fq -s s_trimmed.fq.gz

        paste - - - - < r1_trimmed.fq > r1_temp.tab
        paste - - - - < r2_trimmed.fq > r2_temp.tab
        paste r1_temp.tab r2_temp.tab > interleaved_temp.tab
        rm r1_temp.tab r2_temp.tab r1_trimmed.fq r2_trimmed.fq

        shuf < interleaved_temp.tab > interleaved_shuffled_temp.tab
        rm interleaved_temp.tab

        tr '\t' '\n' < interleaved_shuffled_temp.tab > interleaved_paired_only_trimmed_shuffled_temp.fq
        """
        //paste - - - - - - - - < interleaved_paired_only_trimmed_shuffled_temp.fq | tee >(cut -f 1-4 | tr "\t" "\n" > r1.fastq) | cut -f 5-8 | tr "\t" "\n" > r2.fastq
}




/* 
    
    Assembly metagenomes

*/ 

runs_ch.METASPADES_params
  .combine( METASPADES_input )
  .set{ METASPADES_run_params }

process run_metaspades {
    publishDir "results/1a_metaspades", pattern: "${label}.fasta"
    conda "$baseDir/conda_envs/assembly.yml"
    input:
        set val(tool), val(label), val(run_params), path(reads) from METASPADES_run_params
    output:
        path("${label}.fasta") into METASPADES_assemblies
    script:
        """
        spades.py ${run_params} -o ${label} --12 $reads -t ${task.cpus}
        mv ${label}/contigs.fasta ${label}.fasta
        """
}


runs_ch.IDBA_UD_params 
    .combine( IDBA_UD_input )
    .set{ IDBA_UD_run_params }

process run_idba_ud {
    publishDir "results/1b_idba_ud", pattern: "${label}.fasta"
    conda "$baseDir/conda_envs/assembly.yml"
    input:
        set val(tool), val(label), val(run_params), path(reads) from IDBA_UD_run_params
    output:
        path("${label}.fasta") into IDBAUD_assemblies
    script:
        """
        fq2fa --paired $reads reads.fasta;
        idba_ud ${run_params} -r reads.fasta -o ${label} --num_threads ${task.cpus}
        mv ${label}/contig.fa ${label}.fasta
        """
}

runs_ch.MEGAHIT_params
    .combine( MEGAHIT_input )
    .set{ MEGAHIT_run_params }


process run_megahit {
    publishDir "results/1c_megahit", pattern: "${label}.fasta" 
    conda "$baseDir/conda_envs/assembly.yml"
    input:
        set val(tool), val(label), val(run_params), path(reads) from MEGAHIT_run_params
    output:
        path("${label}.fasta") into MEGAHIT_assemblies
    script:
        """
        megahit --12 $reads ${run_params} -o ${label} -t ${task.cpus}
        mv ${label}/final.contigs.fa ${label}.fasta
        """
}


/*

    Get contig coverage for assemblies

*/
METASPADES_assemblies
    .concat( IDBAUD_assemblies, MEGAHIT_assemblies )
    .combine( read_mapping_input )
    .set{ assembly_coverage_input }

process coverage_analysis {
    publishDir "results/2_coverages", pattern: "*.depth.tsv"
    conda "$baseDir/conda_envs/coverage.yml"
    input:  
        set path(assembly), path(reads) from assembly_coverage_input
    output:
        set path(assembly), path("*.depth.tsv") into (METABAT2_input, MAXBIN2_input, CONCOCT_input)
    script:
        """
        bowtie2-build ${assembly} ${assembly}
        bowtie2 -p ${task.cpus} -x ${assembly} --interleaved ${reads} | samtools view -bS - > assembly.bam
        samtools sort -@ ${task.cpus} -O BAM assembly.bam > assembly_sorted.bam
        jgi_summarize_bam_contig_depths --outputDepth ${assembly.baseName}.depth.tsv assembly_sorted.bam
        """
}

/*

    Binning

*/

//runs_ch.METABAT2_params
//    .combine(METABAT2_input )
//    .set{ METABAT2_run_params }
//
//process run_metabat2{
//    publishDir "results/3a_metabat2/${assembly.baseName}_${label}", pattern: "${assembly.baseName}_${label}_bins.*.fa"
//    conda "$baseDir/conda_envs/binners.yml"
//    input:
//        set val(tool), val(label), val(run_params), path(assembly), path(coverage) from METABAT2_run_params
//    output:
//        path("${assembly.baseName}_${label}*.fa") into METABAT2_bins
//    script:
//        """
//        metabat2 -i ${assembly} -a ${coverage} -o ${assembly.baseName}_${label} --unbinned --seed 42 -t ${task.cpus} ${run_params}
//        mv ${assembly.baseName}_${label}.unbinned.fa ${assembly.baseName}_${label}.unbinned.fasta
//        mv ${assembly.baseName}_${label}.tooShort.fa ${assembly.baseName}_${label}.tooShort.fasta
//        mv ${assembly.baseName}_${label}.lowDepth.fa ${assembly.baseName}_${label}.lowDepth.fasta
//        touch ${assembly.baseName}_${label}.0.fa
//        """
//}
//
//
//runs_ch.MAXBIN2_params
//    .combine( MAXBIN2_input )
//    .set{ MAXBIN2_run_params }
//
//process run_maxbin2{
//    publishDir "results/3b_maxbin2/${assembly.baseName}_${label}", pattern: "${assembly.baseName}_${label}.*.fasta"
//    conda "$baseDir/conda_envs/binners.yml"
//    input:
//        set val(tool), val(label), val(run_params), path(assembly), path(coverage) from MAXBIN2_run_params
//    output:
//        path("${assembly.baseName}_${label}.*.fasta") into MAXBIN2_bins
//    script:
//        """
//        run_MaxBin.pl ${run_params} -thread ${task.cpus} -abund ${coverage} -contig ${assembly} -out ${assembly.baseName}_${label}
//        """
//}

//runs_ch.CONCOCT_params
//    .combine( CONCOCT_input )
//    .set{ CONCOCT_run_params }
//
//process run_concoct {
//    publishDir "results/3c_concoct"
//    conda "$baseDir/conda_envs/binners.yml"
//    input:
//        set val(tool), val(label), val(run_params), path(assembly), path(coverage) from CONCOCT_run_params
//    output: 
//        path(bins) into CONCOCT_bins
//    script:
//        """
//        cut_up_fasta.py -c 10000 -o 0 --merge_last -b ${assembly.baseName}_contigs10k.bed ${assembly} > ${assembly.baseName}_contigs10k.fa
//        samtools index ${coverage}
//        concoct_coverage_table.py ${assembly.baseName}_contigs10k.bed ${coverage} > concoct_coverage.tsv
//        concoct ${run_params} -r 250 -t ${task.cpu} --coverage_file concoct_coverage.tsv -s 42 --composition_file ${assembly.baseName}_contigs10k.fa -b concoct_output
//        merge_cutup_clustering.py concoct_output/*_gt1000.csv > all_clustering_merged.csv
//        extract_fasta_bins.py --output_path ${assembly.baseName}_${label} ${assembly} all_clustering_merged.csv
//        """
//}

//runs_ch.DASTOOL_params
//    .combine(CONCOCT_bins, 
//             MAXBIN2_bins,
//             METABAT2_bins)
//
//process run_DASTOOL {
//    publishDir "results/3c_concoct"
//    conda "$baseDir/conda_envs/binners.yml"
//
//        mkdir -p dastool/${assembler}/${assembler}
//        # to prevent unbinned being included as a bin
//        mv metabat2/${assembler}/${assembler}_all.unbinned.fa metabat2/${assembler}/${assembler}_all.unbinned.fasta
//        Fasta_to_Scaffolds2Bin.sh -i metabat2/${assembler} -e 'fa' > metabat2/${assembler}/metabat2_scaffolds2bin.tsv
//        Fasta_to_Scaffolds2Bin.sh -i maxbin2/${assembler} -e 'fasta' > maxbin2/${assembler}/maxbin2_scaffolds2bin.tsv
//        Fasta_to_Scaffolds2Bin.sh -i concoct/${assembler} -e 'fa' > concoct/${assembler}/concoct_scaffolds2bin.tsv
//
//        DAS_Tool -i metabat2/${assembler}/metabat2_scaffolds2bin.tsv,maxbin2/${assembler}/maxbin2_scaffolds2bin.tsv,concoct/${assembler}/concoct_scaffolds2bin.tsv -c ${assembly} --write_bins --write_unbinned -l metabat2,maxbin2,concoct --search_engine diamond -t 4 -o dastool/${assembler}/${assembler}
//
//        mv metabat2/${assembler}/${assembler}_all.unbinned.fasta metabat2/${assembler}/${assembler}_all.unbinned.fa
//
//
//}
