#!/bin/bash -l
#
#SBATCH --workdir=/home/jjjjia/scratch/magassessment/MAG_gi_plasmid_analysis/mag_recovery/
#SBATCH --account=rrg-fiona-ad
#SBATCH --job-name=MAGassessment
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=8:00:00
#SBATCH --mail-user=<bja20@sfu.ca>
#SBATCH --mail-type=ALL

#!/bin/bash 

#sickle pe -g -f metagenome1.fq.gz -r metagenome2.fq.gz -t sanger -o metagenome_trimmed1.fq.gz -p metagenome_trimmed2.fq.gz -s metagenome_trimmed_s.fq.gz

## don't know if this supports gzip if not add <(zcat file) instead
#fq2fa --merge metagenome_trimmed1.fq.gz metagenome_trimmed2.fq.gz metagenome_trimmed.fna

#
## shuffle
#shuffle.sh in=metagenome_trimmed.fna interleaved=t out=metagenome_trimmed_shuffled.fna
#
## subsets fastas so divisable by 4
#for subset in p5 p10 p25 p35 p66 p100
#do
#    mkdir -p $subset
#    head -n ... metagenome_trimmed_shuffled.fna | tail -n ... > ${subset}/${subset}_metagenome.fna
#    
#    idba_ud -r ${subset}/${subset}_metagenome.fna -o ${subset}
#    ref=${subset}/contig.fa
#
#    bowtie2-build ${ref}
#    bowtie2 -p $THREADS -x ${ref} -1 ${subset}/${subset}_metagenome.fna -S ${subset}/${subset}.sam
#
#    samtools view -bt ${ref} ${subset}/${subset}.sam > ${subset}/${subset}.bam
#    samtools sort ${subset}/${subset}.bam ${subset}/${subset}_sorted.bam
#    samtools index ${subset}/${subset}_sorted.bam
#
#    # not bothering to remove duplicates with picard
#
#    genomeCoverageBed -ibam ${subset}/${subset}_sorted.bam > ${subset}/${subset}_coverage
#    # might need to use concoct script gen_contig_cov_per_bam_table.py
#done
#
#mkdir -p merges
#minimus2 merges -D OVERLAP=200 -D MINID=97 ...
##	-D REFCOUNT=<n>         \       # Number of sequences in the 1st assembly ; (Required)
#
##bowtie2 for merged assembly



# non-gradated
mkdir -p p100/bins
idba_ud --num_threads 48 -r metagenome_trimmed.fna -o p100
ref=p100/contig.fa

bowtie2-build -f --threads 48 $ref $ref

bowtie2 -p 48 -f --sensitive-local -x p100/contigs.fa --interleaved metagenome_trimmed.fna -S p100/p100.sam
samtools view -@ 48 -S -b p100/p100.sam > p100/p100.bam
samtools sort -@ 48 -o p100/p100_sorted.bam p100/p100.bam

runMetaBat.sh -m 1500 -t 48 "$ref" p100/p100_sorted.bam 
##concoct
#
