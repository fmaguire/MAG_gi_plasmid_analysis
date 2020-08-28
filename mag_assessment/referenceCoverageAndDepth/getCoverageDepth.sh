#!/bin/bash -e

#
#SBATCH --chdir=/scratch/jjjjia/mag_bowtie2/genomecov
#SBATCH --account=rrg-fiona-ad
#SBATCH --job-name=run_samtools
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=00:20:00
#SBATCH --mail-user=<bja20@sfu.ca>
#SBATCH --mail-type=ALL

date;

source activate samtools


echo $1;

if [[ -f ../fastas/chromosome/$1 ]]; then
	echo "chromosome"
	reference="../fastas/chromosome/$1"
	outPath="chromosome";
	echo $reference
elif [[ -f ../fastas/plasmid/$1 ]]; then
	echo "plasmid"
	outPath="plasmid"
	reference="../fastas/plasmid/$1";
	echo $reference
else
	reference="../fastas/$1";
	echo "$reference";
	outPath="combined"
fi


echo $outPath;
mkdir -p $outPath;

bam=$1.sam.sorted.bam;
echo $reference;
echo $bam;

echo "depth";

samtools depth -a ../alignments/bam/$bam > $outPath/$bam.depth


echo "coverage"
samtools coverage --reference $reference -o $outPath/$bam.coverage ../alignments/bam/$bam;

echo "coverage histogram"
samtools coverage -m --reference $reference -o $outPath/$bam.coverage.histogram ../alignments/bam/$bam

source deactivate;

date;
