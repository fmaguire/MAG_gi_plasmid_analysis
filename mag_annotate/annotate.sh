#!/bin/bash -l
#
#SBATCH --workdir=/home/jjjjia/scratch/magassessment/MAG_gi_plasmid_analysis/mag_annotate
#SBATCH --account=rrg-fiona-ad
#SBATCH --job-name=mag_annotate
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=23:50:00
#SBATCH --mail-user=<bja20@sfu.ca>
#SBATCH --mail-type=ALL


#annotate existing genomes from refseq in ../data/sequences/*

OIFS="$IFS"
IFS=$'\n'
 
for i in `ls ../data/sequences`;
do
	#$i is the spp name
	echo "$i";
	for k in `ls ../data/sequences/"$i"/*/`;
	do 
		echo "$k" #this is either chr or plasmid
	
		for j in `ls ../data/sequences/"$i"/*/"$k"/*.fasta`;
			#$j is the fasta name
			do
			echo $j;
			prokka --cpu 18 --outdir "$j.annotation" "$j";
		done;
	done;
done;

IFS="$OIFS"
