#!/bin/bash

#bins
for binning_tool in maxbin2 metabat2 concoct dastool; do
	for assembly_tools in idba_ud metaspades megahit; do
		if [ "$binning_tool" == "dastool" ]; then
			echo $binning_tool/$assembly_tools;
			for file in `ls $binning_tool/$assembly_tools/*DASTool_bins/*.fasta`; do
				echo $file;
				bname=`basename $file`;
				#prodigal -i $file -a $file.protein 
				count=`cat $file.protein | grep ">" | wc -l`
				echo "$bname","$count","bin","$binning_tool","$assembly_tools" >> geneCounts.csv
			done;
		else
			echo $binning_tool/$assembly_tools;
			for file in `ls $binning_tool/$assembly_tools/*.fa`; do
				echo $file;
				bname=`basename $file`;
				prodigal -i $file -a $file.protein 
				count=`cat $file.protein | grep ">" | wc -l`
				echo "$bname","$count","bin","$binning_tool","$assembly_tools" >> geneCounts.csv
			done
		fi
	done;
done;


