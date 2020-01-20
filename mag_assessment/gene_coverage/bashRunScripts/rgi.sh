#!/bin/bash -e

for binning_tool in maxbin2 metabat2 concoct dastool; do
        for assembly_tools in idba_ud metaspades megahit; do
                if [ "$binning_tool" == "dastool" ]; then
                        echo $binning_tool/$assembly_tools;
                        for file in `ls $binning_tool/$assembly_tools/*DASTool_bins/*.fa*`; do
                                echo $file;
                                mkdir -p $file.rgi;
                                bname=`basename $file`;
#                               rgi main -i $file -o $file.rgi/$bname -t contig -a DIAMOND -n 16 --include_loose --local;
                        done;
                else
                        echo $binning_tool/$assembly_tools;
                        for file in `ls $binning_tool/$assembly_tools/*.fa*`; do
                                echo $file;
                                mkdir -p $file.rgi;
                                bname=`basename $file`;
                                rgi main -i $file -o $file.rgi/$bname -t contig -a DIAMOND -n 1 --include_loose --local;
                                exit
                        done
                fi
        done;
done;