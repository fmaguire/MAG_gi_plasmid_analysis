mkdir -p bin_blast_dbs/{maxbin2,metabat2,concoct,dastool}/{idba_ud,metaspades,megahit}



for binning_tool in maxbin2 metabat2 concoct dastool; do
    for assembly in idba_ud metaspades megahit; do 
            
        if [ "$binning_tool" == "maxbin2" ]; then
            bins=(../mag_recovery/binning/${binning_tool}/${assembly}/*.fasta)
        elif [ "$binning_tool" == "metabat2" ]; then
            bins=(../mag_recovery/binning/${binning_tool}/${assembly}/*.fa)
        elif [ "$binning_tool" == "concoct" ]; then
            bins=(../mag_recovery/binning/${binning_tool}/${assembly}/*.fa)
        elif [ "$binning_tool" == "dastool" ]; then
            bins=(../mag_recovery/binning/${binning_tool}/${assembly}/${assembly}_DASTool_bins/*.fa)
        else
            echo "Error"
            exit 1
        fi
        
        echo $binning_tool
        echo $assembly
      
        for bin in "${bins[@]}"; do
            makeblastdb -title ${binning_tool}_${assembly}_$(basename $bin) -dbtype nucl -in $bin -out bin_blast_dbs/${binning_tool}/${assembly}/$(basename $bin)
        done

    done
done 


