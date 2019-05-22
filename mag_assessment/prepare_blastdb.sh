mkdir -p bin_blast_dbs
for i in ../mag_recovery/metabat_bins/bins/*.fa; 
do 
    makeblastdb -dbtype nucl -in $i -out bin_blast_dbs/$(basename $i)
done

