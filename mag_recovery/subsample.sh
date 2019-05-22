
seed=42
for i in 0.05 0.1
do
    
    mkdir -p p${i}
    seqtk sample -s $seed metagenome_trimmed1.fq.gz $i > p${i}/p${i}_r1.fq
    seqtk sample -s $seed metagenome_trimmed2.fq.gz $i > p${i}/p${i}_r2.fq
    echo $seed
    ((seed+=1))
done 


