# trim reads
./sickle/sickle pe -t sanger -f ../../data_simulation/metagenome1.fq.gz -r ../../data_simulation/metagenome2.fq.gz -o r1_trimmed.fq.gz -p r2_trimmed.fq.gz -s s_trimmed.fq.gz

# shuffle to make sure there are no biases from the ART generation
paste - - - - < r1_trimmed.fq > r1_temp.fq
paste - - - - < r2_trimmed.fq > r2_temp.fq
paste r1_temp.fq r2_temp.fq > int_temp.fq
rm r1_temp.fq r2_temp.fq
shuf < int_temp.fq > int_shuf_temp.fq
rm int_temp.fq
tr '\t' '\n' < int_shuf_temp.fq > interleaved_paired_only_trimmed_shuffled.fq
rm int_shuf_temp.fq
