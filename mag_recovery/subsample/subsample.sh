total=205461152

head -n 1027296 ../trimmed/interleaved_paired_only_trimmed_shuffled.fq > perc5.fq
head -n 20546115 ../trimmed/interleaved_paired_only_trimmed_shuffled.fq > perc10.fq
head -n 51365288 ../trimmed/interleaved_paired_only_trimmed_shuffled.fq > perc25.fq
head -n 67802180 ../trimmed/interleaved_paired_only_trimmed_shuffled.fq > perc33.fq
head -n 135604360 ../trimmed/interleaved_paired_only_trimmed_shuffled.fq > perc66.fq
cat ../trimmed/interleaved_paired_only_trimmed_shuffled > perc100.fq
