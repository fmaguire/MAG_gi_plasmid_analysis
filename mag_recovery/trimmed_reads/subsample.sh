total=205461152

head -n 1027296 interleaved_paired_only_trimmed_shuffled.fq > perc5.fq
head -n 20546115 interleaved_paired_only_trimmed_shuffled.fq > perc10.fq
head -n 51365288 interleaved_paired_only_trimmed_shuffled.fq > perc25.fq
head -n 67802180 interleaved_paired_only_trimmed_shuffled.fq > perc33.fq
head -n 135604360 interleaved_paired_only_trimmed_shuffled.fq > perc66.fq
cat interleaved_paired_only_trimmed_shuffled > perc100.fq
