#!/bin/bash

input_file=$1

tail -n +2 $input_file | xargs -I@ -n 1 -P 12 bash -c 'wdir=$(echo @ | cut -d " " -f7); name=$(echo @ | cut -d " " -f6); chrom=$(echo @ | cut -d " " -f1); start=$(echo @ | cut -d " " -f2); end=$(echo @ | cut -d " " -f3); samtools view -b results/bwa_aligned/*sorted_filtered.bam $chrom:$start-$end > $wdir$name.bam'

#tail -n +2 $input_file | xargs -I@ -n 1 -P 12 bash -c 'wdir=$(echo @ | cut -d " " -f7); name=$(echo @ | cut -d " " -f6); chrom=$(echo @ | cut -d " " -f1); start=$(echo @ | cut -d " " -f2); end=$(echo @ | cut -d " " -f3); samtools view results/bwa_aligned/*sorted_filtered.bam $chrom:$start-$end | awk -v same="=" "{if (\$7 == same) print}" | cut -f1 | sort | uniq > $wdir"tmp_"$name"_"ids.txt'

