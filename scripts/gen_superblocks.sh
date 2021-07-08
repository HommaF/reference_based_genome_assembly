#!/bin/bash

paired_mapped_int=$1
tmp_blocks_file=$2

zcat $paired_mapped_int | awk '{if ($4>9) print $0 "\tblock"; else print $0 "\tgap"}' | awk -v gap=0 -v len=0 -v chrom=0 -v prev_type=0 -v prev_len=0 '{if (chrom == 0) chrom = $1; else if (chrom != $1 && prev_type == "gap") {chrom = $1; len=0;} else if (chrom != $1 && prev_type == "block") {print chrom "\tblock\t" len "\t" prev_len "\t" prev_len-len; chrom=$1; len=0; gap=0;} if ($5 == "block" && gap == 0 && chrom == $1) {print $1 "\tgap\t" len "\t" $2 "\t" $2-len; gap=1; len=$2;} else if ($5 == "gap" && gap == 1 && chrom == $1) {print $1 "\tblock\t" len "\t" $2 "\t" $2-len; gap=0; len=$2;} if ($5 == "gap") {prev_type="gap"} else if ($5 == "block") {prev_type="block"; prev_len=$3;}}' > $tmp_blocks_file
