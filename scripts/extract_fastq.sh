#!/bin/bash

number=$1
unm_ids=$2
fwd_paired=$3
rev_paired=$4
unm_fwd_out=$5
unm_rev_out=$6

if [ $number == 1 ]; then
	seqkit grep --pattern-file $unm_ids $fwd_paired | bgzip > $unm_fwd_out

elif [ $number == 2 ]; then
	seqkit grep --pattern-file $unm_ids $rev_paired | bgzip > $unm_rev_out
fi
