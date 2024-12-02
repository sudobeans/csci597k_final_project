#!/bin/bash
source seq/bin/activate
prefixes=("0M" "0P" "1a" "2P" "3L" "4C")
for pre in "${prefixes[@]}"; do
  bam_files=()
  for bam in sorted_bam/"$pre"*.bam; do
    bam_files+=($bam)
  done
  echo counting ${bam_files[@]} 
  htseq-count --format=bam --order=name --stranded=yes $bam_files Homo_sapiens.GRCh38.113.gtf > htseqResults/"$pre"_output.csv
  echo output to htseqResults/"$pre"_output.csv
done
