#!/bin/bash
source sid/bin/activate
prefixes=("0M" "0P" "1a" "2P" "3L" "4C")
for pre in "${prefixes[@]}"; do
  for bam in sorted_bam/"$pre"*.bam; do
  echo counting ${bam}
  htseq-count --format=bam --order=name --stranded=yes $bam Homo_sapiens.GRCh38.113.gtf > /home/bhimirs/Documents/csci474/htseqCounts/"$pre"/"$(basename "$bam")"_output.csv
  echo "Output to htseqResults/${pre}/$(basename "$bam")_output.csv"
  done
done
