#!/bin/bash
# convert SAM to BAM and sort by name
for dir in results/*/; do
  echo dir $dir
  if [ -d "$dir" ]; then
    for subdir in "$dir"/*/; do
      subdir_name=$(basename "$subdir")
      if [ -f "$subdir/Aligned.out.sam" ]; then
        sam_file="$subdir/Aligned.out.sam"
        group="${dir:9:3}"
        short="${subdir_name:0:12}"
        output=$group-$short
        ./sambamba view -S --format=bam $sam_file | ./sambamba sort -n -o sorted_bam/${output}.bam /dev/stdin
      fi
    done
  fi
done

# generate one count matrix per group
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
