#!/bin/bash

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
