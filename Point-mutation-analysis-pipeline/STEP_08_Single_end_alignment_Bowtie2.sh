#!/bin/bash

# Loop over all files ending with _filtered.fastq.gz in the current directory
for fqfile in *_filtered.fastq.gz; do
    # Extract the base name by removing the suffix _filtered.fastq.gz
    base=$(basename "$fqfile" _filtered.fastq.gz)

    # Run bowtie2 single-end alignment
    bowtie2 -x /Path_to_index/hg38_bt2_index/hg38 \
        -U "$fqfile" -t -S "${base}_filtered.sam" \
        --no-unal -p 20 2> "${base}_filtered_alignment_stats.txt"
done
