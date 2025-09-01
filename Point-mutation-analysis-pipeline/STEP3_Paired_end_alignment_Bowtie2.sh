#!/bin/bash

# Loop over all files ending with _R1_filtered.fastq.gz in the current directory
for fqfile in *_R1_filtered.fastq.gz; do
    # Extract the base name by removing the suffix _R1_filtered.fastq.gz
    base=$(basename "$fqfile" _R1_filtered.fastq.gz)

    # Run bowtie2 paired-end alignment
    bowtie2 -x /Path_to_index/hg38_bt2_index/hg38 \
        -1 "${base}_R1_filtered.fastq.gz" -2 "${base}_R2_filtered.fastq.gz" \
        -t -S "${base}_filtered.sam" --no-unal -p 20 2> "${base}_filtered_alignment_stats.txt"
done
