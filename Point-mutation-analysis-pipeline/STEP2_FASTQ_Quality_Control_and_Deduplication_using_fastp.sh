#!/bin/bash

# Loop over all files ending with _R1.fastq.gz in the current directory
for sample in *_R1.fastq.gz; do
    # Extract the base name by removing the suffix _R1.fastq.gz
    base=$(basename "$sample" _R1.fastq.gz)

    # Run fastp for paired-end quality control and deduplication
    fastp -i "${base}_R1.fastq.gz" -I "${base}_R2.fastq.gz" \
          -o "${base}_R1_fastp.fastq.gz" -O "${base}_R2_fastp.fastq.gz" \
          -D -w 16
done
