#!/bin/bash

# Loop over all files ending with _filtered.sam in the current directory
for samfile in *_filtered.sam; do
    # Extract the base name by removing the suffix _filtered.sam
    base=$(basename "$samfile" _filtered.sam)
    
    # Convert SAM to BAM using 20 threads
    samtools view -@ 20 -bS "${samfile}" > "${base}.bam"
    
    # Sort BAM file using 20 threads
    samtools sort -@ 20 "${base}.bam" -o "${base}_sorted.bam"
    
    # Index sorted BAM file (indexing is single-threaded)
    samtools index "${base}_sorted.bam"
    
    # Convert sorted BAM to BED (bedtools does not support multithreading)
    bedtools bamtobed -i "${base}_sorted.bam" > "${base}.bed"
done
