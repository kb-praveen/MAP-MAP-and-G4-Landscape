#!/bin/bash

# Loop over all BED files ending with "_sorted.bed"
for bedfile in *_sorted.bed; do
    # Extract base name (without extension)
    base=$(basename "$bedfile" .bed)

    echo "Processing: $bedfile"

    # Run awk to extract columns 1 and 2, save to new file
    awk 'BEGIN {OFS="\t"} {print $1, $2}' "$bedfile" > "${base}_grep.bed"
done
