#!/bin/bash

# Define input/output directories and reference genome
INPUT_DIR="Path_to_input_files/Mutations_segregated/Transition"
OUTPUT_DIR="Path_to_output_files/Mutations_segregated/Transition/1kb_up_down"
GENOME="Path_to_genome_fasta/hg38.fa"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Process each .bed file
for bed in "$INPUT_DIR"/*.bed; do
    base=$(basename "$bed" .bed)
    
    # 1. Expand coordinates by 1000bp up and down (i.e., start - 1000, end + 1000)
    awk 'BEGIN{OFS="\t"} {start=$2-1000; if(start<0) start=0; end=$3+1000; print $1, start, end}' "$bed" > "$OUTPUT_DIR/${base}_1kb_up_down.bed"
    
    # 2. Extract FASTA sequences
    bedtools getfasta -fi "$GENOME" -bed "$OUTPUT_DIR/${base}_1kb_up_down.bed" -fo "$OUTPUT_DIR/${base}_1kb_up_down.fa"
done
