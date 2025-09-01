#!/bin/bash

# Input directory containing BED files
input_dir="Path_to_input_directory/Human_G4_Analysis"

# Output directory (can be the same as input or different)
output_dir="$input_dir"

# Loop through all BED files in the directory
for input_file in "$input_dir"/*.bed; do
    filename=$(basename "$input_file" .bed)
    output_file="${output_dir}/${filename}_single_base.bed"

    # Empty the output file if it exists
    > "$output_file"

    # Read the input BED file line by line
    while IFS=$'\t' read -r chrom start end; do
        # Loop from start to end-1 and write single base entries
        for (( pos=start; pos<end; pos++ )); do
            echo -e "${chrom}\t${pos}\t$((pos + 1))" >> "$output_file"
        done
    done < "$input_file"

    echo "Processed: $input_file → $output_file"
done
