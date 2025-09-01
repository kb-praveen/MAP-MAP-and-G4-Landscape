#!/bin/bash

# Input and output directories (adjust these paths if needed)
input_dir="Path_to_input_files/The_final_raw_output_file"
output_dir="Path_to_output_files/Coverage_filtered_files_4_and_above"

# Loop through all .tsv.gz files in the input directory
for input_file in "$input_dir"/*.tsv.gz; do
    # Get the filename without the path
    filename=$(basename "$input_file")
    
    # Define the output file path
    output_file="$output_dir/filtered_${filename}"
    
    # Use awk to filter rows where Coverage is between 4 and 1000 (inclusive)
    # Assuming Coverage is in the 3rd column (adjust if needed)
    zcat "$input_file" | awk -F'\t' '{ if ($3 >= 4 && $3 <= 1000) print $0 }' | gzip > "$output_file"
    
    echo "Filtered file saved as: $output_file"
done
