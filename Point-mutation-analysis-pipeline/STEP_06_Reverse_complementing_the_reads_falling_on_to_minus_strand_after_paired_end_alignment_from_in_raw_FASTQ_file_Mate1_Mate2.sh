#!/bin/bash

# Set the directory where the .fastq.gz files are located
input_dir="Path_to_input_file_directory/Minus_strand_fastq"
output_dir="Path_to_output_file_directory/Minus_strand_fastq/Reverse_complemented"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through all .fastq.gz files in the input directory
for file in "$input_dir"/*.fastq.gz; do
    # Extract the filename without the path and extension
    filename=$(basename "$file" .fastq.gz)
    
    # Reverse and complement the sequence using seqkit (both reverse and complement)
    seqkit seq -rp "$file" -o "$output_dir/${filename}_reverse_complemented.fastq.gz"
    
    echo "Processed $file and saved reverse complemented file as $output_dir/${filename}_reverse_complemented.fastq.gz"
done
