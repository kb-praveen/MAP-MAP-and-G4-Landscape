#!/bin/bash

# Set the input directory and the output directory
input_dir="Path_to_input_directory/Coverage_4_to_1k_final"
output_dir="Path_to_output_directory/Base_percentage_converted"

# Ensure the output directory exists
mkdir -p "$output_dir"

# Function to process each file
process_file() {
  local file="$1"
  local output_dir="$2"
  
  # Decompress the file, process it, and output modified file to the output directory
  zcat "$file" | awk 'BEGIN {OFS="\t"} 
  NR==1 {print $0, "A", "C", "G", "T", "N"} 
  NR>1 {
    a=0; c=0; g=0; t=0; n=0; 
    # Iterate over each character in the 4th column (Base_calls)
    for(i=1; i<=length($4); i++) { 
      base=substr($4, i, 1); 
      if(base=="A") a++; 
      else if(base=="C") c++; 
      else if(base=="G") g++; 
      else if(base=="T") t++; 
      else if(base=="N") n++; 
    } 
    total=length($4); 
    # Print the original fields along with the calculated percentages for A, C, G, T, N
    print $1, $2, $3, $4, a/total*100, c/total*100, g/total*100, t/total*100, n/total*100
  }' | gzip > "$output_dir/$(basename "$file" .tsv.gz)_unsorted.tsv.gz"
}

# Export the function for parallel to use
export -f process_file

# Find all *.tsv.gz files in the input directory and process them in parallel with 6 cores
find "$input_dir" -name "*.tsv.gz" | parallel -j 11 process_file {} "$output_dir"
