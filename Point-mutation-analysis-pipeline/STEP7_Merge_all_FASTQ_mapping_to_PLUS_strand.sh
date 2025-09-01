#!/bin/bash

# Set the directory where the sample FASTQ files are located
input_dir="Path_to_the_input_directory"
output_dir="$input_dir"  # Output will be in the same directory

# Number of threads (not needed when using cat)
threads=1

# Loop through the sample names from A to K
for sample in {A..K}; do
    # Define the input files for the current sample
    R1_minus="${input_dir}/${sample}_R1_minus_strand_reads_rc.fastq.gz"
    R1_plus="${input_dir}/${sample}_R1_plus_strand_reads.fastq.gz"
    R2_minus="${input_dir}/${sample}_R2_minus_strand_reads_rc.fastq.gz"
    R2_plus="${input_dir}/${sample}_R2_plus_strand_reads.fastq.gz"
    
    # Define the output file name (gzipped)
    output_file="${output_dir}/${sample}_all_reads_on_plus_strand.fastq.gz"
    
    # Use `cat` to concatenate the files (no threading needed here)
    cat "$R1_minus" "$R1_plus" "$R2_minus" "$R2_plus" > "$output_file"
    
    echo "Processed ${sample}, saved to $output_file"
    
    # Remove the input files after the final file has been generated
    rm "$R1_minus" "$R1_plus" "$R2_minus" "$R2_plus"
    echo "Removed input files for ${sample}"
done
