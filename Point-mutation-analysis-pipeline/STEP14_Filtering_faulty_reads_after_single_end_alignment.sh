#!/bin/bash

# File paths
read_ids_file="Path_to_input_file/A_common_reads_with_at_only_ids_parallel.txt"
bam_file="Path_to_input_file/A_R1_R2_sorted.bam"
output_file="Path_to_output_file/A_sample_faulty_reads_alignment.txt"

# Remove "@" from read IDs in the text file and store them in a new file
sed 's/^@//' "$read_ids_file" > /tmp/cleaned_read_ids.txt

# Extract the IDs from the BAM file and filter out the ones that exist in the cleaned read IDs file
samtools view "$bam_file" | \
awk 'BEGIN { while ((getline line < "/tmp/cleaned_read_ids.txt") > 0) read_ids[line]=1 } \
     !($1 in read_ids) { print $1"\t"$3"\t"$4"\t"$2 }' > "$output_file"

# Clean up
rm /tmp/cleaned_read_ids.txt

echo "Faulty reads information has been saved to: $output_file"
