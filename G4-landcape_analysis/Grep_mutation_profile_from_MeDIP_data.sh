#!/bin/bash

# Define file paths
bed_file="Path_of_bed_file/Homo_sapiens_G4_sorted_single_base_grep.bed"
mutation_file="Path_to_mutation_calculated_file/Hs_mutation_profile.tsv.gz"
output_file="Path_to_output_file/Hs_genomewide_G4_mutation_profile.tsv"

# Build associative array from bed file
awk '{
    key = $1 ":" $2
    pos[key] = 1
} END {
    for (k in pos) print k
}' "$bed_file" > keys.txt

# Use the keys to filter from mutation file
zcat "$mutation_file" | awk -F'\t' -v keys_file="keys.txt" '
BEGIN {
    while ((getline line < keys_file) > 0) {
        key[line] = 1
    }
}
NR==1 {
    print; next
}
{
    k = $1 ":" $2
    if (k in key)
        print
}
' > "$output_file"

# Cleanup
rm keys.txt

echo "Output written to $output_file"
