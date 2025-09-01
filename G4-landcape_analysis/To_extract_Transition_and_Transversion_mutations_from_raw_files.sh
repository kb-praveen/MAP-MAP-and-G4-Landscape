#!/bin/bash

input_dir="Path_to_raw_input_mutation_files/Mutation"
output_dir="Path_to_output_mutation_files/"

mkdir -p "$output_dir"

for file in "$input_dir"/*.tsv.gz; do
    filename=$(basename "$file" .tsv.gz)

    # Transition
    zcat "$file" | awk -F'\t' '$11 == "Transition" {print $1, $2}' OFS='\t' | sort -k1,1 -k2,2n > "$output_dir/${filename}_Transition.bed"

    # Transversion
    zcat "$file" | awk -F'\t' '$11 == "Transversion" {print $1, $2}' OFS='\t' | sort -k1,1 -k2,2n > "$output_dir/${filename}_Transversion.bed"

    # Create shuffled +1 BED files (1 million entries)
    for type in Transition Transversion; do
        input_bed="$output_dir/${filename}_${type}.bed"
        sampled_bed="$output_dir/${filename}_${type}_1M_raw.bed"
        final_bed="$output_dir/${filename}_${type}_1M_final.bed"

        # Shuffle and pick 1 million lines, sort
        shuf -n 1000000 "$input_bed" | sort -k1,1 -k2,2n > "$sampled_bed"

        # Add +1 to column 2 to make proper BED format
        awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $2+1}' "$sampled_bed" > "$final_bed"
    done
done
