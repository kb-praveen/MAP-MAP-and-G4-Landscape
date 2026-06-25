#!/bin/bash

input_dir="/media/umash/Extreme SSD/Vertebrates_genome_1kb_up_all_promoters/Human_G4_mutaion_analysis_in_all_1kb_up_promoters/raw_mutation_tsv_files"
output_dir="/media/umash/Extreme SSD/Vertebrates_genome_1kb_up_all_promoters/Human_G4_mutaion_analysis_in_all_1kb_up_promoters/Mutation_segragated"

mkdir -p "$output_dir"

for file in "$input_dir"/*.tsv; do
    filename=$(basename "$file" .tsv)

    # Transition: print chrom, start, start+1 with tab separation
    awk -F'\t' 'BEGIN {OFS="\t"} $11 == "Transition" {print $1, $2, $2+1}' "$file" | sort -k1,1 -k2,2n > "$output_dir/${filename}_Transition.bed"

    # Transversion: print chrom, start, start+1 with tab separation
    awk -F'\t' 'BEGIN {OFS="\t"} $11 == "Transversion" {print $1, $2, $2+1}' "$file" | sort -k1,1 -k2,2n > "$output_dir/${filename}_Transversion.bed"

    # C -> T mutations: extract rows where column 10 contains "C -> T", print chrom, start, start+1 with tab separation
    awk -F'\t' 'BEGIN {OFS="\t"} $10 ~ /C -> T/ {print $1, $2, $2+1}' "$file" | sort -k1,1 -k2,2n > "$output_dir/${filename}_C_to_T.bed"
done

echo "✅ BED files generated in: $output_dir"
