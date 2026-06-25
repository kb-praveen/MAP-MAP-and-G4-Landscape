#!/bin/bash

# Fixed path to the input directory
INPUT_DIR="Path_to_the_bed_files/bed_files/Raw_bed"

for bedfile in "$INPUT_DIR"/*.bed; do
  output_file="${bedfile%.bed}_G4.bed"

  awk -F'\t' '
  {
    split($1, a, "[:-]");
    chr = a[1];         # e.g., chr1
    anchor = a[2] + 0;  # e.g., 9107
    start = anchor + $2;  # offset from column 2
    end = anchor + $3;    # offset from column 3
    print chr "\t" start "\t" end;
  }' "$bedfile" | sort --parallel=22 -k1,1 -k2,2n > "$output_file"

  echo "✅ Processed: $(basename "$bedfile") → $(basename "$output_file")"
done
