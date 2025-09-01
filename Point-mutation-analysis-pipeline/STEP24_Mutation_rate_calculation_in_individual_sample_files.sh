#!/bin/bash

# Input directory containing all .tsv.gz files
input_dir="Path_to_input_file_directory/Muatation_rate_raw_files"

# List all .tsv.gz files in the input directory
files=($(find "$input_dir" -name "*mutation_profile.tsv.gz"))

# Loop through each file and process them one by one
for file_path in "${files[@]}"; do
  # Extract the sample name (the base name of the file without the extension)
  sample_name=$(basename "$file_path" | sed -r 's/_mutation_profile.tsv.gz$//')

  # Output file for this specific sample
  output_file="${input_dir}/${sample_name}_mutation_summary.tsv"

  # Initialize counters
  total_entries=0
  no_mutations=0
  transversions=0
  transitions=0
  others=0

  # Check if the output file already exists; if not, add header
  if [ ! -f "$output_file" ]; then
    echo -e "Samples\tTotal_entries\tNo_mutations\tTransversion\tTransition\tOthers\tNo mutation %\tTransversion %\tTransition %\tOthers %" > "$output_file"
  fi

  # Process the file with zcat and awk
  zcat "$file_path" | awk '
    BEGIN {
      FS = "\t"
      # Skip the header row
      getline header
    }
    NR > 1 {
      total_entries++
      mutation_type = $11  # Mutation type is in the 11th column (Mutation_type_Ac)

      # Count mutations based on type
      if (mutation_type == "No mutation") no_mutations++
      else if (mutation_type == "Transversion") transversions++
      else if (mutation_type == "Transition") transitions++
      else if (mutation_type == "Others") others++
    }
    END {
      # Avoid division by zero if no valid entries
      if (total_entries > 0) {
        no_mutations_pct = (no_mutations / total_entries) * 100
        transversions_pct = (transversions / total_entries) * 100
        transitions_pct = (transitions / total_entries) * 100
        others_pct = (others / total_entries) * 100

        # Print results for this file
        printf "%s\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n", 
          "'$sample_name'", total_entries, no_mutations, transversions, transitions, others, 
          no_mutations_pct, transversions_pct, transitions_pct, others_pct
      }
    }
  ' >> "$output_file"

  # Output to indicate completion of this file processing
  echo "Processed: $file_path -> $output_file"
done
