#!/bin/bash

# Function to process each file
process_file() {
    local input_file=$1
    local output_directory=$2
    local output_file="${output_directory}/$(basename ${input_file%.tsv.gz}_processed.tsv.gz)"
    
    echo "Processing $input_file"
    
    # Process the file and create a new one with additional columns
    zcat "$input_file" | awk -F'\t' 'BEGIN {OFS="\t"}
    NR==1 {print $0, "Base_conversion", "Mutation"}  # Add headers for new columns
    NR>1 {
        split($3, bases, ",")
        
        # Skip rows where more than 2 bases are reported in Base_called column
        if (length(bases) > 2) {
            next
        }
        
        # For single base reported, classify Base_conversion and Mutation as NA
        if (length(bases) == 1) {
            base_conversion = "NA"
            mutation = "NA"
            print $0, base_conversion, mutation
            next
        }
        
        # Determine the base percentages from columns A, C, G, T
        max_percent = 0
        min_percent = 100
        max_base = ""
        min_base = ""
        
        for (i=4; i<=7; i++) {
            if ($i > max_percent) {
                max_percent = $i
                if (i == 4) max_base = "A"
                if (i == 5) max_base = "C"
                if (i == 6) max_base = "G"
                if (i == 7) max_base = "T"
            }
            if ($i < min_percent && $i > 0) {
                min_percent = $i
                if (i == 4) min_base = "A"
                if (i == 5) min_base = "C"
                if (i == 6) min_base = "G"
                if (i == 7) min_base = "T"
            }
        }
        
        # Reorder Base_called based on the percentages
        if (max_percent > min_percent) {
            $3 = max_base "," min_base
        } else {
            $3 = min_base "," max_base
        }
        
        # Determine the base conversion direction (A->G, G->A, C->T, T->C)
        base_conversion = max_base " -> " min_base
        
        # Determine mutation type (Transition or Transversion)
        if ((max_base == "A" && min_base == "G") || (max_base == "G" && min_base == "A") || 
            (max_base == "C" && min_base == "T") || (max_base == "T" && min_base == "C")) {
            mutation = "Transition"
        } else {
            mutation = "Transversion"
        }
        
        print $0, base_conversion, mutation
    }' | gzip > "$output_file"
    
    echo "Processed and saved: $output_file"
}

# Main function to process all files
process_all_files() {
    local input_directory=$1
    local output_directory=$2
    
    # Create output directory if it doesn't exist
    mkdir -p "$output_directory"
    
    # Process each .tsv.gz file in the input directory
    for file in "$input_directory"/*.tsv.gz; do
        process_file "$file" "$output_directory" &
    done
    
    # Wait for all background processes to finish
    wait
}

# Path to the input directory (containing the .tsv.gz files)
input_directory="Path_to_input_file_diretory/Converted_to_string_to_uniq_base"

# Path to the output directory (where processed files will be saved)
output_directory="Path_to_output_file_diretory/Base_conversion"

# Run the processing function
process_all_files "$input_directory" "$output_directory"
