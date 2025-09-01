#!/bin/bash

# Set the input and output directories
INPUT_DIR="Path_to_input_file_directory/Base_percentage_converted"
OUTPUT_DIR="Path_to_output_file_directory/Base_calls_to_unique_strings/"

# Check if both input and output directories are provided
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Both input and output directories must be set."
    exit 1
fi

# Function to process a single file
process_file() {
    INPUT_FILE=$1
    OUTPUT_DIR=$2

    # Extract the file name from the input file path
    FILENAME=$(basename "$INPUT_FILE")
    OUTPUT_FILE="$OUTPUT_DIR${FILENAME%.gz}_processed.gz"

    # Prepare the header
    HEADER="Chromosome\tPosition\tCoverage\tBase_called\tA\tC\tG\tT\tN"

    # Process the input file and output to a gzipped file, skipping the header line
    {
        # Print the header
        echo -e "$HEADER"
        
        # Read and process the file (skip the header line)
        zcat "$INPUT_FILE" | awk 'BEGIN {OFS="\t"} 
            NR == 1 {next}  # Skip the header line (first line)
            {
                # Get the base calls from the 4th column
                base_calls = $4

                # Create a set of unique nucleotides
                n = split(base_calls, arr, "")
                delete unique
                for (i = 1; i <= n; i++) {
                    unique[arr[i]] = 1
                }

                # Count the number of unique bases
                unique_count = 0
                for (base in unique) {
                    unique_count++
                }

                # If more than 2 unique bases, skip the line
                if (unique_count > 2) {
                    next
                }

                # Prepare the unique nucleotides as a comma-separated string
                base_seq = ""
                for (base in unique) {
                    if (base_seq != "") base_seq = base_seq "," base
                    else base_seq = base
                }

                # Replace the 4th column with the new base sequence
                $4 = base_seq

                # Print the modified line
                print $0
            }'
    } | gzip > "$OUTPUT_FILE"

    # Notify user that processing is complete for this file
    echo "Processing complete for $INPUT_FILE. Output saved to $OUTPUT_FILE"
}

export -f process_file

# Loop through each .tsv.gz file in the input directory (no subdirectories) and process them in parallel
find "$INPUT_DIR" -maxdepth 1 -name "*.tsv.gz" | parallel -j 11 process_file {} "$OUTPUT_DIR"
