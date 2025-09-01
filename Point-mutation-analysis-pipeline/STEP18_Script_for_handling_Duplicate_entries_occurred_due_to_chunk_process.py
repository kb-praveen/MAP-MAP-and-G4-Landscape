import os
import pandas as pd
from collections import defaultdict
import gzip

# Define input and output directories
input_dir = 'Path_to_the_input_file_directory/Raw_Files'
output_dir = 'Path_to_the_output_file_directory/Error_correction_raw_file_generation'

# Step 1: Get all .tsv.gz files from the input directory
input_files = [f for f in os.listdir(input_dir) if f.endswith('.tsv.gz')]

# Step 2: Process each file
for file in input_files:
    file_path = os.path.join(input_dir, file)
    print(f"Processing {file_path}")

    # Step 2.1: Load the data into a pandas dataframe
    df = pd.read_csv(file_path, sep='\t', header=None, compression='gzip', names=['chr', 'position', 'count', 'sequence', 'samples'])

    # Step 2.2: Determine the corresponding duplicated regions file
    duplicated_file = os.path.join(input_dir, 'Duplicate_entries', file.replace('.tsv.gz', '_duplicated_regions.txt'))
    print(f"Loading duplicated regions from: {duplicated_file}")
    
    # Read duplicated positions into a set for fast lookup
    duplicated_positions = set()
    with open(duplicated_file, 'r') as f:
        for line in f:
            chrom, pos = line.strip().split()
            duplicated_positions.add((chrom, pos))

    # Step 2.3: Create a dictionary to store merged results for duplicates
    merged_data = defaultdict(lambda: {'sequence': '', 'samples': [], 'count': 0})

    # Step 2.4: Process all rows to either merge duplicates or print unique ones
    non_duplicated_df = []
    for idx, row in df.iterrows():
        chrom, pos, count, sequence, samples = row
        if (chrom, pos) in duplicated_positions:
            # Merge duplicate rows
            merged_data[(chrom, pos)]['sequence'] += sequence  # Merge sequences
            merged_data[(chrom, pos)]['samples'].append(samples)  # Merge sample IDs
            merged_data[(chrom, pos)]['count'] += count  # Sum counts
        else:
            # Keep non-duplicate rows
            non_duplicated_df.append(row)

    # Step 2.5: Prepare the final data by combining non-duplicates and merged duplicates
    output_data = non_duplicated_df

    for (chrom, pos), data in merged_data.items():
        merged_line = [chrom, pos, data['count'], data['sequence'], ', '.join(data['samples'])]
        output_data.append(merged_line)

    # Step 2.6: Write the merged data to a new gzipped file
    output_file = os.path.join(output_dir, f"merged_{file}")
    
    with gzip.open(output_file, 'wt') as out_file:
        # Write the header
        out_file.write("chr\tposition\tcount\tsequence\tsamples\n")
        
        # Write the data
        for row in output_data:
            out_file.write("\t".join(map(str, row)) + "\n")
            print("\t".join(map(str, row)))  # Print to screen

    print(f'Merged file saved to {output_file}')

print("All files processed successfully.")
