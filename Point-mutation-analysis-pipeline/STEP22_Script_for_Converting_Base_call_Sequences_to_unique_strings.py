import csv
import gzip
import os
from multiprocessing import Pool

# Input directory and output directory
input_dir = 'Path_to_input_file_directory/Base_calls_to_unique_strings'
output_dir = 'Path_to_output_file_directory/Base_calls_to_unique_strings'

# Function to rearrange bases based on frequencies
def rearrange_bases(base_called, freq_A, freq_C, freq_G, freq_T, freq_N):
    if ',' in base_called:  # Check if there are two bases
        # Split base_called into individual bases
        bases = base_called.split(',')
        
        # Create a dictionary of frequencies
        base_freqs = {
            'A': freq_A,
            'C': freq_C,
            'G': freq_G,
            'T': freq_T,
            'N': freq_N
        }
        
        # Sort the bases by frequency in descending order
        base_freqs_sorted = sorted(
            [(base, base_freqs[base]) for base in bases],
            key=lambda x: x[1],
            reverse=True
        )
        
        # Return the bases in the order of their frequencies
        return ','.join([base_freqs_sorted[0][0], base_freqs_sorted[1][0]])
    
    return base_called  # No change if only one base

# Function to process each file
def process_file(input_file):
    output_file = os.path.join(output_dir, os.path.basename(input_file).replace('_final_percentage.tsv_processed.gz', '_modified.tsv.gz'))
    
    # Open input file for reading and output file for writing
    with gzip.open(input_file, 'rt') as infile, gzip.open(output_file, 'wt', newline='') as outfile:
        # Create CSV reader and writer
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        # Read the header
        header = next(reader)
        writer.writerow(header)  # Write the header to the output file
        
        # Process each line in the input file
        for row in reader:
            # Extract the necessary columns (Base_called and frequencies)
            chromosome, position, coverage, base_called, freq_A, freq_C, freq_G, freq_T, freq_N = row
            
            # Convert frequencies to floats
            freq_A = float(freq_A)
            freq_C = float(freq_C)
            freq_G = float(freq_G)
            freq_T = float(freq_T)
            freq_N = float(freq_N)
            
            # Rearrange the base called based on frequencies
            new_base_called = rearrange_bases(base_called, freq_A, freq_C, freq_G, freq_T, freq_N)
            
            # Update the row with the new Base_called value
            row[3] = new_base_called  # Base_called is the 4th column (index 3)
            
            # Write the updated row to the output file
            writer.writerow(row)
    
    print(f"Modified file saved to {output_file}")

# Function to get all input files
def get_input_files(directory):
    return [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('_final_percentage.tsv_processed.gz')]

# Main function to handle parallel processing
def main():
    # Get all input files
    input_files = get_input_files(input_dir)

    # Create a Pool of 11 processes
    with Pool(processes=11) as pool:
        # Process files in parallel
        pool.map(process_file, input_files)

if __name__ == '__main__':
    main()
