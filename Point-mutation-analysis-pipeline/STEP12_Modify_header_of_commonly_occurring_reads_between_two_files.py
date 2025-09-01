import os
import multiprocessing

# Function to prepend "@" to read IDs in a chunk of the text file
def process_chunk(chunk, output_file):
    with open(output_file, 'a') as outfile:
        for line in chunk:
            # Split the line to get the first column (read ID)
            read_id = line.split()[0]
            # Prepend '@' to the read ID and write it to the output file
            outfile.write(f"@{read_id}\n")

# Function to divide the input file into chunks and parallel process them
def parallel_process_text_file(input_file, output_file, num_workers):
    # Read the input file and divide it into chunks for parallel processing
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    # Calculate the chunk size for each worker
    chunk_size = len(lines) // num_workers
    chunks = [lines[i:i + chunk_size] for i in range(0, len(lines), chunk_size)]
    
    # Use multiprocessing to process each chunk in parallel
    with multiprocessing.Pool(processes=num_workers) as pool:
        # For each chunk, we append the results to the output file
        pool.starmap(process_chunk, [(chunk, output_file) for chunk in chunks])

# File paths
text_file = "Path_to_the_input_file/J_common_reads.txt"  # Original text file containing read IDs
output_file = "Path_to_the_output_file/J_common_reads_with_at_only_ids_parallel.txt"  # Output file with modified read IDs

# Number of CPU cores to use for parallel processing
num_workers = 20

# Remove the output file if it already exists to avoid appending to old data
if os.path.exists(output_file):
    os.remove(output_file)

# Run the parallel processing function
parallel_process_text_file(text_file, output_file, num_workers)

print(f"Modified read IDs (with '@') saved to {output_file}")
