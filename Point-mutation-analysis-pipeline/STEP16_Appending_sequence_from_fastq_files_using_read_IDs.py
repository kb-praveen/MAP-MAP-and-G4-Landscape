import gzip
import os
import concurrent.futures

def read_fastq(file_path):
    """Reads the fastq file and returns a dictionary of read ids and their sequences."""
    fastq_dict = {}
    with gzip.open(file_path, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # Skip plus sign (+)
            f.readline()  # Skip quality line
            read_id = header.split()[0][1:]  # Removing the '@' symbol
            fastq_dict[read_id] = seq
    return fastq_dict

def process_bed_chunk(bed_chunk, fastq_dict):
    """Process a chunk of the bed file and return the lines with sequences appended."""
    result = []
    for line in bed_chunk:
        columns = line.strip().split("\t")
        read_id = columns[3]
        if read_id in fastq_dict:
            sequence = fastq_dict[read_id]
            result.append("\t".join(columns) + "\t" + sequence + "\n")
        else:
            result.append("\t".join(columns) + "\t" + "NOT_FOUND" + "\n")
    return result

def split_bed_file(bed_file_path, chunk_size):
    """Splits the BED file into chunks."""
    chunks = []
    with open(bed_file_path, 'r') as bed_file:
        chunk = []
        for line in bed_file:
            chunk.append(line)
            if len(chunk) >= chunk_size:
                chunks.append(chunk)
                chunk = []
        if chunk:
            chunks.append(chunk)
    return chunks

def process_bed_file_parallel(bed_file_path, fastq_dict, output_file_path, num_workers):
    """Process the bed file in parallel and save the output to a gzip file."""
    # Split the BED file into chunks
    bed_chunks = split_bed_file(bed_file_path, chunk_size=1000000)  # Adjust chunk size as needed
    
    # Open the gzip output file
    with gzip.open(output_file_path, 'wt') as output_file:
        with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
            # Submit the tasks to the executor
            futures = [executor.submit(process_bed_chunk, chunk, fastq_dict) for chunk in bed_chunks]
            
            # Collect and write the results
            for future in concurrent.futures.as_completed(futures):
                for result in future.result():
                    output_file.write(result)

# File paths
bed_file_path = 'Path_to_bed_file/A_all_reads_on_plus_strand_filtered_sorted_with_cigar.bed'
fastq_file_path = 'Path_to_fastq_file/A_all_reads_on_plus_strand_filtered.fastq.gz'

# Define the output file path (compressed .gz)
output_file_path = 'Path_to_output_file/A_output_with_sequences.gz'

# Number of worker threads (cores available)
num_workers = 22

# Ensure the output directory exists
os.makedirs(os.path.dirname(output_file_path), exist_ok=True)

# Step 1: Read the fastq file into a dictionary
fastq_dict = read_fastq(fastq_file_path)

# Step 2: Process the bed file in parallel and write the output with sequences (to a gzip-compressed file)
process_bed_file_parallel(bed_file_path, fastq_dict, output_file_path, num_workers)

print(f"Processing complete. Output saved to '{output_file_path}'.")
