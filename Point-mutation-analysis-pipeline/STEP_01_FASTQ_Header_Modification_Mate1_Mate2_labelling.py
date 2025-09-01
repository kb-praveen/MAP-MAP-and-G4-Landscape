import gzip
import os
import concurrent.futures

def modify_fastq_header(input_file, output_file):
    """Modify the header by appending _1 to each header line."""
    with gzip.open(input_file, 'rt') as infile, gzip.open(output_file, 'wt') as outfile:
        while True:
            # Read the four lines corresponding to a single FASTQ entry
            header = infile.readline()
            if not header:
                break  # EOF reached
            sequence = infile.readline()
            plus_line = infile.readline()
            quality = infile.readline()

            # Modify the header: append "_1" before the first space (keeping everything else)
            if header.startswith('@'):
                # Find the position of the first space (end of the identifier part)
                space_pos = header.find(' ')
                if space_pos != -1:
                    header = header[:space_pos] + "_1" + header[space_pos:]  # Append "_1" before the rest of the metadata

            # Write the modified header and the rest of the lines to the output file
            outfile.write(header)
            outfile.write(sequence)
            outfile.write(plus_line)
            outfile.write(quality)

def process_single_file(input_file, output_file):
    """Process a single FASTQ file."""
    print(f"Processing: {input_file} -> {output_file}")
    modify_fastq_header(input_file, output_file)

def process_files(input_dir, output_dir):
    """Process all Mate 1 files in the input directory concurrently."""
    files_to_process = [
        os.path.join(input_dir, f) 
        for f in os.listdir(input_dir) 
        if f.endswith(".fastq.gz")
    ]
    
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Use ThreadPoolExecutor to process files concurrently, with 22 threads
    with concurrent.futures.ThreadPoolExecutor(max_workers=22) as executor:
        # For each file, start a new thread to process it
        for input_file in files_to_process:
            # Prepare output file name by appending _1
            output_file = os.path.join(output_dir, os.path.basename(input_file).replace(".fastq.gz", "_1.fastq.gz"))
            # Submit a task to process the file
            executor.submit(process_single_file, input_file, output_file)

if __name__ == "__main__":
    input_dir = "Path_to_input_files/Mate_1 or Mate_2"
    output_dir = "Path_to_output_files/Mate_1 or Mate_2"
    
    # Process Mate 1 or Mate 2 files with "_1" or "_2" suffix
    process_files(input_dir, output_dir)
