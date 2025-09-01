import gzip
import re
import os
import time
from collections import defaultdict

# Function to parse CIGAR string and adjust base coordinates
def parse_cigar(cigar, start_position, seq):
    """Parse CIGAR string and return the adjusted positions and bases"""
    position = start_position
    bases = []
    
    # CIGAR parsing
    i = 0
    seq_index = 0  # To keep track of the current position in the sequence
    while i < len(cigar):
        # Match a sequence of digits (the length) and the operation
        match = re.match(r'(\d+)([MDIX=])', cigar[i:])
        if match:
            length = int(match.group(1))  # Length is the first captured group (digits)
            op = match.group(2)  # Operation is the second captured group (letter)
            
            if op == 'M' or op == '=' or op == 'X':  # Match
                for _ in range(length):
                    base = seq[seq_index]  # Get the base from the sequence
                    bases.append((base, position))
                    position += 1
                    seq_index += 1
            elif op == 'I':  # Insertion, skip
                seq_index += length  # Skip the inserted bases
            elif op == 'D':  # Deletion, just skip position
                position += length
            
            # Move past this match in the string
            i += len(match.group(0))
        else:
            # If no match is found, skip a character (this shouldn't happen in well-formed CIGAR strings)
            i += 1
    
    return bases

# Function to process the input file in chunks
def process_reads_in_chunks(input_file, output_file, chunk_size=10000):
    read_data = defaultdict(lambda: defaultdict(list))  # stores (chrom, pos) -> (list of base, read ids)
    
    with gzip.open(input_file, 'rt') as f, gzip.open(output_file, 'wt') as out:
        # Write headers to output file
        out.write("Chromosome\tPosition\tCoverage\tBase_called\tRead_IDs\n")
        
        batch = []
        line_count = 0
        total_lines_processed = 0
        start_time = time.time()  # Start time for the entire processing
        
        for line in f:
            fields = line.strip().split('\t')
            chrom, start, end, read_id, _, strand, cigar, seq = fields[:8]
            start, end = int(start), int(end)
            
            # Process the CIGAR string to get base positions and actual base from seq
            base_positions = parse_cigar(cigar, start, seq)
            
            # Fill the read data structure for the current line
            for base, pos in base_positions:
                read_data[chrom][pos].append((base, read_id))
            
            batch.append(line)  # Store the line in the batch
            line_count += 1
            total_lines_processed += 1
            
            # When we reach the chunk size, process the batch
            if line_count >= chunk_size:
                # Process the batch and write to the output file
                chunk_start_time = time.time()  # Start time for the current chunk
                
                for chrom in sorted(read_data):
                    for pos in sorted(read_data[chrom]):
                        bases = ''.join([base for base, _ in read_data[chrom][pos]])  # Join bases
                        read_ids = ', '.join(sorted(set(read_id for _, read_id in read_data[chrom][pos])))  # Unique read IDs
                        coverage = len(read_data[chrom][pos])  # Coverage is the number of occurrences at the position
                        out.write(f"{chrom}\t{pos}\t{coverage}\t{bases}\t{read_ids}\n")
                
                # Print or log how many lines were processed and the time taken
                chunk_elapsed_time = time.time() - chunk_start_time
                print(f"Processed {total_lines_processed} lines in this batch. Time taken: {chunk_elapsed_time:.2f} seconds.")
                
                # Reset data for the next batch
                read_data.clear()
                batch = []
                line_count = 0
        
        # Process any remaining lines after the last batch
        if batch:
            for chrom in sorted(read_data):
                for pos in sorted(read_data[chrom]):
                    bases = ''.join([base for base, _ in read_data[chrom][pos]])  # Join bases
                    read_ids = ', '.join(sorted(set(read_id for _, read_id in read_data[chrom][pos])))  # Unique read IDs
                    coverage = len(read_data[chrom][pos])  # Coverage is the number of occurrences at the position
                    out.write(f"{chrom}\t{pos}\t{coverage}\t{bases}\t{read_ids}\n")
        
        # Print or log the total processing time
        total_elapsed_time = time.time() - start_time
        print(f"Total processing time: {total_elapsed_time:.2f} seconds.")
        
# Example usage
input_file = 'Path_to_the_input_file/A_output_with_sequences.gz'  # Path to your zipped input file
output_file = 'Path_to_the_output_file/A_file_final.tsv.gz'    # Desired output path for zipped output

process_reads_in_chunks(input_file, output_file)
