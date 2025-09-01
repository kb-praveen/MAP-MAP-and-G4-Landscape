import gzip
import concurrent.futures

# Load the read IDs from the text file into a set for fast lookup
def load_read_ids(text_file):
    read_ids = set()
    with open(text_file, 'r') as f:
        for line in f:
            read_id = line.strip()  # Read ID is the entire line
            read_ids.add(read_id)
    print(f"Loaded {len(read_ids)} read IDs.")
    return read_ids

# Process a block of reads and filter the ones with matching IDs
def process_block(fastq_lines, read_ids):
    filtered_reads = []
    for i in range(0, len(fastq_lines), 4):
        header = fastq_lines[i].strip()
        sequence = fastq_lines[i + 1].strip()
        plus = fastq_lines[i + 2].strip()
        quality = fastq_lines[i + 3].strip()

        # Extract the read ID from the header (the part before the first space)
        read_id = header.split()[0]

        if read_id in read_ids:
            filtered_reads.append(f"{header}\n{sequence}\n{plus}\n{quality}\n")
    return filtered_reads

# Process the gzipped FASTQ file and filter the reads
def filter_fastq(fastq_gz_file, read_ids, output_gz_file, num_workers=20):
    with gzip.open(fastq_gz_file, 'rt') as fastq, gzip.open(output_gz_file, 'wt') as out:
        # Read the entire FASTQ file into memory in chunks
        fastq_lines = []
        count_filtered = 0

        # Read in blocks of 1000 lines (250 reads)
        for line in fastq:
            fastq_lines.append(line)
            if len(fastq_lines) == 4000:  # 1000 reads = 4000 lines (4 lines per read)
                # Process the block in parallel
                filtered_reads = process_block(fastq_lines, read_ids)
                out.writelines(filtered_reads)
                count_filtered += len(filtered_reads) // 4
                fastq_lines = []  # Clear the buffer

        # Process any remaining lines
        if fastq_lines:
            filtered_reads = process_block(fastq_lines, read_ids)
            out.writelines(filtered_reads)
            count_filtered += len(filtered_reads) // 4

        print(f"Total reads matching: {count_filtered}")
        print(f"Filtered gzipped FASTQ file saved as {output_gz_file}")

# File paths
text_file = "Path_to_the_text_file/K_common_reads_with_at_only_ids_parallel.txt"  # Your text file containing read IDs
fastq_gz_file = "Path_to_the_fastq_file/K_all_reads_on_plus_strand.fastq.gz"  # Your gzipped input FASTQ file
output_gz_file = "Path_to_the_output_file/K_all_reads_on_plus_strand_filtered.fastq.gz"  # The output gzipped FASTQ file with matching reads

# Load the read IDs
read_ids = load_read_ids(text_file)

# Filter the gzipped FASTQ file
filter_fastq(fastq_gz_file, read_ids, output_gz_file)
