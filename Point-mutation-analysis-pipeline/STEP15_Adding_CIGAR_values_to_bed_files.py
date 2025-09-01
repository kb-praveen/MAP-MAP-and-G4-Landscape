import pysam
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

# Function to extract the CIGAR value for a given read
def extract_cigar(read):
    if not read.is_unmapped:
        return read.query_name, read.cigarstring
    return None

# Step 1: Open BAM file and initialize thread pool
bam_file = "Path_to_bam_file/A_all_reads_on_plus_strand_filtered_sorted.bam"
samfile = pysam.AlignmentFile(bam_file, "rb")  # Open BAM file

# Step 2: Use ThreadPoolExecutor to extract CIGAR values in parallel
read_cigar_dict = {}

with ThreadPoolExecutor(max_workers=22) as executor:  # 22 threads (or set as needed)
    futures = []
    for read in samfile.fetch():
        futures.append(executor.submit(extract_cigar, read))
    
    for future in futures:
        result = future.result()
        if result:
            read_cigar_dict[result[0]] = result[1]

samfile.close()

# Step 3: Read the existing BED file
bed_file = "Path_to_bed_file/A_all_reads_on_plus_strand_filtered_sorted.bed"
bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=["chrom", "start", "end", "read_name", "score", "strand"])

# Step 4: Add CIGAR values to the BED dataframe by matching the read IDs
bed_df["cigar"] = bed_df["read_name"].map(read_cigar_dict)

# Step 5: Save the modified BED file with the CIGAR column
output_bed_file = "Path_to_output_file/A_all_reads_on_plus_strand_filtered_sorted_with_cigar.bed"
bed_df.to_csv(output_bed_file, sep="\t", header=False, index=False)

print(f"Modified BED file with CIGAR column saved to {output_bed_file}")
