# Function to read the file and return a set of tuples for comparison
def read_file(filename):
    data = set()  # To store unique rows as tuples
    with open(filename, 'r') as f:
        for line in f:
            # Assuming columns are separated by spaces or tabs
            parts = line.strip().split()  
            if len(parts) == 3:
                read_id, chromosome, coordinate = parts
                data.add((read_id, chromosome, coordinate))
    return data

# Function to compare the two files and write the matching rows
def compare_files(file1, file2, output_file):
    # Read the data from both files
    data1 = read_file(file1)
    data2 = read_file(file2)
    
    # Find the intersection (common rows)
    common_data = data1.intersection(data2)
    
    # Write the common rows to the output file
    with open(output_file, 'w') as f:
        for row in common_data:
            f.write('\t'.join(row) + '\n')

# Input and output files
file1 = 'Path_to_input_file/Everything_on_plus_strand/File_1_sorted.txt'
file2 = 'Path_to_input_file/Everything_on_plus_strand/File_2_sorted.txt'
output_file = 'Path_to_output_file/Everything_on_plus_strand/A_common_reads.txt'

# Compare the files and write the matching rows to the output file
compare_files(file1, file2, output_file)

print(f"Matching rows written to {output_file}")
