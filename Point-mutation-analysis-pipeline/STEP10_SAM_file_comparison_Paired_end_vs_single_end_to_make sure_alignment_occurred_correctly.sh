#!/bin/bash

# Define the directory for temporary files and output
TEMP_DIR="Path_to_input_directory/Everything_on_plus_strand"

# Define the input SAM files
SAM_FILE_1="PATH_to_the_paired_end_alignment_sam_file/A_all_reads_on_plus_strand.sam"
SAM_FILE_2="PATH_to_the_single_end_alignment_sam_file/A_R1_R2.sam"

# Define the output SAM file
OUTPUT_FILE="${TEMP_DIR}/A_matching_reads.sam"
MISALIGNED_FILE="${TEMP_DIR}/misaligned_reads.txt"

# Temporary files to store extracted information from both SAM files
TMP_FILE1="${TEMP_DIR}/file1_tmp.txt"
TMP_FILE2="${TEMP_DIR}/file2_tmp.txt"

# Matching read IDs and coordinates
MATCHING_READS="${TEMP_DIR}/matching_reads.txt"
MATCHING_READ_IDS="${TEMP_DIR}/matching_read_ids.txt"
MISALIGNED_READS="${TEMP_DIR}/misaligned_reads.txt"

# Step 1: Extract read ID, chromosome, and position from both SAM files
# Exclude header lines (lines starting with @), then extract the first column (read ID), third column (chromosome), and fourth column (position)
awk '!/^@/ {print $1"\t"$3"\t"$4}' $SAM_FILE_1 > $TMP_FILE1
awk '!/^@/ {print $1"\t"$3"\t"$4}' $SAM_FILE_2 > $TMP_FILE2

# Step 2: Use `join` to find matching read IDs and coordinates (only keep matching read ID and chromosome+position)
# Sort the temporary files first
sort $TMP_FILE1 > ${TMP_FILE1}.sorted
sort $TMP_FILE2 > ${TMP_FILE2}.sorted

# Join the files on read ID and coordinates, output matching entries
join -1 1 -2 1 -o 1.1,1.2,1.3 ${TMP_FILE1}.sorted ${TMP_FILE2}.sorted > $MATCHING_READS

# Step 3: Extract the matching read IDs and coordinates from the joined file
# This will give us the read IDs that exist in both files with matching coordinates
awk '{print $1"\t"$2"\t"$3}' $MATCHING_READS > $MATCHING_READ_IDS

# Step 4: Extract matching reads from the original SAM file (file1.sam) and write them to the new SAM file
awk '/^@/ {print}' $SAM_FILE_1 > $OUTPUT_FILE  # Copy header lines

# Extract reads from SAM_FILE_1 with matching read IDs and coordinates
awk 'NR==FNR {match[$1"\t"$2"\t"$3]; next} ($1"\t"$3"\t"$4) in match' $MATCHING_READ_IDS $SAM_FILE_1 >> $OUTPUT_FILE

# Step 5: Extract misaligned reads from SAM_FILE_1 and write them to the misaligned file
# First, copy the header lines from SAM_FILE_1 to misaligned reads file
awk '/^@/ {print}' $SAM_FILE_1 > $MISALIGNED_FILE

# Extract non-matching reads from SAM_FILE_1 (reads that do not have matching coordinates or read IDs)
awk 'NR==FNR {match[$1"\t"$2"\t"$3]; next} !($1"\t"$3"\t"$4) in match' $MATCHING_READ_IDS $SAM_FILE_1 >> $MISALIGNED_FILE

# Step 6: Keep the temporary matching reads file
echo "Matching reads saved in $MATCHING_READS"
echo "Misaligned reads saved in $MISALIGNED_FILE"

# Step 7: Clean up temporary files
rm $TMP_FILE1 $TMP_FILE2 ${TMP_FILE1}.sorted ${TMP_FILE2}.sorted

echo "Matching reads have been saved to $OUTPUT_FILE"
echo "Misaligned reads have been saved to $MISALIGNED_FILE"
