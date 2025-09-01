#!/bin/bash

# Set directory for Fastp files
cd /PATH_to_the_directory:where files are present/ || { echo "Directory not found! Exiting..."; exit 1; }

# Function to process each sample
process_sample() {
    SAMPLE=$1

    # Set file paths for the current sample
    SAM="${SAMPLE}_R1_R2.sam"
    FASTQ_R1="${SAMPLE}_R1_fastp_1.fastq.gz"
    FASTQ_R2="${SAMPLE}_R2_fastp_2.fastq.gz"

    # Output file names with extensions based on the sample name
    PLUS_FILE_R1="${SAMPLE}_R1_plus_strand_reads.fastq"
    MINUS_FILE_R1="${SAMPLE}_R1_minus_strand_reads.fastq"
    PLUS_FILE_R2="${SAMPLE}_R2_plus_strand_reads.fastq"
    MINUS_FILE_R2="${SAMPLE}_R2_minus_strand_reads.fastq"

    # Create/empty the output files before appending
    > $PLUS_FILE_R1
    > $MINUS_FILE_R1
    > $PLUS_FILE_R2
    > $MINUS_FILE_R2

    echo "Processing SAM file: $SAM"

    # Step 1: Extract read IDs for plus (+) and minus (-) strand from SAM file, ignoring header lines, and remove duplicates
    echo "Extracting plus strand read IDs (flags 0, 99, 163)..."
    # Flag 0, 99: Plus strand, Mate 1
    # Flag 163: Plus strand, Mate 2 (properly paired)
    awk '!/^@/ {if ($2 == 0 || $2 == 99 || $2 == 163) print $1}' $SAM | sort | uniq > plus_strand_read_ids.txt

    echo "Extracting minus strand read IDs (flags 16, 147, 83)..."
    # Flag 16, 147: Minus strand, Mate 2
    # Flag 83: Minus strand, Mate 1 (properly paired)
    awk '!/^@/ {if ($2 == 16 || $2 == 147 || $2 == 83) print $1}' $SAM | sort | uniq > minus_strand_read_ids.txt

    # Step 2: Extract plus strand reads from FASTQ files using seqtk
    echo "Extracting plus strand reads from FASTQ files..."
    # For Mate 1 (_1) plus strand reads
    seqtk subseq $FASTQ_R1 plus_strand_read_ids.txt >> $PLUS_FILE_R1
    # For Mate 2 (_2) plus strand reads
    seqtk subseq $FASTQ_R2 plus_strand_read_ids.txt >> $PLUS_FILE_R2

    # Step 3: Extract minus strand reads from FASTQ files using seqtk
    echo "Extracting minus strand reads from FASTQ files..."
    # For Mate 1 (_1) minus strand reads
    seqtk subseq $FASTQ_R1 minus_strand_read_ids.txt >> $MINUS_FILE_R1
    # For Mate 2 (_2) minus strand reads
    seqtk subseq $FASTQ_R2 minus_strand_read_ids.txt >> $MINUS_FILE_R2

    # Step 4: Check if seqtk produced any output
    echo "Checking output files:"
    ls -l $PLUS_FILE_R1 $MINUS_FILE_R1 $PLUS_FILE_R2 $MINUS_FILE_R2

    # Step 5: Clean up intermediate files
    rm plus_strand_read_ids.txt minus_strand_read_ids.txt

    # Step 6: Compress the final output files
    echo "Compressing output files..."
    gzip $PLUS_FILE_R1
    gzip $MINUS_FILE_R1
    gzip $PLUS_FILE_R2
    gzip $MINUS_FILE_R2

    echo "Finished processing sample: $SAMPLE"
}

# Step 7: Loop over the samples A to K and process each sample sequentially
for SAMPLE in {A..K}; do
    echo "Processing sample: $SAMPLE"
    process_sample $SAMPLE
done

echo "All samples processed successfully!"
