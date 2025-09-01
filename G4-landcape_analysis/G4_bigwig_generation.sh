#!/bin/bash

# --- Set the file paths below ---
WINDOW_BED="/path/to/your/10bp_windows.bed"           # BED file of 10bp genome intervals
G4_BED="/path/to/your/g4_regions.bed"                 # BED file of G4-forming regions
BEDGRAPH="/path/to/output_g4_signal.bedGraph"         # Output bedGraph
BIGWIG="/path/to/output_g4_signal.bw"                 # Output bigWig
CHROMSIZES="/path/to/hg38.chrom.sizes"                # Chromosome sizes (tab-separated)

# --- Calculate coverage: fraction of each window overlapped by G4s ---
bedtools coverage -a "$WINDOW_BED" -b "$G4_BED" > "${BEDGRAPH}.tmp"

# --- Create bedGraph with percentage signal (values 0.00–100.00) ---
awk 'BEGIN {OFS="\t"} {printf "%s\t%s\t%s\t%.2f\n", $1, $2, $3, $7*100}' "${BEDGRAPH}.tmp" > "$BEDGRAPH"
rm -f "${BEDGRAPH}.tmp"

# --- Sort (multi-threaded) for bedGraphToBigWig compatibility ---
sort --parallel=22 -k1,1 -k2,2n "$BEDGRAPH" -o "$BEDGRAPH"

# --- Convert bedGraph to bigWig ---
bedGraphToBigWig "$BEDGRAPH" "$CHROMSIZES" "$BIGWIG"

echo "Done."
echo "  bedGraph: $BEDGRAPH"
echo "  bigWig:   $BIGWIG"
