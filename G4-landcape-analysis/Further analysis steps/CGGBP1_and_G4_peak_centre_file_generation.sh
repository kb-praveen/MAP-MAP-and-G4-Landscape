#!/bin/bash

# 1. Define File Paths
G4_PRED="/Users/praveenbishnoi/Desktop/Thesis_rebutal/G4_analysis/predicted_G4s.bed"
G4_CHIP="/Users/praveenbishnoi/Desktop/Thesis_rebutal/G4_analysis/GSM6122747_G4_ChIP_broadpeaks.bed"
PEAK_FILE="/Users/praveenbishnoi/Desktop/Thesis_rebutal/G4_analysis/CGGBP1_narrow_peaks_all_merged.bed"
OUTPUT_DIR="/Users/praveenbishnoi/Desktop/Thesis_rebutal/G4_analysis"

# Define Output Filenames
CENTERS_BED="${OUTPUT_DIR}/CGGBP1_peak_centers.bed"
DIST_PRED_OUT="${OUTPUT_DIR}/CGGBP1_to_Predicted_G4_distances.txt"
DIST_CHIP_OUT="${OUTPUT_DIR}/CGGBP1_to_G4_ChIP_distances.txt"

echo "Starting dual spatial analysis..."

# 2. Step 1: Calculate CGGBP1 Peak Centers
# Midpoints ensure the distance reflects the core binding event 
awk 'BEGIN {OFS="\t"} {midpoint = int(($2 + $3) / 2); print $1, midpoint, midpoint + 1}' "$PEAK_FILE" > "$CENTERS_BED"

# 3. Step 2: Calculate Distances
if command -v bedtools &> /dev/null; then
    # Comparison 1: CGGBP1 to Predicted G4s
    echo "Processing: CGGBP1 vs Predicted G4s..."
    bedtools closest -a "$CENTERS_BED" -b "$G4_PRED" -d > "$DIST_PRED_OUT"
    
    # Comparison 2: CGGBP1 to Experimental G4 ChIP peaks
    echo "Processing: CGGBP1 vs G4 ChIP peaks..."
    bedtools closest -a "$CENTERS_BED" -b "$G4_CHIP" -d > "$DIST_CHIP_OUT"
    
    echo "Done. Results saved to:"
    echo "  1. $DIST_PRED_OUT"
    echo "  2. $DIST_CHIP_OUT"
else
    echo "Error: bedtools not found."
    exit 1
fi