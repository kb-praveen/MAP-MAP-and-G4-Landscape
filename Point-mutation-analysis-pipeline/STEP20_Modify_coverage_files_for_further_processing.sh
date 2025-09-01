#!/bin/bash

# Loop over all files ending with .tsv.gz in the current directory
for file in *.tsv.gz; do
    # Unzip, modify the header and columns using awk, then gzip the output
    zcat "$file" | awk 'NR==1 {print "Chromosome\tPosition\tCoverage\tBase_calls"} NR>1 {print $1 "\t" $2 "\t" $3 "\t" $4}' | gzip > "modified_$(basename "$file")"
done
