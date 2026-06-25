#!/bin/bash

# Root path tracking your decoupled multi-tier directories framework
master_dir="Path/FIMO/TFBS_in_pG4/Motif_specific_tsv_file"

declare -A mutation_files
mutation_files=(
    [Ac]="Path/Genome_wide_C_to_T/Fastp_files_all/Everything_on_plus_strand/Merged_fastq_files/Coverage_between_4_and_1k/Base_calls_to_unique_strings/Lc_Ac_Gg_Hs_Ev_Input/Mutation_rate_raw_files/TSV_file/Ac_mutation_profile.tsv.gz"
    [Ev]="Path/Genome_wide_C_to_T/Fastp_files_all/Everything_on_plus_strand/Merged_fastq_files/Coverage_between_4_and_1k/Base_calls_to_unique_strings/Lc_Ac_Gg_Hs_Ev_Input/Mutation_rate_raw_files/TSV_file/Ev_mutation_profile.tsv.gz"
    [Gg]="Path/Genome_wide_C_to_T/Fastp_files_all/Everything_on_plus_strand/Merged_fastq_files/Coverage_between_4_and_1k/Base_calls_to_unique_strings/Lc_Ac_Gg_Hs_Ev_Input/Mutation_rate_raw_files/TSV_file/Gg_mutation_profile.tsv.gz"
    [Hs]="Path/Genome_wide_C_to_T/Fastp_files_all/Everything_on_plus_strand/Merged_fastq_files/Coverage_between_4_and_1k/Base_calls_to_unique_strings/Lc_Ac_Gg_Hs_Ev_Input/Mutation_rate_raw_files/TSV_file/Hs_mutation_profile.tsv.gz"
    [Input]="Path/Genome_wide_C_to_T/Fastp_files_all/Everything_on_plus_strand/Merged_fastq_files/Coverage_between_4_and_1k/Base_calls_to_unique_strings/Lc_Ac_Gg_Hs_Ev_Input/Mutation_rate_raw_files/TSV_file/Input_mutation_profile.tsv.gz"
    [Lc]="Path/Genome_wide_C_to_T/Fastp_files_all/Everything_on_plus_strand/Merged_fastq_files/Coverage_between_4_and_1k/Base_calls_to_unique_strings/Lc_Ac_Gg_Hs_Ev_Input/Mutation_rate_raw_files/TSV_file/Lc_mutation_profile.tsv.gz"
)

# Configured parallelization limits to explicitly deploy across 22 processing cores
max_jobs=22 

# MODIFICATION: Restricted finding/counting targets strictly to the nested 'Promoter' path
total_motifs=$(find "$master_dir" -mindepth 2 -maxdepth 2 -type d -path "*/Promoter/*" | wc -l | awk '{print $1}')
echo "[INFO] Found a total of $total_motifs unique Promoter motif targets inside workspace."

for key in Ac Ev Gg Hs Input Lc; do
    mutation_file_compressed="${mutation_files[$key]}"
    echo "---------------------------------------------------------"
    echo "Processing Organism Array for Key: $key"
    echo "Source Track: $mutation_file_compressed"
    echo "---------------------------------------------------------"

    if [ ! -f "$mutation_file_compressed" ]; then
        echo "[ERROR] File not found: $mutation_file_compressed. Skipping key $key."
        continue
    fi

    # Create a temporary file for the uncompressed mutation data
    temp_mutation_file="$(mktemp)"
    zcat "$mutation_file_compressed" > "$temp_mutation_file"

    # Fetch the exact native header string from line 1 of the uncompressed file
    native_header="$(head -n 1 "$temp_mutation_file")"

    job_count=0
    finished_count=0 # Progress tracker initialization for current organism loop pass
    
    # MODIFICATION: Changed path iteration step to look exclusively within the 'Promoter' subdirectory
    for subdir in "$master_dir"/Promoter/*/; do
        subdir="${subdir%/}"
        
        (
            bed_file="${subdir}/single_base_positions.bed"
            output_file="${subdir}/${key}_mutation_profile_single_base.tsv"

            if [ ! -f "$bed_file" ]; then
                exit
            fi

            # Clean and sanitize the BED file on-the-fly to explicitly drop invisible \r characters
            # Then use a highly secure formatting regex check inside awk to ensure matches map precisely
            awk -F'\t' '
            NR==FNR {
                # Process cleaned single_base_positions.bed file
                gsub(/\r/, "", $1)
                gsub(/\r/, "", $2)
                if ($1 != "" && $2 != "") {
                    # Strip spaces to ensure crisp matching keys
                    gsub(/[ \t]+$/, "", $1)
                    gsub(/[ \t]+$/, "", $2)
                    idx = $1 "::" $2
                    key[idx] = 1
                }
                next
            }
            FNR==1 {
                # Ignore the header row of the mutation profile
                next
            }
            {
                # Process temp_mutation_file row-by-row safely
                m_idx = $1 "::" $2
                if (m_idx in key) {
                    print $0
                }
            }
            ' <(tr -d '\r' < "$bed_file") "$temp_mutation_file" > "$output_file"

            # Check if any matching entries were actually found before appending headers
            if [ -s "$output_file" ]; then
                tmp_file="${output_file}.tmp"
                mv "$output_file" "$tmp_file"
                echo -e "$native_header" > "$output_file"
                cat "$tmp_file" >> "$output_file"
                rm "$tmp_file"
            else
                # Keep the file clean with just the header if no exact chromosome position overlaps exist
                echo -e "$native_header" > "$output_file"
            fi
        ) &
        
        ((job_count++))
        ((finished_count++))
        
        # Real-time verbose notifications showing exact motif process completions inside the Promoter context
        echo "Progress Update: [Processed $finished_count / $total_motifs Promoter locations] -> Current: $(basename "$subdir")"
        
        if (( job_count >= max_jobs )); then
            wait -n
            ((job_count--))
        fi
    done

    wait  # Sync current threads step before processing next species array block
    rm "$temp_mutation_file"
    echo ">>> Completed all Promoter intersections for organism track: $key <<<"
done

echo "All Promoter contexts, mutation profiles, and subfolders processed completely across 22 cores."