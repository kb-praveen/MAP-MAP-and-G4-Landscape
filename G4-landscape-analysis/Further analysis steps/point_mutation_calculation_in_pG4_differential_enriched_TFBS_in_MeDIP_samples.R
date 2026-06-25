# 1. LOAD REQUIRED LIBRARIES
suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(data.table) # Essential for high-speed table reading
})

# --- SETUP WORKING DIRECTORIES ---
# MODIFICATION: Configured output to write directly into your Promoter folder
master_dir  <- "Path/Promoter"
output_file <- file.path(master_dir, "Final_all_mutations_summary_promoter.xlsx")
samples     <- c("Lc", "Ac", "Gg", "Hs", "Ev", "Input")

sig_prom_file <- "Path/Significant_Motifs_Summary_Promoter.tsv"

# List of all point mutations
all_mutations <- c(
  "A -> C", "A -> G", "A -> T",
  "C -> A", "C -> G", "C -> T",
  "G -> A", "G -> C", "G -> T",
  "T -> A", "T -> C", "T -> G"
)

# Transitions and transversions sets
transition_set <- c("A -> G", "G -> A", "C -> T", "T -> C")
transversion_set <- c(
  "A -> C", "A -> T",
  "C -> A", "C -> G",
  "G -> C", "G -> T",
  "T -> A", "T -> G"
)

# --- 2. LOAD LOOKUP INFORMATION ---
if (!file.exists(sig_prom_file)) {
  stop("[ERROR] Significant Promoter summary sheet is missing from your Desktop directory.")
}

# Read lookup data and select columns for mapping
metadata_df <- fread(sig_prom_file, sep = "\t", header = TRUE, showProgress = FALSE) %>%
  dplyr::select(
    `Motif ID`, 
    `Motif Name`, 
    `Enrichment status`
  ) %>%
  distinct()

# Fetch all unique motif subdirectories natively present in the Promoter folder
motif_dirs <- list.dirs(master_dir, full.names = TRUE, recursive = FALSE)
motif_names <- basename(motif_dirs)

if (length(motif_dirs) == 0) {
  stop("[ERROR] No motif subdirectories found inside your specified Promoter path.")
}

# Prepare final results repository list
all_mutation_results <- list()

# Initialize empty result matrices for each mutation type to compile data blocks smoothly
for (mutation in all_mutations) {
  mutation_matrix <- matrix(0, nrow = length(motif_dirs), ncol = length(samples))
  rownames(mutation_matrix) <- motif_names
  colnames(mutation_matrix) <- samples
  all_mutation_results[[mutation]] <- mutation_matrix
}

cat(sprintf("[INFO] Scanning through %d motif subdirectories for %d samples...\n", 
            length(motif_dirs), length(samples)))

# --- 3. LOOP THROUGH DIRECTORIES AND CALCULATE FRACTIONS ---
for (d in seq_along(motif_dirs)) {
  current_dir  <- motif_dirs[d]
  current_motif <- motif_names[d]
  
  for (sample in samples) {
    sample_file <- file.path(current_dir, paste0(sample, "_mutation_profile_single_base.tsv"))
    
    if (!file.exists(sample_file) || file.info(sample_file)$size == 0) next
    
    dt <- fread(sample_file, sep = "\t", header = TRUE, showProgress = FALSE)
    if (nrow(dt) == 0) next
    
    mut_col  <- grep("^Mutation(_|$)", colnames(dt), value = TRUE)[1]
    type_col <- grep("^Mutation_type(_|$)", colnames(dt), value = TRUE)[1]
    
    if (is.na(mut_col) || is.na(type_col)) next
    
    mut_counts  <- table(dt[[mut_col]])
    type_counts <- table(dt[[type_col]])
    
    total_transitions  <- sum(type_counts[names(type_counts) == "Transition"], na.rm = TRUE)
    total_transversions <- sum(type_counts[names(type_counts) == "Transversion"], na.rm = TRUE)
    
    for (mutation in all_mutations) {
      mutation_count <- sum(mut_counts[names(mut_counts) == mutation], na.rm = TRUE)
      
      is_transition <- mutation %in% transition_set
      denom_count <- ifelse(is_transition, total_transitions, total_transversions)
      
      fraction <- ifelse(denom_count == 0, 0, mutation_count / denom_count)
      
      all_mutation_results[[mutation]][current_motif, sample] <- fraction
    }
  }
}

# --- 4. APPLY NORMALIZATION, MERGE METADATA, AND EXPORT ---
wb <- createWorkbook()

for (mutation in all_mutations) {
  mat <- all_mutation_results[[mutation]]
  
  # Row-wise normalization: Subtract motif mean from each sample (across motifs row)
  for (r in 1:nrow(mat)) {
    avg <- mean(mat[r, ], na.rm = TRUE)
    if (is.nan(avg)) avg <- 0
    mat[r, ] <- mat[r, ] - avg
  }
  
  # Format into initial data frame structure layout
  mat_df <- data.frame(
    `Motif ID` = rownames(mat),
    mat,
    row.names = NULL,
    check.names = FALSE
  )
  
  # MODIFICATION: Join with the metadata lookup table to bring in Motif Name and Enrichment status
  final_df <- left_join(mat_df, metadata_df, by = "Motif ID") %>%
    # Re-order columns beautifully so annotations sit directly next to the Motif ID
    relocate(`Motif Name`, `Enrichment status`, .after = `Motif ID`)
  
  # Replace characters to fit safe Excel sheet nomenclature requirements
  sheet_name <- gsub(" -> ", "_to_", mutation)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, final_df)
}

# Save spreadsheet
saveWorkbook(wb, output_file, overwrite = TRUE)
cat("\n✅ SUCCESS! Promoter summary successfully calculated, annotated, and written inside folder to:\n", output_file, "\n")