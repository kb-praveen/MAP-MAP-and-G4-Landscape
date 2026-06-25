library(openxlsx)
library(data.table) # Optimized for processing high-throughput genomic data tables

# --- SETUP WORKING DIRECTORIES ---
master_dir  <- "Path/Promoter"
output_file <- file.path(master_dir, "Final_raw_mutation_counts_with_summaries_promoter.xlsx")
samples     <- c("Lc", "Ac", "Gg", "Hs", "Ev", "Input")

# Structural list of items to compute as final table rows
all_mutations   <- c(
  "A -> C", "A -> G", "A -> T",
  "C -> A", "C -> G", "C -> T",
  "G -> A", "G -> C", "G -> T",
  "T -> A", "T -> C", "T -> G"
)
summary_metrics <- c("Transition", "Transversion", "No mutation", "Others", "Total_bases")
all_target_rows <- c(all_mutations, summary_metrics)

# Fetch all unique motif subdirectories natively present in the Promoter folder
motif_dirs  <- list.dirs(master_dir, full.names = TRUE, recursive = FALSE)
motif_names <- basename(motif_dirs)

if (length(motif_dirs) == 0) {
  stop("[ERROR] No motif subdirectories found inside your specified Promoter path.")
}

# Create the master workbook container
wb <- createWorkbook()

cat(sprintf("[INFO] Scanning through %d motif subdirectories for raw counts parsing...\n", 
            length(motif_dirs)))

# --- PROCESS EACH MOTIF SUBDIRECTORY AS A SEPARATE SHEET ---
for (d in seq_along(motif_dirs)) {
  current_dir  <- motif_dirs[d]
  current_motif <- motif_names[d]
  
  # Ensure the Excel tab name doesn't exceed 31 characters (Excel limitation)
  sheet_tab_name <- substr(current_motif, 1, 31)
  
  # Initialize an isolated blank matrix for this specific motif
  raw_counts_matrix <- matrix(0, nrow = length(all_target_rows), ncol = length(samples))
  rownames(raw_counts_matrix) <- all_target_rows
  colnames(raw_counts_matrix) <- samples
  
  for (sample in samples) {
    sample_file <- file.path(current_dir, paste0(sample, "_mutation_profile_single_base.tsv"))
    
    # Skip if file is missing or completely empty
    if (!file.exists(sample_file) || file.info(sample_file)$size == 0) next
    
    # Read the data file at maximum throughput execution speeds
    dt <- fread(sample_file, sep = "\t", header = TRUE, showProgress = FALSE)
    if (nrow(dt) == 0) next
    
    # Matches columns starting precisely with "Mutation_[Sample]" or "Mutation_type_[Sample]"
    mut_col_name  <- grep(paste0("^Mutation_", sample, "$"), colnames(dt), value = TRUE)[1]
    type_col_name <- grep(paste0("^Mutation_type_", sample, "$"), colnames(dt), value = TRUE)[1]
    
    # Fallback to general patterns if sample-specific suffix match fails
    if (is.na(mut_col_name)) {
      mut_col_name  <- grep("^Mutation", colnames(dt), value = TRUE)[1]
    }
    if (is.na(type_col_name)) {
      type_col_name <- grep("^Mutation_type", colnames(dt), value = TRUE)[1]
    }
    
    if (is.na(mut_col_name) || is.na(type_col_name)) next
    
    # Extract clean vectors as characters to guarantee strings map accurately
    mut_vector  <- trimws(as.character(dt[[mut_col_name]]))
    type_vector <- trimws(as.character(dt[[type_col_name]]))
    
    # Generate frequency tables
    mut_counts  <- table(mut_vector)
    type_counts <- table(type_vector)
    
    # 1. Populate the 12 primary single-base point mutations
    for (mutation in all_mutations) {
      count_val <- sum(mut_counts[names(mut_counts) == mutation], na.rm = TRUE)
      raw_counts_matrix[mutation, sample] <- count_val
    }
    
    # 2. Populate Structural Summaries
    raw_counts_matrix["Transition", sample]   <- sum(type_counts[names(type_counts) == "Transition"], na.rm = TRUE)
    raw_counts_matrix["Transversion", sample] <- sum(type_counts[names(type_counts) == "Transversion"], na.rm = TRUE)
    
    # MODIFICATION: Dynamic sum evaluating both "No mutation" text AND "-" character entries as No Mutation
    no_mut_text_count <- sum(mut_counts[names(mut_counts) == "No mutation"], na.rm = TRUE)
    no_mut_dash_count <- sum(mut_counts[names(mut_counts) == "-"], na.rm = TRUE)
    
    raw_counts_matrix["No mutation", sample]  <- no_mut_text_count + no_mut_dash_count
    
    # Extract "Others" summary count
    raw_counts_matrix["Others", sample]       <- sum(type_counts[names(type_counts) == "Others"], na.rm = TRUE)
    
    # 3. Calculate Total Bases (Matches total rows in the input file)
    raw_counts_matrix["Total_bases", sample]  <- nrow(dt)
  }
  
  # Convert the finalized matrix into a cleanly formatted dataframe for export
  motif_df <- data.frame(
    `Metric Type` = rownames(raw_counts_matrix),
    raw_counts_matrix,
    row.names = NULL,
    check.names = FALSE
  )
  
  # Append this motif as an isolated tab inside the excel file
  addWorksheet(wb, sheet_tab_name)
  writeData(wb, sheet_tab_name, motif_df)
}

# Save structured workbook
saveWorkbook(wb, output_file, overwrite = TRUE)
cat("\n✅ SUCCESS! Suffix-aware script updated to capture hyphens (-) as unmutated records.\nOutput file written to:\n", output_file, "\n")