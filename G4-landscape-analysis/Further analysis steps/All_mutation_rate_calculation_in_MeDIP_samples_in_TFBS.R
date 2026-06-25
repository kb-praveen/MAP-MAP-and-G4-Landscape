library(openxlsx)
library(readxl)

master_dir <- "PathGC_rich_TFBS_mutation_profile_genome_wide/Mutation_summary_of_GC_rich_TFBS_genomewide"
output_file <- "PathGC_rich_TFBS_mutation_profile_genome_wide/Final_all_mutations_summary_genomewide.xlsx"
samples <- c("Lc", "Ac", "Gg", "Hs", "Ev", "Input")

# List of all point mutations
all_mutations <- c(
  "A -> C", "A -> G", "A -> T",
  "C -> A", "C -> G", "C -> T",
  "G -> A", "G -> C", "G -> T",
  "T -> A", "T -> C", "T -> G"
)

# Transitions and transversions explicitly defined
transition_set <- c("A -> G", "G -> A", "C -> T", "T -> C")
transversion_set <- c(
  "A -> C", "A -> T",
  "C -> A", "C -> G",
  "G -> C", "G -> T",
  "T -> A", "T -> G"
)

# Prepare result store
all_mutation_results <- list()
motif_files <- list.files(master_dir, pattern = "\\.xlsx$", full.names = TRUE)

for (mutation in all_mutations) {
  mutation_matrix <- matrix(0, nrow = length(motif_files), ncol = length(samples))
  rownames(mutation_matrix) <- tools::file_path_sans_ext(basename(motif_files))
  colnames(mutation_matrix) <- samples
  
  for (f in seq_along(motif_files)) {
    motif_file <- motif_files[f]
    motif_name <- tools::file_path_sans_ext(basename(motif_file))
    
    for (i in seq_along(samples)) {
      sample <- samples[i]
      
      if (!(sample %in% excel_sheets(motif_file))) next
      
      df <- read.xlsx(motif_file, sheet = sample, startRow = 2, colNames = FALSE)
      if (nrow(df) == 0) next
      
      mut_idx <- which(df[, 1] == "Mutation")
      type_idx <- which(df[, 1] == "Mutation_type")
      if (length(mut_idx) == 0 || length(type_idx) == 0) next
      
      mut_table <- df[(mut_idx + 1):(type_idx - 2), 1:2, drop = FALSE]
      colnames(mut_table) <- c("Mutation", "Count")
      mut_table <- mut_table[!is.na(mut_table$Mutation) & mut_table$Mutation != "", ]
      mut_table$Count <- as.numeric(mut_table$Count)
      
      type_table <- df[(type_idx + 1):nrow(df), 1:2, drop = FALSE]
      colnames(type_table) <- c("Mutation_type", "Count")
      type_table <- type_table[!is.na(type_table$Mutation_type) & type_table$Mutation_type != "", ]
      type_table$Count <- as.numeric(type_table$Count)
      
      # Get mutation count
      mutation_count <- mut_table$Count[mut_table$Mutation == mutation]
      mutation_count <- ifelse(length(mutation_count) == 0, 0, mutation_count)
      
      # Denominator: transition or transversion
      is_transition <- mutation %in% transition_set
      is_transversion <- mutation %in% transversion_set
      denominator_type <- ifelse(is_transition, "Transition", 
                                 ifelse(is_transversion, "Transversion", NA))
      denom_count <- type_table$Count[type_table$Mutation_type == denominator_type]
      denom_count <- ifelse(length(denom_count) == 0, 0, denom_count)
      fraction <- ifelse(denom_count == 0, 0, mutation_count / denom_count)
      mutation_matrix[motif_name, sample] <- fraction
    }
  }
  
  # Subtract motif mean from each sample (across motifs)
  for (r in 1:nrow(mutation_matrix)) {
    avg <- mean(mutation_matrix[r, ])
    mutation_matrix[r, ] <- mutation_matrix[r, ] - avg
  }
  
  mutation_df <- data.frame(
    Motif = rownames(mutation_matrix),
    mutation_matrix,
    row.names = NULL,
    check.names = FALSE
  )
  all_mutation_results[[mutation]] <- mutation_df
}

# Write to Excel, each mutation in its own sheet
wb <- createWorkbook()
for (mutation in names(all_mutation_results)) {
  sheet_name <- gsub(" -> ", "_to_", mutation)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, all_mutation_results[[mutation]])
}

saveWorkbook(wb, output_file, overwrite = TRUE)
cat("✅ Done! Final mutation summary written to:\n", output_file, "\n")
