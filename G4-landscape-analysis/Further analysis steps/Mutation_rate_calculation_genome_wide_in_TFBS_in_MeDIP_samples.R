library(openxlsx)
library(readxl)

master_dir <- "Path/GC_rich_TFBS_mutation_profile_genome_wide/Mutation_summary_of_GC_rich_TFBS_genomewide"
output_file <- "Path/GC_rich_TFBS_mutation_profile_genome_wide/final_CtoT_fraction_subtracted.xlsx"
samples <- c("Lc", "Ac", "Gg", "Hs", "Ev", "Input")

motif_files <- list.files(master_dir, pattern = "\\.xlsx$", full.names = TRUE)
result_list <- list()

for (motif_file in motif_files) {
  motif_name <- tools::file_path_sans_ext(basename(motif_file))
  fractions <- numeric(length(samples))
  names(fractions) <- samples
  for (i in seq_along(samples)) {
    sample <- samples[i]
    if (!(sample %in% excel_sheets(motif_file))) {
      fractions[i] <- 0
      next
    }
    # Read the sheet, skip the first row (empty), no headers
    df <- read.xlsx(motif_file, sheet = sample, startRow = 2, colNames = FALSE)
    if (nrow(df) == 0) {
      fractions[i] <- 0
      next
    }
    # Find where Mutation and Mutation_type tables start
    mut_idx <- which(df[,1] == "Mutation")
    type_idx <- which(df[,1] == "Mutation_type")
    if (length(mut_idx) == 0 || length(type_idx) == 0) {
      fractions[i] <- 0
      next
    }
    # Mutation table: from mut_idx+1 up to type_idx-2 (skip blank row before next table)
    mut_table <- df[(mut_idx+1):(type_idx-2), 1:2, drop=FALSE]
    colnames(mut_table) <- c("Mutation", "Count")
    mut_table <- mut_table[!is.na(mut_table$Mutation) & mut_table$Mutation != "", ]
    mut_table$Count <- as.numeric(mut_table$Count)
    c_to_t <- mut_table$Count[mut_table$Mutation == "C -> T"]
    c_to_t <- ifelse(length(c_to_t) == 0, 0, c_to_t)
    # Mutation_type table: from type_idx+1 to end
    type_table <- df[(type_idx+1):nrow(df), 1:2, drop=FALSE]
    colnames(type_table) <- c("Mutation_type", "Count")
    type_table <- type_table[!is.na(type_table$Mutation_type) & type_table$Mutation_type != "", ]
    type_table$Count <- as.numeric(type_table$Count)
    transition <- type_table$Count[type_table$Mutation_type == "Transition"]
    transition <- ifelse(length(transition) == 0, 0, transition)
    # Calculate fraction
    if (transition == 0) {
      fractions[i] <- 0
    } else {
      fractions[i] <- c_to_t / transition
    }
  }
  avg <- mean(fractions)
  avg_subtracted <- fractions - avg
  result_list[[motif_name]] <- avg_subtracted
}

# Build final data frame
final_df <- data.frame(
  Motif = names(result_list),
  do.call(rbind, result_list),
  row.names = NULL,
  check.names = FALSE
)
colnames(final_df)[-1] <- samples

# Replace any NA with 0
final_df[is.na(final_df)] <- 0

# Write to Excel
write.xlsx(final_df, output_file, rowNames = FALSE)
cat("Done! Final summary written to:\n", output_file, "\n")
