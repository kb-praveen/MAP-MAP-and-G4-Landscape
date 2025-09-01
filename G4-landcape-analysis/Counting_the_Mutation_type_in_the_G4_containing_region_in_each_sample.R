# Set working directory
setwd("Path_to_the_directory")

# Get list of all TSV files
input_files <- list.files(pattern = "\\.tsv$")

# Process each file
for (input_file in input_files) {
  # Read data (ensure column names match your file)
  data <- read.delim(input_file, header = TRUE, sep = "\t", check.names = FALSE)
  
  # Generate counts (using Mutation column as per your original data)
  mutation_counts <- as.data.frame(table(data$Mutation))
  colnames(mutation_counts) <- c("Mutation", "Count")
  
  mutation_type_counts <- as.data.frame(table(data$Mutation_type))
  colnames(mutation_type_counts) <- c("Mutation_type", "Count")
  
  # Create output filename
  output_file <- sub("\\.tsv$", "_counts.tsv", input_file)
  
  # Create file connection
  con <- file(output_file, open = "w")
  
  # Write header and content
  writeLines(paste("## Analysis for:", input_file), con)
  writeLines("\n## Mutation Counts", con)
  write.table(mutation_counts, con, sep = "\t", row.names = FALSE, quote = FALSE)
  
  writeLines("\n## Mutation Type Counts", con)
  write.table(mutation_type_counts, con, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Close connection
  close(con)
}

# Verification message
cat("Processed", length(input_files), "files. Output files created:\n", paste(list.files(pattern = "_counts\\.tsv$"), collapse = "\n "))
