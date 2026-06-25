suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

# --- PATH CONFIGURATION ---
base_path    <- "Path/G4_motif_analysis"
results_path <- file.path(base_path, "Results")
length_file  <- file.path(results_path, "Promoter_Search_Space_Summary.tsv")

# Input Raw Count Files
g4_raw_file   <- file.path(base_path, "G4_in_promoters/Motif_counts_in_G4_containing_promoters_raw.tsv")
free_raw_file <- file.path(base_path, "Fimo_on_G4_free_promoters/Motif_counts_in_G4_free_promoters_raw.tsv")

# Define Output Directories
g4_out_dir   <- file.path(results_path, "G4_in_promoter")
free_out_dir <- file.path(results_path, "G4_free_promoter")

# Create directories if they don't exist
dir.create(g4_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(free_out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load the Exact Lengths
cat("[INFO] Loading search space summary...\n")
len_df <- read_tsv(length_file, show_col_types = FALSE)

# 2. Normalization Function
normalize_to_mb <- function(count_path, dataset_label, out_dir, out_filename) {
  if(!file.exists(count_path)) {
    message("[ERROR] File not found: ", count_path)
    return(NULL)
  }
  
  # Load raw counts - Removed check.names = FALSE
  counts <- read_tsv(count_path, show_col_types = FALSE)
  
  motif_meta <- counts[, 1:2]
  counts_mat <- as.matrix(counts[, -c(1,2)])
  
  # Extract lengths for this specific dataset
  size_subset <- len_df %>% 
    filter(Dataset == dataset_label) %>%
    select(Species, Search_Space_Mb)
  
  # Align species: ensure the divisor vector matches the matrix columns
  species_in_mat <- colnames(counts_mat)
  divisors <- size_subset$Search_Space_Mb[match(species_in_mat, species_in_mat)]
  # Correction: match against size_subset$Species
  divisors <- size_subset$Search_Space_Mb[match(species_in_mat, size_subset$Species)]
  
  # Handle potential zeros or NAs to avoid Inf
  divisors[is.na(divisors) | divisors == 0] <- NA
  
  # Normalization: Count / Mb
  cat("[PROCESS] Normalizing", dataset_label, "for", ncol(counts_mat), "species...\n")
  density_mat <- sweep(counts_mat, 2, divisors, "/")
  
  # Combine with metadata
  final_df <- cbind(motif_meta, as.data.frame(density_mat))
  
  # Save
  write_tsv(final_df, file.path(out_dir, out_filename))
  cat("[SUCCESS] Saved to:", file.path(out_dir, out_filename), "\n")
}

# 3. Execute for both datasets
normalize_to_mb(
  count_path   = g4_raw_file, 
  dataset_label = "G4_in_promoter", 
  out_dir      = g4_out_dir, 
  out_filename = "Motif_Density_per_Mb_G4_Promoters.tsv"
)

normalize_to_mb(
  count_path   = free_raw_file, 
  dataset_label = "G4_free_promoter", 
  out_dir      = free_out_dir, 
  out_filename = "Motif_Density_per_Mb_G4_free_Promoters.tsv"
)

cat("\n[FINISH] All promoter normalization tasks complete.\n")