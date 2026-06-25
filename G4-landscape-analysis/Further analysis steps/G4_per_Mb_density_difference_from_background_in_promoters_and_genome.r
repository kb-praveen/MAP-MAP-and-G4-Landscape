suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

# --- PATH CONFIGURATION ---
base_path <- "path/G4_motif_analysis/Results"

# 1. Promoter Paths
prom_g4_file   <- file.path(base_path, "G4_in_promoter/Motif_Density_per_Mb_G4_Promoters.tsv")
prom_free_file <- file.path(base_path, "G4_free_promoter/Motif_Density_per_Mb_G4_free_Promoters.tsv")

# 2. Genome-wide Paths
# Note: Using your previous naming conventions for the genomewide density files
gen_g4_file    <- file.path(base_path, "Genomewide_comparison/Motif_Density_per_Mb_G4s_Matrix.tsv")
gen_free_file  <- file.path(base_path, "Genomewide_G4_free/Motif_Density_per_Mb_G4free_Matrix.tsv")

# --- CALCULATION FUNCTION ---
calculate_diff <- function(g4_path, free_path, output_name) {
  if(!file.exists(g4_path) | !file.exists(free_path)) {
    message("[ERROR] Missing files for: ", output_name)
    return(NULL)
  }
  
  cat("[PROCESS] Calculating difference for:", output_name, "\n")
  
  # Load Data
  g4_df   <- read_tsv(g4_path, show_col_types = FALSE)
  free_df <- read_tsv(free_path, show_col_types = FALSE)
  
  # Ensure Motifs match
  common_motifs <- intersect(g4_df$Motif_ID, free_df$Motif_ID)
  g4_df   <- g4_df %>% filter(Motif_ID %in% common_motifs) %>% arrange(Motif_ID)
  free_df <- free_df %>% filter(Motif_ID %in% common_motifs) %>% arrange(Motif_ID)
  
  # Extract Metadata and Matrices
  motif_meta <- g4_df[, 1:2]
  g4_mat     <- as.matrix(g4_df[, -c(1,2)])
  free_mat   <- as.matrix(free_df[, -c(1,2)])
  
  # Ensure Species (Columns) match exactly
  common_species <- intersect(colnames(g4_mat), colnames(free_mat))
  g4_mat   <- g4_mat[, common_species]
  free_mat <- free_mat[, common_species]
  
  # Calculate Difference: G4 - Background
  diff_mat <- g4_mat - free_mat
  
  # Combine and Write
  final_df <- cbind(motif_meta, as.data.frame(diff_mat))
  write_tsv(final_df, file.path(base_path, output_name))
  cat("[SUCCESS] Saved to:", file.path(base_path, output_name), "\n")
}

# --- EXECUTE ---

# 1. Promoter Difference
calculate_diff(prom_g4_file, prom_free_file, "Motif_Density_Difference_Promoters_G4_minus_Free.tsv")

# 2. Genome-wide Difference
calculate_diff(gen_g4_file, gen_free_file, "Motif_Density_Difference_Genomewide_G4_minus_Free.tsv")

cat("\n[FINISH] Difference calculations complete.\n")