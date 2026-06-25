# ==============================================================================
# PAIRWISE AMINO ACID PERCENTAGE IDENTITY CALCULATOR (NATIVE APE IMPLEMENTATION)
# ==============================================================================
library(ape)
library(dplyr)

# 1. SETUP YOUR PATHS
alignment_path   <- "/Users/praveenbishnoi/Downloads/CGGBP1_alignment.fasta"
output_sim_path  <- "/Users/praveenbishnoi/Downloads/CGGBP1_percent_identity_human.tsv"
reference_anchor <- "Homo_sapiens"

# 2. READ THE ALIGNMENT FILE NATIVELY WITH APE
# read.FASTA loads the file as a specialized biological list matrix
align_raw <- read.FASTA(alignment_path, type = "AA")
align_matrix <- as.matrix(align_raw)

# Extract raw names present in the sequence dataset rows
raw_names <- rownames(align_matrix)

# Quick sanity check: verify that your reference anchor is present
if (!(reference_anchor %in% raw_names)) {
  # If there is a subspecies or database tag in the file headers, catch it dynamically
  matching_row <- grep(reference_anchor, raw_names, value = TRUE)
  if(length(matching_row) > 0) {
    reference_anchor <- matching_row[1]
  } else {
    stop("Could not find the reference anchor species in your alignment file headers!")
  }
}

# 3. ISOLATE THE HUMAN REFERENCE SEQUENCE VECTOR
human_seq <- align_matrix[reference_anchor, ]

# 4. LOOP OVER ROWS TO CALCULATE MATCH PERCENTAGES
identity_results <- data.frame(
  TreeNames = raw_names,
  Percent_Identity = NA,
  stringsAsFactors = FALSE
)

# Convert character raw bytes to clear text matrices if ape read format shifts
human_chars <- as.character(human_seq)

for (i in seq_along(raw_names)) {
  current_sp   <- raw_names[i]
  current_chars <- as.character(align_matrix[current_sp, ])
  
  # EXCLUSION CRITERIA FOR GAPS:
  # Filter out structural alignment position columns where both rows have an introduced gap ('-')
  valid_positions <- !(human_chars == "-" & current_chars == "-")
  
  filt_human   <- human_chars[valid_positions]
  filt_current <- current_chars[valid_positions]
  
  # Calculate exact amino acid matches across valid columns
  matches <- sum(filt_human == filt_current)
  total_length <- length(filt_human)
  
  # Convert to a clean percentage value (0% - 100%)
  identity_results$Percent_Identity[i] <- (matches / total_length) * 100
}

# 5. TAXONOMY CLEAN-UP FOR TREEDATA COMPATIBILITY
# Strip out subspecies name variations to ensure 1:1 matching with your tree tips
identity_results$TreeNames <- identity_results$TreeNames %>%
  gsub("Canis_lupus_familiaris", "Canis_lupus", .) %>%
  gsub("Bison_bison_bison", "Bison_bison", .) %>%
  gsub("Gorilla_gorilla_gorilla", "Gorilla_gorilla", .) %>%
  gsub("Chrysemys_picta_bellii", "Chrysemys_picta", .) %>%
  gsub("Aquila_chrysaetos_chrysaetos", "Aquila_chrysaetos", .) %>%
  gsub("Strix_occidentalis_caurina", "Strix_occidentalis", .)

# 6. EXPORT LOGGED VALUES TO TAB-SEPARATED TSV FILE
write.table(identity_results, 
            output_sim_path, 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

cat("\n*** Success! Percentage identity metrics saved cleanly to:", output_sim_path, "***\n")