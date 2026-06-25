# 1. LOAD REQUIRED LIBRARIES
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(data.table) # Crucial for high-speed, row-wise coordinate expansions
})

# --- SETUP WORKING DIRECTORIES ---
base_output_dir <- "Path/G4_analysis/G4_motif_analysis/G4_in_promoters/Mammals/Homo_sapiens/Motif_specific_tsv_file"
prom_base_dir   <- file.path(base_output_dir, "Promoter")

# Input file paths (Exclusively for Promoter)
fimo_prom_file  <- "Path/G4_analysis/G4_motif_analysis/G4_in_promoters/Mammals/Homo_sapiens/fimo.tsv"
sig_prom_file   <- "Path/Significant_Motifs_Summary_Promoter.tsv"

# --- 2. LOAD SIGNIFICANT MOTIF LISTS ---
if (!file.exists(sig_prom_file)) {
  stop("[ERROR] Significant Promoter summary sheet is missing from your Desktop directory.")
}

sig_prom_ids <- read_tsv(sig_prom_file, show_col_types = FALSE) %>% 
  pull(`Motif ID`) %>% unique()

cat(sprintf("[INFO] Loaded %d significant Promoter motifs.\n", length(sig_prom_ids)))

# --- 3. HIGH-SPEED BED EXPORT CONVERSION FUNCTION ---
convert_and_save_single_base_bed <- function(dt_chunk, target_ids, context_base_dir) {
  if (nrow(dt_chunk) == 0) return(NULL)
  
  # Filter chunk to keep ONLY motifs that are significant for this context
  filtered_dt <- dt_chunk[motif_id %in% target_ids]
  if (nrow(filtered_dt) == 0) return(NULL)
  
  # Fail-safe coordinate parsing block
  if (any(str_detect(filtered_dt$sequence_name, ":"))) {
    # If colons exist, parse out Chromosome and Absolute Chunk Start
    filtered_dt[, c("raw_chr", "chunk_start") := tstrsplit(sequence_name, ":", keep = c(1, 2))]
    filtered_dt[, chunk_start := as.numeric(tstrsplit(chunk_start, "-", keep = 1)[[1]])]
  } else {
    # If NO colons exist, the sequence name is treated as the raw chromosome directly
    filtered_dt[, raw_chr := sequence_name]
    filtered_dt[, chunk_start := 1] 
  }
  
  # Standardize chromosome strings to include the native "chr" prefix safely
  filtered_dt[, chr := if_else(str_detect(raw_chr, "^chr"), raw_chr, paste0("chr", raw_chr))]
  
  # Calculate absolute genomic start and end positions on the chromosome
  filtered_dt[, absolute_start := chunk_start + start - 1]
  filtered_dt[, absolute_end   := chunk_start + stop - 1]
  
  unique_motifs <- unique(filtered_dt$motif_id)
  
  for (id in unique_motifs) {
    motif_subset <- filtered_dt[motif_id == id]
    
    # Fast single row-by-row sequence calculations
    single_base_dt <- motif_subset[, .(
      chr = chr,
      Genomic_Position = seq(absolute_start, absolute_end)
    ), by = 1:nrow(motif_subset)]
    
    # FIX: Swapped out the illegal standalone outside operator ':=' for standard assignment arrow '<-'
    output_dt <- single_base_dt[, .(chr, Genomic_Position)]
    
    # Sort data frame alphanumerically by Chromosome and numerically by Position
    # This prevents any downstream sorting problems during file intersections
    setorder(output_dt, chr, Genomic_Position)
    
    # Dynamic creation of a standalone nested directory named after the Motif ID
    motif_dir <- file.path(context_base_dir, id)
    dir.create(motif_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Define file output path as a .bed configuration
    output_file <- file.path(motif_dir, "single_base_positions.bed")
    
    # Output as clean tab-separated lines containing no text header rows
    fwrite(
      output_dt, 
      file = output_file, 
      sep = "\t", 
      append = file.exists(output_file), 
      col.names = FALSE
    )
  }
}

# --- 4. STREAM AND BLOCK-PROCESS STREAM FUNCTION ---
process_fimo_stream <- function(fimo_path, target_ids, context_base_dir, label) {
  if (!file.exists(fimo_path)) {
    stop(sprintf("[ERROR] Could not find fimo file for %s at path:\n%s", label, fimo_path))
  }
  
  cat(sprintf("[INFO] Processing %s stream from: %s\n", label, fimo_path))
  
  # Establish clear text file connection stream
  con <- file(fimo_path, "rt")
  
  # Read headers to preserve schema definitions
  header_line <- readLines(con, n = 1)
  
  chunk_size <- 500000
  has_data <- TRUE
  counter <- 0
  
  while (has_data) {
    lines <- readLines(con, n = chunk_size)
    
    if (length(lines) == 0) {
      has_data = FALSE
      next
    }
    
    # Load chunk block dynamically into memory as a data.table matrix
    chunk_dt <- fread(text = c(header_line, lines), sep = "\t", header = TRUE, showProgress = FALSE)
    
    # Process, sort, and save 
    convert_and_save_single_base_bed(chunk_dt, target_ids, context_base_dir)
    
    counter <- counter + length(lines)
    cat(sprintf("  Processed %s lines smoothly for %s context...\n", format(counter, big.mark=","), label))
  }
  
  close(con)
}

# --- 5. EXECUTE PIPELINE PROCESSING PASS ---
# Clean up old structural runs inside the Promoter folder to prevent cross-contamination
if (dir.exists(prom_base_dir)) {
  unlink(list.files(prom_base_dir, full.names = TRUE), recursive = TRUE)
} else {
  dir.create(prom_base_dir, showWarnings = FALSE, recursive = TRUE)
}

# Run the Promoter processing block exclusively
process_fimo_stream(fimo_prom_file, sig_prom_ids, prom_base_dir, "Promoter")

cat("\n[SUCCESS] Promoter-only pipeline execution complete! Alphanumeric-sorted single-base profiles generated successfully.\n")