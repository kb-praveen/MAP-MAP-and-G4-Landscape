library(pqsfinder)
library(Biostrings)
library(BiocParallel)
library(stringr)
library(tools)
library(openxlsx)

# ---------------------- USER INPUTS -----------------------
input_dir <- "Path_to_input_genomes_directory/Taxa"
output_dir <- "Path_to_output_directory/Taxa"
num_cores <- 20

# ---------------------- DIRECTORIES SYNC ------------------
bed_dir    <- file.path(output_dir, "bed_files")
tsv_dir    <- file.path(output_dir, "tsv_files")
fasta_dir  <- file.path(output_dir, "fasta_files")
summary_csv <- file.path(output_dir, "G4_summary.csv")

dir.create(bed_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tsv_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fasta_dir, showWarnings = FALSE, recursive = TRUE)

# ------- Initialize summary file with header if it doesn't exist ------
if (!file.exists(summary_csv)) {
  write.table(data.frame(
    Species=character(), 
    Genome_Size_bp=numeric(), 
    GC_content_of_Genome=numeric(),
    G4_Predicted=integer(),
    Percent_Genome_Covered_by_G4=numeric(), 
    G4_per_Mb=numeric(),
    Average_G4_Length=numeric(), 
    Average_G4_Score=numeric(),
    GC_Content_Percent=numeric(), 
    Num_Chromosomes=integer()
  ),
  file = summary_csv, sep = ",", row.names = FALSE, col.names = TRUE
  )
}

# ---------------------- PARALLEL SETUP -------------------
param <- MulticoreParam(workers = num_cores)
fasta_files <- list.files(input_dir, pattern = "\\.fa(sta)?$", full.names = TRUE)

# ---------------------- G4 PROCESSING FUNCTION -----------
process_fasta <- function(fasta_file) {
  sample_name <- file_path_sans_ext(basename(fasta_file))
  sample_name <- sub("\\.fa(sta)?$", "", sample_name)
  
  output_bed <- file.path(bed_dir, paste0(sample_name, ".bed"))
  output_tsv <- file.path(tsv_dir, paste0(sample_name, ".tsv"))
  output_fasta <- file.path(fasta_dir, paste0(sample_name, ".fa"))
  
  if (file.exists(output_bed) && file.exists(output_tsv) && file.exists(output_fasta)) {
    cat("✅ Already Processed:", sample_name, "\n")
    return(NULL)
  }
  
  cat("🔍 Processing:", sample_name, "\n")
  
  sequences <- readDNAStringSet(fasta_file)
  if (length(sequences) == 0) return(NULL)
  
  genome_size <- sum(width(sequences))
  num_chroms <- length(sequences)
  
  # --- Calculate GC content of the entire genome
  genome_cat_seq <- unlist(sequences)
  base_counts <- alphabetFrequency(genome_cat_seq, baseOnly=TRUE)
  gc_genome <- sum(base_counts[c("G", "C")])
  acgt_total <- sum(base_counts[c("A", "C", "G", "T")])
  gc_genome_pct <- round(gc_genome / acgt_total * 100, 2)
  
  # Open file connections for line-wise writing
  bed_conn <- file(output_bed, open = "w")
  tsv_conn <- file(output_tsv, open = "w")
  fasta_conn <- file(output_fasta, open = "w")
  
  # Write TSV header
  writeLines("chrom\tstart\tend\tscore\tstrand\tsequence\tGC_content", tsv_conn)
  
  total_g4s <- 0
  total_g4_bp <- 0
  all_scores <- c()
  all_gcs <- c()
  all_widths <- c()
  
  for (i in seq_along(sequences)) {
    res <- pqsfinder(sequences[[i]], min_score = 52)
    if (length(res) == 0) next
    
    starts <- start(res)
    ends <- end(res)
    scores <- score(res)
    widths_ <- width(res)
    seq_id <- names(sequences)[i]
    g4_seqs <- as.character(Views(sequences[[i]], starts, ends))
    
    gc_vals <- round((str_count(g4_seqs, "G") + str_count(g4_seqs, "C")) / nchar(g4_seqs) * 100, 2)
    
    for (j in seq_along(starts)) {
      bed_line <- sprintf("%s\t%d\t%d\tG4_%d\t%d\t.", seq_id, starts[j]-1, ends[j], j, scores[j])
      writeLines(bed_line, bed_conn)
      
      tsv_line <- sprintf("%s\t%d\t%d\t%d\t.\t%s\t%.2f", seq_id, starts[j], ends[j], scores[j], g4_seqs[j], gc_vals[j])
      writeLines(tsv_line, tsv_conn)
      
      fasta_header <- sprintf(">%s:%d-%d", seq_id, starts[j], ends[j])
      writeLines(fasta_header, fasta_conn)
      writeLines(g4_seqs[j], fasta_conn)
    }
    
    total_g4s <- total_g4s + length(res)
    total_g4_bp <- total_g4_bp + sum(widths_)
    all_scores <- c(all_scores, scores)
    all_gcs <- c(all_gcs, gc_vals)
    all_widths <- c(all_widths, widths_)
  }
  
  close(bed_conn); close(tsv_conn); close(fasta_conn)
  
  if (total_g4s == 0) {
    cat("⚠️ No G4s found in", sample_name, "\n")
    return(NULL)
  }
  
  summary_row <- data.frame(
    Species = sample_name,
    Genome_Size_bp = genome_size,
    GC_content_of_Genome = gc_genome_pct,
    G4_Predicted = total_g4s,
    Percent_Genome_Covered_by_G4 = round((total_g4_bp / genome_size) * 100, 4),
    G4_per_Mb = round(total_g4s / (genome_size / 1e6), 2),
    Average_G4_Length = round(mean(all_widths), 2),
    Average_G4_Score = round(mean(all_scores), 2),
    GC_Content_Percent = round(mean(all_gcs), 2),
    Num_Chromosomes = num_chroms,
    stringsAsFactors = FALSE
  )
  
  # Append to summary file immediately
  write.table(
    summary_row, file = summary_csv,
    sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE
  )
  
  cat("✅ Done:", sample_name, "\n")
  return(NULL)
}

# ---------------------- EXECUTE PARALLEL ------------------
bplapply(fasta_files, FUN = process_fasta, BPPARAM = MulticoreParam(num_cores))

cat("\n🎉 All FASTA files processed! Summary CSV:\n", summary_csv, "\n")

# ---------------------- CONVERT CSV TO EXCEL --------------
summary_df <- read.csv(summary_csv)
summary_xlsx <- sub("\\.csv$", ".xlsx", summary_csv)
write.xlsx(summary_df, file = summary_xlsx, rowNames = FALSE)
cat("✅ Excel file saved at:", summary_xlsx, "\n")
