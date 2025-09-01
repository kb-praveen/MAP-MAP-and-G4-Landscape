library(data.table)
library(R.utils)

input_dir <- "/home/umash/MeDIP_raw_data/Genome_wide_C_to_T/Fastp_files_all/Everything_on_plus_strand/Merged_fastq_files/Coverage_between_4_and_1k/Base_calls_to_unique_strings/Lc_Ac_Gg_Hs_Ev_Input/Mutation_rate_raw_files/Mutation"
output_dir <- file.path(input_dir, "C_to_T_extracted")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Get all .tsv.gz files
files <- list.files(input_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)

for (f in files) {
  # Read the gzipped TSV (assuming tab-separated, no header)
  dat <- fread(cmd = paste("zcat", shQuote(f)), header = FALSE, sep = "\t", fill = TRUE)
  
  # Keep only rows with 'C -> T' in *exactly* the 10th column
  c2t_rows <- dat[V10 == "C -> T"]
  
  # Write filtered TSV
  if (nrow(c2t_rows) > 0) {
    output_base <- sub("\\.tsv\\.gz$", "_C_to_T.tsv", basename(f))
    output_path <- file.path(output_dir, output_base)
    fwrite(c2t_rows, output_path, sep = "\t", col.names = FALSE)
    
    # Create BED file (use V1, V2, V3=V2+1)
    bed <- c2t_rows[, .(V1, as.integer(V2), as.integer(V2) + 1)]
    colnames(bed) <- c("chrom", "start", "end")
    
    # Sort BED and write
    bed <- bed[order(chrom, start)]
    bed_path <- file.path(output_dir, sub("\\.tsv\\.gz$", ".bed", basename(f)))
    fwrite(bed, bed_path, sep = "\t", col.names = FALSE)
  }
}
cat("✅ Extraction and BED creation complete!\n")
