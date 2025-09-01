library(data.table)
library(R.utils)

# --------- USER EDIT: Specify paths here ---------
input_dir <- "Path_to_mutation_file_directory/Mutation"
output_dir <- "Path_to_output_file_directory/A_to_G"      # <<<--- EDIT THIS PATH
# -------------------------------------------------

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Get all .tsv.gz files
files <- list.files(input_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)

for (f in files) {
  # Read the gzipped TSV (tab-separated, no header)
  dat <- fread(cmd = paste("zcat", shQuote(f)), header = FALSE, sep = "\t", fill = TRUE)
  
  # Keep only rows with 'A -> G' in exactly the 10th column
  c2t_rows <- dat[V10 == "A -> G"]
  
  # Write filtered TSV
  if (nrow(c2t_rows) > 0) {
    output_base <- sub("\\.tsv\\.gz$", "_A_to_G.tsv", basename(f))
    output_path <- file.path(output_dir, output_base)
    fwrite(c2t_rows, output_path, sep = "\t", col.names = FALSE)
    
    # Create BED file (use V1, V2, V3 = V2 + 1)
    bed <- c2t_rows[, .(V1, as.integer(V2), as.integer(V2) + 1)]
    colnames(bed) <- c("chrom", "start", "end")
    
    # Sort BED and write
    bed <- bed[order(chrom, start)]
    bed_path <- file.path(output_dir, sub("\\.tsv\\.gz$", ".bed", basename(f)))
    fwrite(bed, bed_path, sep = "\t", col.names = FALSE)
  }
}
