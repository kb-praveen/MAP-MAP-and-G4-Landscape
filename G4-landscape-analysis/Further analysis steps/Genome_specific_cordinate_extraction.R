#!/usr/bin/env Rscript

#This script was applied to several clades separately
# Automated GTF → BED6 extractor for multiple species (handles .gtf.gz)
# Input: taxa_dir/*gtf* → species-specific BED6 files
# Output: Genome_specific_coordinates/Mammals/species_name/*.bed

library(dplyr)

# === CONFIGURE PATHS ===
taxa_dir <- " path to the directory /Gtf_files/Mammals"
output_base <- "path to the directory /Genome_specific_coordinates"

# Create output base directory
dir.create(output_base, showWarnings = FALSE, recursive = TRUE)

# Find GTF files (matches .gtf, .gtf.gz, _gtf.gz, etc.)
gtf_files <- Sys.glob(file.path(taxa_dir, "*gtf*"))
cat("Found", length(gtf_files), "GTF files:\n")
print(gtf_files)
cat("\n")

if (length(gtf_files) == 0) {
  stop("No GTF files found in ", taxa_dir)
}

# === HELPER FUNCTIONS ===
parse_key <- function(attr_str, key) {
  if (is.na(attr_str) || attr_str == "." || nchar(attr_str) < 10) return(NA)
  pat <- paste0(key, ' "([^"]+)"')
  match <- regmatches(attr_str, regexpr(pat, attr_str, perl = TRUE))
  if (length(match) > 0) gsub(pat, "\\1", match) else NA
}

read_gtf_safe <- function(gtf_file) {
  # Read all lines, filter to those with exactly 9 tab-separated fields (skip headers/comments)
  lines <- readLines(gtf_file)
  valid_lines <- lines[sapply(strsplit(lines, "\t"), length) == 9]
  if (length(valid_lines) == 0) stop("No valid GTF data lines found in ", basename(gtf_file))
  
  # Parse valid lines into data.frame
  con <- textConnection(valid_lines)
  gtf_raw <- read.table(con, sep = "\t", quote = "", stringsAsFactors = FALSE,
                        col.names = c("seqnames", "source", "type", "start", "end", 
                                      "score", "strand", "phase", "attributes"))
  close(con)
  return(gtf_raw)
}

make_bed6 <- function(df) {
  data.frame(
    chrom = df$seqnames,
    start = df$start - 1,  # BED 0-based
    end = df$end,
    name = ifelse(!is.na(df$gene_id), df$gene_id, 
                  ifelse(!is.na(df$transcript_id), df$transcript_id, "NA")),
    score = ".",
    strand = df$strand,
    stringsAsFactors = FALSE
  )
}

# === PROCESS EACH GTF ===
for (gtf_file in gtf_files) {
  cat("=== Processing", basename(gtf_file), "===\n")
  
  # Extract species name (remove .gtf.gz or _gtf.gz)
  species_name <- sub(".*([A-Z][a-z]+_[a-z]+).*gtf.*", "\\1", basename(gtf_file))
  species_dir <- file.path(output_base, basename(dirname(gtf_file)), species_name)
  dir.create(species_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("Species dir:", species_dir, "\n")
  
  # Read GTF safely (handles .gz, headers, comments)
  gtf_raw <- read_gtf_safe(gtf_file)
  
  cat("Total lines:", nrow(gtf_raw), "\n")
  
  # Parse gene/transcript IDs
  gtf_raw$gene_id <- sapply(gtf_raw$attributes, parse_key, "gene_id")
  gtf_raw$transcript_id <- sapply(gtf_raw$attributes, parse_key, "transcript_id")
  
  ## 1. TSS coordinates
  transcripts <- gtf_raw[gtf_raw$type == "transcript", ]
  if (nrow(transcripts) == 0) transcripts <- gtf_raw[gtf_raw$type == "gene", ]
  tss_bed <- transcripts
  tss_bed$start <- ifelse(tss_bed$strand == "+", tss_bed$start, tss_bed$end)
  tss_bed$end <- tss_bed$start
  tss_bed6 <- make_bed6(tss_bed)
  write.table(tss_bed6, file.path(species_dir, "TSS_coordinates.bed"), sep = "\t", 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  ## 2. Gene Body coordinates
  genes <- gtf_raw[gtf_raw$type == "gene", ]
  if (nrow(genes) > 0) {
    genes_bed6 <- make_bed6(genes)
    write.table(genes_bed6, file.path(species_dir, "Gene_Body_coordinates.bed"), sep = "\t", 
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  ## 3. 5' UTR regions
  five_utr <- gtf_raw[gtf_raw$type %in% c("five_prime_utr", "5UTR"), ]
  if (nrow(five_utr) > 0) {
    five_bed6 <- make_bed6(five_utr)
    write.table(five_bed6, file.path(species_dir, "5UTR_regions.bed"), sep = "\t", 
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  ## 4. 3' UTR regions
  three_utr <- gtf_raw[gtf_raw$type %in% c("three_prime_utr", "3UTR"), ]
  if (nrow(three_utr) > 0) {
    three_bed6 <- make_bed6(three_utr)
    write.table(three_bed6, file.path(species_dir, "3UTR_regions.bed"), sep = "\t", 
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  ## 5. Intron coordinates (derived from exons)
  exons <- gtf_raw[gtf_raw$type == "exon" & !is.na(gtf_raw$transcript_id), ]
  if (nrow(exons) > 0) {
    exons_by_tx <- split(exons, exons$transcript_id)
    introns_rows <- do.call(rbind, lapply(names(exons_by_tx), function(tx) {
      tx_exons <- exons_by_tx[[tx]][order(exons_by_tx[[tx]]$start), ]
      if (nrow(tx_exons) > 1) {
        starts <- tx_exons$end[-nrow(tx_exons)] + 1
        ends <- tx_exons$start[-1] - 1
        valid <- ends > starts
        if (any(valid)) {
          data.frame(seqnames = tx_exons$seqnames[1], start = starts[valid], 
                     end = ends[valid], strand = tx_exons$strand[1], 
                     gene_id = tx_exons$gene_id[1], transcript_id = tx)
        } else data.frame()
      } else data.frame()
    }))
    if (nrow(introns_rows) > 0) {
      introns_bed6 <- make_bed6(introns_rows)
      write.table(introns_bed6, file.path(species_dir, "Intron_coordinates.bed"), sep = "\t", 
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
  }
  
  cat("✅ Outputs saved to:", species_dir, "\n")
}

cat("\n🎉 ALL SPECIES PROCESSED!\n")
