input_file <- "Path_to_input_bed_file/Homo_sapiens.bed"
output_file <- "Path_to_input_bed_file/Homo_sapiens_G4.bed"

# Read the file line by line
lines <- readLines(input_file)

# Open a connection to write output line by line
con_out <- file(output_file, open = "w")

for (line in lines) {
  fields <- strsplit(line, "\t")[[1]]
  if (length(fields) < 3) next  # Skip invalid lines

  # Extract only the first token from the first column (up to first space)
  chr_raw <- strsplit(fields[1], " ")[[1]][1]

  # Format chromosome
  if (grepl("^[0-9]+$", chr_raw)) {
    chr <- paste0("chr", chr_raw)
  } else if (chr_raw %in% c("X", "Y", "M")) {
    chr <- paste0("chr", chr_raw)
  } else {
    chr <- chr_raw
  }

  # Extract start and end
  start <- as.integer(fields[2])
  end   <- as.integer(fields[3])
  if (is.na(start) || is.na(end)) next

  # Write as BED line
  cat(chr, start, end, sep = "\t", file = con_out, append = TRUE)
  cat("\n", file = con_out, append = TRUE)

  # Optional: also print each line to screen
  cat(paste(chr, start, end, sep = "\t"), "\n")
}

close(con_out)

cat(paste("✅ BED file written to:", output_file, "\n"))
