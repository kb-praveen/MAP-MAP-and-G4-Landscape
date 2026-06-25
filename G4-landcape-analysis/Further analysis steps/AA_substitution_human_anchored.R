library(ape)
library(phangorn)

# 1. Read the FASTA alignment natively using phangorn's primary importer
align_raw <- read.phyDat("/Users/praveenbishnoi/Downloads/CGGBP1_alignment.fasta", format = "fasta", type = "AA")

# 2. Compute the amino acid substitution distance matrix using the JTT model
dist_matrix <- as.matrix(dist.ml(align_raw, model = "JTT"))

# 3. Extract divergence relative to your reference anchor
reference_species <- "Homo_sapiens" 

cggbp1_metrics <- data.frame(
  TreeNames = rownames(dist_matrix),
  CGGBP1_divergence = dist_matrix[, reference_species],
  stringsAsFactors = FALSE
)

# 4. Save this intermediate file as a clean tab-separated spreadsheet
write.table(cggbp1_metrics, 
            "/Users/praveenbishnoi/Downloads/CGGBP1_alignment_aa_processed_rates.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

cat("\n*** Distance matrix generated successfully! ***\n")