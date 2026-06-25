suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(eulerr) 
  library(ggplot2)
})

# --- PATH CONFIGURATION ---
base_results <- "/scratch/kumarpraveen/G4_motif_analysis/Results"
out_dir      <- file.path(base_results, "Genome_vs_Promoter_comparison")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load your data
prom_file <- file.path(base_results, "Motif_Density_Difference_Promoters_G4_minus_Free.tsv")
gen_file  <- file.path(base_results, "Motif_Density_Difference_Genomewide_G4_minus_Free.tsv")

prom_df <- read_tsv(prom_file, show_col_types = FALSE)
gen_df  <- read_tsv(gen_file, show_col_types = FALSE)

# 1. Prepare Comprehensive Comparison Data
# We include Motif_Name here so we can save it in the text files later
full_comp <- inner_join(
  prom_df %>% 
    mutate(Mean_Prom = rowMeans(select(., -Motif_ID, -Motif_Name), na.rm = TRUE)) %>% 
    select(Motif_ID, Motif_Name, Mean_Prom),
  gen_df %>% 
    mutate(Mean_Gen = rowMeans(select(., -Motif_ID, -Motif_Name), na.rm = TRUE)) %>% 
    select(Motif_ID, Mean_Gen),
  by = "Motif_ID"
)

# 2. Define Categories for Text Files
# This ensures Motif_ID and Motif_Name are saved together
categories <- list(
  Universal_Enriched       = full_comp %>% filter(Mean_Prom > 0 & Mean_Gen > 0),
  Universal_Depleted       = full_comp %>% filter(Mean_Prom < 0 & Mean_Gen < 0),
  Promoter_Unique_Enriched = full_comp %>% filter(Mean_Prom > 0 & Mean_Gen <= 0),
  Genome_Unique_Enriched   = full_comp %>% filter(Mean_Prom <= 0 & Mean_Gen > 0),
  Promoter_Unique_Depleted = full_comp %>% filter(Mean_Prom < 0 & Mean_Gen >= 0),
  Genome_Unique_Depleted   = full_comp %>% filter(Mean_Prom >= 0 & Mean_Gen < 0)
)

# Save text files with both columns
for(name in names(categories)){
  write_tsv(categories[[name]] %>% select(Motif_ID, Motif_Name), 
            file.path(out_dir, paste0(name, ".tsv"))) # Saved as TSV for clarity
}

# --- AREA-PROPORTIONAL PLOTTING FUNCTION ---
save_euler_plot <- function(prom_ids, gen_ids, filename, main_title, fill_colors) {
  sets <- list(Promoter = prom_ids, Genome = gen_ids)
  fit <- euler(sets)
  
  png(file.path(out_dir, filename), width = 2400, height = 2400, res = 300)
  print(plot(fit,
       quantities = list(type = c("counts", "percent"), font = 1, cex = 0.7),
       fills = list(fill = fill_colors, alpha = 0.5),
       edges = list(col = "black", lwd = 1.5),
       labels = list(col = "black", font = 1, cex = 1.0),
       adjust_labels = TRUE,
       main = list(label = main_title, cex = 1.2, font = 1)))
  dev.off()
}

# --- EXECUTE PLOTS ---
# Enriched
save_euler_plot(
  prom_ids = categories$Universal_Enriched$Motif_ID, # Logic: Promoter > 0 includes Universal + Unique
  gen_ids  = c(categories$Universal_Enriched$Motif_ID, categories$Genome_Unique_Enriched$Motif_ID),
  filename = "Euler_Enriched_Clean_300dpi.png",
  main_title = "Motifs Enriched in G4s",
  fill_colors = c("#377eb8", "#e41a1c")
)

# Depleted
save_euler_plot(
  prom_ids = c(categories$Universal_Depleted$Motif_ID, categories$Promoter_Unique_Depleted$Motif_ID),
  gen_ids  = categories$Universal_Depleted$Motif_ID,
  filename = "Euler_Depleted_Clean_300dpi.png",
  main_title = "Motifs Depleted in G4s",
  fill_colors = c("#4daf4a", "#984ea3")
)

cat("\n[SUCCESS] Euler diagrams and detailed category files saved in:", out_dir, "\n")