library(dplyr)
library(tidyr)
library(ggplot2)

# --- 1. SETUP ---
base_dir <- "/scratch/kumarpraveen/G4_motif_analysis"
results_dir <- file.path(base_dir, "Results/Homeotherm_vs_Poikilotherm")
out_dir <- file.path(results_dir, "GC_content_analysis")
df_path <- file.path(out_dir, "Motif_Delta_GC_Per_Species.tsv")

if (!file.exists(df_path)) stop("Missing input file! Ensure Python script ran successfully.")
df_wide <- read.delim(df_path)

# --- 2. DATA CLEANING & AGGREGATION (MEAN-BASED) ---
motif_summary <- df_wide %>%
  group_by(Motif_ID, Enrichment_Tier, Thermal_Group) %>%
  summarise(
    # Using Mean to capture the full evolutionary shift across species
    P_Delta = mean(Delta_Promoter, na.rm = TRUE),
    G_Delta = mean(Delta_Genome, na.rm = TRUE),
    .groups = 'drop'
  )

# --- 3. CREATE PAIRED DATA (HOMEOTHERM VS POIKILOTHERM) ---
motif_paired <- motif_summary %>%
  pivot_longer(cols = c(P_Delta, G_Delta), 
               names_to = "Region", values_to = "Delta_Val") %>%
  mutate(Region = ifelse(Region == "P_Delta", "Promoter", "Genome")) %>%
  pivot_wider(names_from = Thermal_Group, values_from = Delta_Val) %>%
  filter(!is.na(Homeotherm) & !is.na(Poikilotherm))

# Set Tier order
tier_levels <- c("Highly_Enriched_Homeotherm", "Moderately_Enriched_Homeotherm", 
                 "Moderately_Reduced_Homeotherm", "Highly_Reduced_Homeotherm")
motif_paired$Enrichment_Tier <- factor(motif_paired$Enrichment_Tier, levels = tier_levels)

# --- 4. GENERATE MOTIF IDENTITY PLOT ---
plot <- ggplot(motif_paired, aes(x = Poikilotherm, y = Homeotherm)) +
  # Diagonal reference line (y=x) indicating no change between groups
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  # Zero reference lines to show background equilibrium
  geom_vline(xintercept = 0, color = "gray85") +
  geom_hline(yintercept = 0, color = "gray85") +
  # Points for individual motifs - No Labels as requested
  geom_point(aes(color = Enrichment_Tier), size = 3, alpha = 0.6) +
  # Facet by Region (Promoter vs Genome)
  facet_wrap(~Region) +
  # Styling
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 14) +
  labs(
    title = "Evolutionary Trajectory of Motif ΔGC (Mean-based)",
    subtitle = "Comparing individual motif GC-specialization: Homeotherm vs Poikilotherm",
    x = "Poikilotherm Mean ΔGC (%)",
    y = "Homeotherm Mean ΔGC (%)",
    color = "Motif Category"
  ) +
  theme(
    legend.position = "bottom", 
    legend.direction = "vertical",
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold")
  )

# --- 5. SAVE ---
output_file <- file.path(out_dir, "Motif_Identity_Mean_Scatter_Plot_NoLabels.png")
ggsave(output_file, plot, width = 12, height = 10, dpi = 300)

cat(paste0("Success: Clean Mean-based Motif Identity plot saved to ", output_file, "\n"))