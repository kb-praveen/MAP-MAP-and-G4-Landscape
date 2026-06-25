# 1. LOAD LIBRARIES
suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr) 
})

# --- SETUP ---
master_file <- "/Users/praveenbishnoi/Desktop/G4_Final_Analysis_Complete/GC_bias_analysis_Final/Master_Retention_Analysis_Full.tsv"
output_dir  <- "/Users/praveenbishnoi/Desktop/G4_Final_Analysis_Complete/GC_bias_analysis_Final/Plots_Ranked_Shift/New"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 2. DATA LOADING & CLEANING ---
df <- read_tsv(master_file, show_col_types = FALSE)

g4_colors <- c(
  "Highly Enriched"     = "#ca0020", 
  "Moderately Enriched" = "#f4a582", 
  "Moderately Reduced"  = "#92c5de", 
  "Highly Reduced"      = "#0571b0"
)

enrichment_levels <- c("Highly Enriched", "Moderately Enriched", "Moderately Reduced", "Highly Reduced")

df_clean <- df %>%
  mutate(Enrichment_Tier = case_when(
    grepl("Highly_Enriched", Enrichment_Tier) ~ "Highly Enriched",
    grepl("Moderately_Enriched", Enrichment_Tier) ~ "Moderately Enriched",
    grepl("Moderately_Reduced", Enrichment_Tier) ~ "Moderately Reduced",
    grepl("Highly_Reduced", Enrichment_Tier) ~ "Highly Reduced",
    TRUE ~ Enrichment_Tier
  )) %>%
  mutate(Enrichment_Tier = factor(Enrichment_Tier, levels = enrichment_levels))

# --- 3. PLOTTING FUNCTION (ULTRA-SCALE RANKED BARS) ---
create_global_top50_plot <- function(full_df, metric_col, title, sub_expr, filename) {
  
  top_data <- full_df %>%
    arrange(desc(!!sym(metric_col))) %>%
    head(50) %>%
    mutate(Motif_Name = str_wrap(Motif_Name, width = 30))
  
  p <- ggplot(top_data, aes(x = reorder(Motif_Name, !!sym(metric_col)), 
                            y = !!sym(metric_col), 
                            fill = Enrichment_Tier)) +
    geom_vline(xintercept = 0, color = "gray30", linetype = "solid") +
    geom_col(color = "white", linewidth = 0.2, alpha = 0.9) +
    coord_flip() +
    facet_grid(Region ~ ., scales = "free_y", space = "free_y") + 
    scale_fill_manual(values = g4_colors, name = "Enrichment Status", drop = FALSE) + 
    theme_bw() + 
    labs(title = title, subtitle = sub_expr, x = "Motifs", y = "GC-retention magnitude") +
    theme(
      legend.position = "bottom",
      # --- FONT SIZE 42/40 LOGIC ---
      plot.title = element_text(size = 42, margin = margin(b=25)),
      plot.subtitle = element_text(size = 40, color = "black", margin = margin(b=20)),
      strip.text.y = element_text(size = 40, angle = -90),
      axis.title = element_text(size = 40),
      legend.title = element_text(size = 40),
      legend.text = element_text(size = 35),
      axis.text.y = element_text(size = 32, color = "black", lineheight = 0.8),
      axis.text.x = element_text(size = 32, color = "black"),
      # ----------------------------
      strip.background = element_rect(fill = "gray95"),
      panel.grid.major.y = element_blank(),
      panel.spacing = unit(3, "lines"),
      plot.margin = margin(2, 2, 2, 2, "cm")
    )
  
  ggsave(file.path(output_dir, filename), p, width = 26, height = 32, dpi = 300)
}

# --- 4. DEFINE UPDATED SUBTITLE FORMULAS ---

# Absolute Subtitle: Delta GC = pG4_homeotherm(motif mean GC) - pG4_poikilotherm(motif mean GC)
sub_abs <- bquote(Delta*GC == pG4[homeotherm]~(motif~mean~GC) - pG4[poikilotherm]~(motif~mean~GC))

# Normalized Subtitle: Delta GC_(homeo-poikilo) = pG4(motif mean GC) - pG4-free(motif mean GC)
sub_norm <- bquote(Delta*GC[homeo-poikilo] == pG4~(motif~mean~GC) - "pG4-free"~(motif~mean~GC))

# --- 5. EXECUTE ---

create_global_top50_plot(df_clean, "Absolute_G4_Shift", 
                         "Absolute GC-retention pG4s TFBS", 
                         sub_abs, "Global_Top50_Absolute_FINAL_UpdatedSub.png")

create_global_top50_plot(df_clean, "Delta_Delta_Shift", 
                         "Background-normalized GC-retention in pG4s TFBS", 
                         sub_norm, "Global_Top50_Background_Normalized_FINAL_UpdatedSub.png")

cat("\n[SUCCESS] Ranked bar plots with updated subtitles saved to:", output_dir)