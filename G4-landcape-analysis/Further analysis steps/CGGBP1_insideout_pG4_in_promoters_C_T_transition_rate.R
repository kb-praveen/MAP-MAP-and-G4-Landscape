# ==========================================
# 1. LOAD LIBRARIES
# ==========================================
suppressPackageStartupMessages({
  library(openxlsx)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(ggpattern) 
})

# ==========================================
# 2. SETUP DATA & PATHWAYS
# ==========================================
file_path  <- "/Users/praveenbishnoi/Desktop/G4/G4_Final_Analysis_Complete/C_to_T_transition_in_pG4s_and_promoters.xlsx"
output_fig <- "/Users/praveenbishnoi/Desktop/G4/G4_Final_Analysis_Complete/C_to_T_transition_plot.png"

df_raw <- read.xlsx(file_path, sheet = 1)

# ==========================================
# 3. TIDY TRANSFORMATION LAYER
# ==========================================
df_tidy <- df_raw %>%
  pivot_longer(
    cols = everything(),
    names_to = "Raw_Header",
    values_to = "Transition_Rate"
  ) %>%
  mutate(
    Sample = case_when(
      startsWith(Raw_Header, "Lc") ~ "Lc",
      startsWith(Raw_Header, "Ac") ~ "Ac",
      startsWith(Raw_Header, "Gg") ~ "Gg",
      startsWith(Raw_Header, "Hs") ~ "Hs",
      startsWith(Raw_Header, "Ev") ~ "Ev"
    ),
    Genomic_Feature = case_when(
      grepl("pG4", Raw_Header) ~ "pG4 in promoter",
      TRUE                     ~ "Entire promoter"
    )
  )

df_tidy$Sample          <- factor(df_tidy$Sample, levels = c("Lc", "Ac", "Gg", "Hs", "Ev"))
df_tidy$Genomic_Feature <- factor(df_tidy$Genomic_Feature, levels = c("Entire promoter", "pG4 in promoter"))

# ==========================================
# 4. PALETTE DEFINITIONS
# ==========================================
custom_colors <- c(
  "Lc" = "#A5A3C7",
  "Ev" = "#B4B1BD",
  "Hs" = "#CF3A3A",
  "Gg" = "#4681AB",
  "Ac" = "#629E6F"
)

# ==========================================
# 5. CHART RENDER ENGINE (GGPLOT2 + GGPATTERN)
# ==========================================
p_bar <- ggplot(df_tidy, aes(x = Sample, y = Transition_Rate, fill = Sample)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  
  geom_col_pattern(
    aes(
      pattern = Genomic_Feature,
      pattern_spacing = Genomic_Feature
    ),
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    linewidth = 0.5,
    pattern_fill = "black",       
    pattern_color = "transparent", 
    pattern_density = 0.15,       
    pattern_angle = 0             
  ) +
  
  scale_pattern_manual(
    values = c("Entire promoter" = "none", "pG4 in promoter" = "circle"),
    name = "Genomic Feature: "
  ) +
  scale_pattern_spacing_manual(
    values = c("Entire promoter" = 0, "pG4 in promoter" = 0.025),
    guide = "none"
  ) +
  
  scale_fill_manual(values = custom_colors, name = "Sample: ") +
  scale_y_continuous(
    limits = c(-0.025, 0.04), 
    breaks = seq(-0.02, 0.04, by = 0.01),
    labels = sprintf("%.2f", seq(-0.02, 0.04, by = 0.01))
  ) +
  theme_classic(base_size = 18) +
  theme(
    text = element_text(family = "sans", face = "plain", color = "black"),
    axis.title.x = element_text(size = 16, margin = margin(t = 12)),
    axis.title.y = element_text(size = 16, margin = margin(r = 12)),
    axis.text.x  = element_text(size = 15, color = "black"),
    axis.text.y  = element_text(size = 14, color = "black"),
    axis.line.x  = element_blank(), 
    axis.line.y  = element_line(linewidth = 0.7, color = "black"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.7, color = "black"),
    panel.grid.major.y = element_line(color = "#EAEAEA", linewidth = 0.5),
    
    # CHANGED: Font sizes expanded across all legend parameters
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 16, face = "plain"),
    legend.text = element_text(size = 15, face = "plain"),
    legend.key.size = unit(1.2, "cm"), 
    legend.margin = margin(t = 12)
  ) +
  labs(
    x = "CGGBP1 evolutionary variants",
    y = "C-to-T transition rate (mean & input corrected)"
  ) +
  guides(
    fill = guide_legend(override.aes = list(pattern = "none")),
    pattern = guide_legend(override.aes = list(fill = "grey80", pattern_spacing = 0.03))
  )

# ==========================================
# 6. FILE EXPORT LAYER
# ==========================================
ggsave(output_fig, p_bar, width = 14.0, height = 6.0, dpi = 300, bg = "white")
cat(paste0("\n*** Success! Plot with larger legend metrics written to:\n", output_fig, "\n"))