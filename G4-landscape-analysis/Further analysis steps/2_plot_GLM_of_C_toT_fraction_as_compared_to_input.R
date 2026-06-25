# 1. LOAD REQUIRED LIBRARIES
suppressPackageStartupMessages({
  library(openxlsx)
  library(tidyverse)
  library(ggplot2)
})

# --- SETUP WORKING DIRECTORIES ---
master_dir <- "Path/Promoter"
excel_file <- file.path(master_dir, "Final_CtoT_GLM_Regression_Significance.xlsx")
plot_file  <- file.path(master_dir, "CtoT_Global_Significance_Bar_Plot.png")

# --- CONFIGURABLE PLOT SIZING INTERFACE ---
plot_width_inches  <- 16   
plot_height_inches <- 8  
plot_dpi           <- 300   # Sharp publication resolution standard
# --------------------------------------------------------

if (!file.exists(excel_file)) {
  stop(sprintf("[ERROR] Pre-calculated regression statistics workbook missing at:\n%s\nRun your GLM engine script first.", excel_file))
}

# --- 2. IMPORT METRIC DATA & ATTACH SPECIFIC CUSTOM COLORS ---
cat("[INFO] Loading statistical metrics from Excel summary...\n")
plot_data <- read.xlsx(excel_file, sheet = "Global_Regression_Significance")

custom_colors <- c(
  "Lc" = "#A5A3C7",
  "Ev" = "#B4B1BD",
  "Hs" = "#CF3A3A",
  "Gg" = "#4681AB",
  "Ac" = "#629E6F"
)

plot_data <- plot_data %>%
  mutate(
    Log2_FC      = as.numeric(Log2_Fold_Change),
    Log2_Error   = as.numeric(Log2_SE),
    FDR_Value    = as.numeric(FDR),
    
    # Formatted string to match plotmath expression rules for proper P(adj) subscripting
    FDR_Text     = paste0("P[adj] == '", formatC(FDR_Value, format = "e", digits = 2), "'"),
    
    # Calculate exact 95% Confidence Intervals mapped to the Log2 scale
    Lower_Bound  = Log2_FC - (1.96 * Log2_Error),
    Upper_Bound  = Log2_FC + (1.96 * Log2_Error),
    
    # Explicit label coordinate pinning relative to the dashed line (0)
    Text_Position = case_when(
      Sample %in% c("Hs", "Gg", "Ev") ~ 0.03,  
      Sample %in% c("Lc", "Ac")       ~ -0.03, 
      TRUE                            ~ 0.03
    ),
    
    # Inverted horizontal alignment to point cleanly outward from the center
    Text_Hjust = case_when(
      Sample %in% c("Hs", "Gg", "Ev") ~ 0, # Left-aligned (extends rightward)
      Sample %in% c("Lc", "Ac")       ~ 1, # Right-aligned (extends leftward)
      TRUE                            ~ 0
    )
  )

# --- 3. GENERATE THE COLOR-CODED BAR PLOT ---
cat("[INFO] Generating bar plot structure with synchronized text matrix...\n")

plot_obj <- ggplot(plot_data, aes(x = Log2_FC, y = reorder(Sample, Log2_FC), fill = Sample)) +
  # Draw solid horizontal bar elements spanning cleanly from the 0 reference line
  geom_col(width = 0.6, alpha = 0.9, color = NA) +
  
  # Draw the 95% Confidence Interval error bar whiskers horizontally
  geom_errorbarh(aes(xmin = Lower_Bound, xmax = Upper_Bound), 
                 height = 0.15, size = 0.8, color = "black") +
  
  # Reference line at Log2FC = 0 
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.6) +
  
  # Text parameters to clear the axes bounds smoothly
  geom_text(aes(x = Text_Position, label = FDR_Text, hjust = Text_Hjust), 
            vjust = 0.5, 
            size = 8, 
            parse = TRUE, 
            show.legend = FALSE, 
            color = "black") +
  
  # Force the application of manual custom color hex map vector
  scale_fill_manual(values = custom_colors) +
  
  # Inject padding buffers so up-sized labels do not stretch outside plot margins
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  
  # Overrides canvas boundaries to prevent text from being truncated
  coord_cartesian(clip = "off") +
  
  # Transitioned to clean minimal theme profile layout
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray95", linewidth = 0.4),
    
    # MODIFICATION: Synchronized and expanded font sizes cleanly across labels
    plot.title    = element_blank(), # MODIFICATION: Completely dropped title element mapping
    axis.title.x  = element_text(size = 26, face = "plain", hjust = 0.5, margin = margin(t=24), lineheight = 1.2), # MODIFICATION: Enforced center-alignment and lineheight wrapping
    axis.title.y  = element_text(size = 26, face = "plain", margin = margin(r=24)),
    axis.text.x   = element_text(color = "black", size = 26, face = "plain"),
    axis.text.y   = element_text(color = "black", size = 26, face = "plain"),
    
    legend.position = "none",
    plot.margin     = margin(1, 1, 1, 1, "cm")
  ) +
  labs(
    y = "MeDIP samples",
    # MODIFICATION: Configured exact mathematical plotmath expression schema for subscript and arrow rendering
    x = bquote(atop(log[2] ~ "fold change C" %->% "T transition rates", "(sample versus input)"))
  )

# --- 4. EXPORT PLOT VIA GGSAVE ---
cat(sprintf("[INFO] Exporting polished image file (Size: %.2fin x %.2fin @ %d DPI)...\n", 
            plot_width_inches, plot_height_inches, plot_dpi))

ggsave(
  filename = plot_file,
  plot     = plot_obj,
  width    = plot_width_inches,
  height   = plot_height_inches,
  dpi      = plot_dpi,
  units    = "in"
)

cat(sprintf("\n🚀 SUCCESS! Publication-ready bar plot saved cleanly to:\n%s\n", plot_file))