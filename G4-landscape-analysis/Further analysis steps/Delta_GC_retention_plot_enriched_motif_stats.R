# 1. LOAD LIBRARIES
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

# --- SETUP PATHS ---
gc_summary_file <- "Path/Motif_GC_content_summary.xlsx"
name_map_file   <- "Path/Checking/Motif_ID_and_name.xlsx"
output_dir      <- "Path/Checking/GC_retention_analysis"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- CONFIGURATION (TARGET DIRECTORY SIZING MATRIX) ---
plot_width_inches  <- 14   
plot_height_inches <- 20  
plot_dpi           <- 300   # Hardcoded sharp publication resolution standard
# --------------------------------------------------------

# --- 2. LOAD DATA & CLEAN STRUCTURE ---
if (!file.exists(name_map_file) | !file.exists(gc_summary_file)) {
  stop("[ERROR] One or more input files are missing from your Desktop directory.")
}

names_df   <- read_excel(name_map_file)
id_to_name <- deframe(names_df %>% dplyr::select(Motif_ID, Motif_Name))
gc_data    <- read_excel(gc_summary_file)

# Standardize names and fix sorting factors for the plots
gc_clean <- gc_data %>%
  mutate(
    Motif_Name = id_to_name[Motif_ID],
    Enrichment_Tier = gsub("_", " ", Enrichment_Tier),
    Enrichment_Tier = factor(Enrichment_Tier, levels = c("Highly Enriched", "Moderately Enriched", 
                                                         "Moderately Reduced", "Highly Reduced")),
    Thermal_Group = factor(Thermal_Group, levels = c("Poikilotherm", "Homeotherm"))
  ) %>%
  filter(!is.na(Enrichment_Tier)) 

# --- 3. DELTA INDEPENDENT PLOTTING & ANNOTATION FUNCTION ---
plot_gc_deltas <- function(df, context_keyword, display_subtitle, filename_base) {
  
  # Clean and isolate structural context layers 
  context_df <- df %>% 
    filter(str_detect(Category, regex(context_keyword, ignore_case = TRUE))) %>%
    mutate(Region_Type = if_else(str_detect(Category, "_G4"), "pG4-forming regions", "pG4-free regions"),
           Region_Type = factor(Region_Type, levels = c("pG4-forming regions", "pG4-free regions")))
  
  if (nrow(context_df) == 0) return(NULL)
  
  # --- STEP A: RESHAPE MATRIX & CALC DELTAS ---
  delta_data <- context_df %>%
    dplyr::select(Motif_ID, Motif_Name, Enrichment_Tier, Region_Type, Thermal_Group, Mean_GC) %>%
    pivot_wider(
      id_cols = c(Motif_ID, Motif_Name, Enrichment_Tier, Region_Type),
      names_from = Thermal_Group,
      values_from = Mean_GC,
      values_fn = list(Mean_GC = ~ mean(.x, na.rm = TRUE))
    ) %>%
    filter(!is.na(Homeotherm), !is.na(Poikilotherm)) %>%
    mutate(
      Individual_GC_Delta = Homeotherm - Poikilotherm
    )
  
  # --- STEP B: PROCESS INDEPENDENT SEPARATE REGION LOOPS ---
  unique_regions <- unique(delta_data$Region_Type)
  
  for (current_region in unique_regions) {
    # Slice vector tables for standalone plot output panels
    region_delta_data <- delta_data %>% filter(Region_Type == current_region)
    
    # Calculate Two-Sided Wilcoxon significance numbers
    paired_shift_stats <- region_delta_data %>%
      group_by(Enrichment_Tier) %>%
      filter(n() >= 3) %>% 
      summarise(
        Motif_Count = n(),
        Median_GC_Delta = median(Individual_GC_Delta, na.rm = TRUE),
        Mean_GC_Delta   = mean(Individual_GC_Delta, na.rm = TRUE),
        Wilcox_Shift_P  = suppressWarnings(wilcox.test(Individual_GC_Delta, mu = 0, alternative = "two.sided")$p.value),
        
        # Position tracking line placed cleanly above highest point cluster
        Y_Position = max(Individual_GC_Delta, na.rm = TRUE) + (max(abs(Individual_GC_Delta), na.rm = TRUE) * 0.15),
        .groups = "drop"
      ) %>%
      mutate(Wilcox_Shift_Padj_BH = p.adjust(Wilcox_Shift_P, method = "BH")) %>%
      # Format labels to display numerical values or 'ns'
      mutate(Significance_Label = case_when(
        Wilcox_Shift_Padj_BH < 0.001 ~ "p < 0.001",
        Wilcox_Shift_Padj_BH < 0.05  ~ sprintf("p = %.3f", Wilcox_Shift_Padj_BH),
        TRUE                         ~ "ns"
      ))
    
    sanitized_region_name <- gsub("-", "_", gsub(" ", "_", current_region))
    
    # FIXED: Swapped out %in= for a proper %in% operator check
    background_points <- region_delta_data %>% 
      filter(!toupper(Motif_Name) %in% c("CGGBP1", "NFIX"))
    
    # FIXED: Swapped out %in= for a proper %in% operator check
    target_highlights <- region_delta_data %>%
      filter(toupper(Motif_Name) %in% c("CGGBP1", "NFIX"))
    
    # Harmonized palette values matching volcano plots exactly
    tier_palette <- c(
      "Highly Enriched"     = "#ca0020", 
      "Moderately Enriched" = "#f4a582", 
      "Highly Reduced"      = "#0571b0", 
      "Moderately Reduced"  = "#92c5de"
    )
    
    # Build complete un-faceted ggplot panel
    p <- ggplot(mapping = aes(x = Enrichment_Tier, y = Individual_GC_Delta)) +
      theme_classic(base_size = 15) +
      theme(
        panel.grid.major.y = element_line(color = "gray93", linewidth = 0.4),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        
        # Suppressed ticks below the horizontal line
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        
        # Consistent text layout properties mapped to exactly 42 points
        plot.title    = element_text(size = 44, face = "plain", hjust = 0.5, margin = margin(b=12)),
        plot.subtitle = element_text(size = 44, face = "plain", hjust = 0.5, margin = margin(b=25)),
        
        # Spacious x-axis title padding distance
        axis.title.x  = element_text(size = 44, face = "plain", margin = margin(t=45)),
        axis.title.y  = element_text(size = 44, face = "plain", margin = margin(r=24)),
        axis.text.y   = element_text(color = "black", size = 44, face = "plain"),
        
        legend.position = "none", 
        plot.margin      = margin(1, 1, 1, 1, "cm")
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
      geom_boxplot(
        data = region_delta_data,
        aes(fill = Enrichment_Tier),
        outlier.shape = NA, 
        color = "black",
        linewidth = 0.6,
        width = 0.4,
        alpha = 0.15,
        fill = "gray80" 
      ) +
      
      # Background standard data points layer
      geom_jitter(
        data = background_points,
        aes(color = Enrichment_Tier),
        width = 0.12,
        height = 0, 
        size = 8.5,
        alpha = 0.6
      ) +
      
      # Targeted highlight layer (CGGBP1 and NFIX)
      geom_jitter(
        data = target_highlights,
        aes(fill = Enrichment_Tier),
        color = "black",       
        shape = 21,            
        width = 0.12,          
        height = 0, 
        size = 8.5,            
        stroke = 3,          
        alpha = 1.0            
      ) +
      
      # Significance labels value tracking configuration
      geom_text(
        data = paired_shift_stats,
        aes(x = Enrichment_Tier, y = Y_Position, label = Significance_Label),
        color = "black",
        size = 15, 
        vjust = 0,
        inherit.aes = FALSE
      ) +
      
      scale_color_manual(values = tier_palette) +
      scale_fill_manual(values = tier_palette) +
      labs(
        title = as.character(current_region),     # Primary Title row gets Region name (e.g. pG4-forming regions)
        subtitle = display_subtitle,              # Secondary Subtitle row gets Scale context (e.g. Genome-wide)
        x = "Enrichment Categories", 
        y = bquote(plain("Individual mean motif ") * Delta * plain("GC % (Homeotherm - Poikilotherm)"))
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.08, 0.18)))
    
    # Save optimized elements safely out to disk
    ggsave(
      filename = file.path(output_dir, paste0("Plot_", filename_base, "_", sanitized_region_name, ".png")), 
      plot     = p, 
      width    = plot_width_inches, 
      height   = plot_height_inches, 
      dpi      = plot_dpi,
      units    = "in"
    )
  }
}

# --- 4. EXECUTE VISUAL GENERATION PASS ---
cat("[INFO] Running header-flipped graphics engine pipeline pass across datasets...\n")
plot_gc_deltas(gc_clean, "Genome", "Genome-wide", "Genomewide")
plot_gc_deltas(gc_clean, "Promoter", "Promoters", "Promoters")

cat(sprintf("\n🚀 SUCCESS! Script fixed. Output plots generated seamlessly inside:\n%s\n", output_dir))