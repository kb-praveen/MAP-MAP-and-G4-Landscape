# --- 1. LOAD LIBRARIES ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

# --- SETUP PATHS ---
raw_species_file <- "Path/Motif_GC_content_summary.xlsx"
name_map_file    <- "Path/Checking/Motif_ID_and_name.xlsx"
output_dir       <- "Path/Checking/GC_retention_analysis/Individual_Species_Plots"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- CONFIGURABLE PLOT SIZING INTERFACE ---
plot_width_inches  <- 12.0   
plot_height_inches <- 14.0  
plot_dpi           <- 300   # Hardcoded sharp publication resolution standard
# --------------------------------------------------------

# --- 2. LOAD DATA & CLEAN STRUCTURE ---
if (!file.exists(name_map_file) | !file.exists(raw_species_file)) {
  stop("[ERROR] One or more input files are missing from your Desktop directory.")
}

names_df   <- read_excel(name_map_file)
id_to_name <- deframe(names_df %>% dplyr::select(Motif_ID, Motif_Name))
raw_data   <- read_excel(raw_species_file)

# Standardize names and filter exclusively for the requested target motifs
gc_clean <- raw_data %>%
  mutate(
    Motif_Name = id_to_name[Motif_ID],
    Thermal_Group = factor(Thermal_Group, levels = c("Poikilotherm", "Homeotherm"))
  ) %>%
  filter(toupper(Motif_Name) %in% c("CGGBP1", "NFIX")) %>%
  # Determine structural Context Type (Genome-wide vs Promoters) and Region Layer Type (pG4 vs free)
  mutate(
    Context_Type = if_else(str_detect(Category, regex("Promoter", ignore_case = TRUE)), "Promoters", "Genome-wide"),
    Region_Type  = if_else(str_detect(Category, "_G4"), "pG4-forming regions", "pG4-free regions"),
    Context_Type = factor(Context_Type, levels = c("Genome-wide", "Promoters")),
    Region_Type  = factor(Region_Type, levels = c("pG4-forming regions", "pG4-free regions"))
  )

# --- 3. EXHAUSTIVE SEPARATED AUTOMATIC RANGE PLOTTING PIPELINE ---
plot_comprehensive_species_distributions <- function(df) {
  
  unique_contexts <- unique(df$Context_Type)
  all_stats_list  <- list()
  
  for (current_context in unique_contexts) {
    context_df <- df %>% filter(Context_Type == current_context)
    unique_regions <- unique(context_df$Region_Type)
    
    for (current_region in unique_regions) {
      region_df <- context_df %>% filter(Region_Type == current_region)
      unique_motifs <- unique(region_df$Motif_Name)
      
      for (current_motif in unique_motifs) {
        motif_data <- region_df %>% filter(Motif_Name == current_motif)
        
        homeo_vec <- motif_data$Mean_GC[motif_data$Thermal_Group == "Homeotherm"] %>% na.omit()
        poiki_vec <- motif_data$Mean_GC[motif_data$Thermal_Group == "Poikilotherm"] %>% na.omit()
        
        n_h <- length(homeo_vec)
        n_p <- length(poiki_vec)
        
        if (n_h < 2 | n_p < 2) next
        
        w_test <- suppressWarnings(wilcox.test(x = homeo_vec, y = poiki_vec, alternative = "two.sided"))
        p_val  <- w_test$p.value
        
        y_max <- max(motif_data$Mean_GC, na.rm = TRUE)
        y_min <- min(motif_data$Mean_GC, na.rm = TRUE)
        y_pos <- y_max + (abs(y_max - y_min) * 0.05)
        
        p_label <- if_else(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val))
        
        stats_row <- tibble(
          Motif_Name = current_motif,
          Context_Type = as.character(current_context),
          Region_Type = as.character(current_region),
          Homeotherm_Species_Count = n_h,
          Poikilotherm_Species_Count = n_p,
          Homeotherm_Mean_GC = mean(homeo_vec, na.rm = TRUE),
          Poikilotherm_Mean_GC = mean(poiki_vec, na.rm = TRUE),
          Unpaired_Wilcox_P = p_val
        )
        all_stats_list[[length(all_stats_list) + 1]] <- stats_row
        
        thermal_palette <- c("Poikilotherm" = "#8D85A1", "Homeotherm" = "#8EBEBF")
        
        # Dynamically constructed the combined title strings
        constructed_title <- paste(current_context, current_region)
        constructed_subtitle <- paste(toupper(current_motif), "motif")
        
        p <- ggplot(motif_data, aes(x = Thermal_Group, y = Mean_GC)) +
          theme_classic(base_size = 15) +
          theme(
            panel.grid.major.y = element_line(color = "gray93", linewidth = 0.4),
            panel.grid.major.x = element_blank(),
            panel.grid.minor   = element_blank(),
            
            # FIXED: Updated hjust to 0.5 to center-align both title and subtitle panels perfectly
            plot.title    = element_text(size = 42, face = "plain", hjust = 0.5, margin = margin(b=15)),
            plot.subtitle = element_text(size = 42, face = "plain", hjust = 0.5, margin = margin(b=25)),
            axis.title.x  = element_blank(), 
            axis.title.y  = element_text(size = 42, face = "plain", margin = margin(r=24)),
            axis.text.x   = element_text(color = "black", size = 42, face = "plain", margin = margin(t=12)),
            axis.text.y   = element_text(color = "black", size = 42, face = "plain"),
            
            legend.position = "none",
            plot.margin      = margin(1, 1, 1, 1, "cm")
          ) +
          geom_boxplot(
            aes(fill = Thermal_Group),
            outlier.shape = NA,
            color = "black",
            linewidth = 0.6,
            width = 0.35,
            alpha = 0.15,
            fill = "gray80"
          ) +
          geom_jitter(
            aes(fill = Thermal_Group),
            color = "black",
            shape = 21,
            stroke = 1.0,
            width = 0.10,
            height = 0,
            size = 7.5,
            alpha = 0.85
          ) +
          
          annotate(
            "text",
            x = 1.5, y = y_pos,
            label = p_label,
            color = "black",
            size = 13,
            vjust = 0
          ) +
          
          scale_color_manual(values = thermal_palette) +
          scale_fill_manual(values = thermal_palette) +
          
          scale_x_discrete(labels = c(
            "Poikilotherm" = sprintf("Poikilotherms\n(n=%d)", n_p),
            "Homeotherm"   = sprintf("Homeotherms\n(n=%d)", n_h)
          )) +
          
          labs(
            title    = constructed_title,                                   
            subtitle = constructed_subtitle, 
            y        = "Mean motif GC% of individual species"
          ) +
          scale_y_continuous(expand = expansion(mult = c(0.06, 0.15)))
        
        sanitized_context <- gsub("-", "_", gsub(" ", "_", current_context))
        sanitized_region  <- gsub("-", "_", gsub(" ", "_", current_region))
        
        ggsave(
          filename = file.path(output_dir, paste0("Plot_", sanitized_context, "_Species_", toupper(current_motif), "_", sanitized_region, ".png")),
          plot     = p, 
          width    = plot_width_inches, 
          height   = plot_height_inches, 
          dpi      = plot_dpi,
          units    = "in"
        )
      }
    }
  }
  
  combined_stats_df <- bind_rows(all_stats_list)
  write_tsv(combined_stats_df, file.path(output_dir, "Species_Stats_Combined_All_Contexts.tsv"))
}

# --- 4. EXECUTE PIPELINE ---
cat("[INFO] Running graphics execution engine pass with center-aligned headers...\n")
plot_comprehensive_species_distributions(gc_clean)

cat(sprintf("\n🚀 SUCCESS! Plots rendered with centered labels inside:\n%s\n", output_dir))
