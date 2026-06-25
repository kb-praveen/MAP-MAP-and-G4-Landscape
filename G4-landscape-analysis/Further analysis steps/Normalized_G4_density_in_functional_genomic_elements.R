# ==========================================
# 1. LOAD MASTER LIBRARIES
# ==========================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggtext)
})

# ==========================================
# 2. GLOBAL PATHS & CONFIGURATIONS
# ==========================================
base_dir   <- "Path/G4/G4_analysis/Output"
output_dir <- "Path/G4/G4_analysis/Output/New_for_fig_1"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

clade_colors <- c("Mammals"="#E41A1C", "Aves"="#377EB8", "Reptiles"="#4DAF4A", "Non-amniotes"="#984EA3")

prepare_feature_df <- function(path, label) {
  if(!file.exists(path)) {
    message(paste("⚠️ Warning: Data source file not found at", path))
    return(NULL)
  }
  read_tsv(path, show_col_types = FALSE) %>%
    mutate(
      normalized_density = as.numeric(as.character(normalized_density)),
      clade = case_when(
        clade %in% c("Amphibians", "Pisces") ~ "Non-amniotes",
        TRUE ~ clade
      ),
      clade = factor(clade, levels = c("Mammals", "Aves", "Reptiles", "Non-amniotes")),
      feature = factor(label)
    ) %>%
    filter(!is.na(normalized_density))
}

format_plain_p <- function(p) {
  if (is.na(p) || p >= 0.05) return("= ns")
  if (p < 0.001) return("< 0.001")
  return(paste0("= ", sprintf("%.3f", p)))
}

# ==========================================
# 3. DATA HARMONIZATION INGESTION LAYER
# ==========================================
cat("📖 Ingesting all 6 genomic tracking data tables...\n")

df_tss    <- prepare_feature_df(file.path(base_dir, "TSS/All_Vertebrates_G4_Density_GC_Normalized.tsv"), "TSS Flank")
df_5utr   <- prepare_feature_df(file.path(base_dir, "UTRs/5UTR/All_Vertebrates_5UTR_Scaled_G4.tsv"), "5' UTR")
df_exons  <- prepare_feature_df(file.path(base_dir, "CDS/All_Vertebrates_CDS_Scaled_G4.tsv"), "Exons")
df_3utr   <- prepare_feature_df(file.path(base_dir, "UTRs/3UTR/All_Vertebrates_3UTR_Scaled_G4.tsv"), "3' UTR")
df_gene   <- prepare_feature_df(file.path(base_dir, "Gene/All_Vertebrates_Gene_Scaled_G4.tsv"), "Gene body")
df_intron <- prepare_feature_df(file.path(base_dir, "Intron/All_Vertebrates_Intron_Scaled_G4.tsv"), "Intron")

individual_tasks <- list(
  "TSS"      = df_tss,
  "5UTR"     = df_5utr,
  "Exons"    = df_exons,
  "3UTR"     = df_3utr,
  "GeneBody" = df_gene,
  "Intron"   = df_intron
)
individual_tasks <- compact(individual_tasks)

# ==========================================
# 4. MASTER ANALYSIS & GRAPHICS ENGINE
# ==========================================
cat("    Executing Unified Graphics Pipeline Core...\n")

walk2(individual_tasks, names(individual_tasks), function(df_sub, name) {
  display_name <- case_when(
    name == "TSS" ~ "TSS Flank",
    name == "5UTR" ~ "5' UTR",
    name == "3UTR" ~ "3' UTR",
    name == "GeneBody" ~ "Gene body",
    name == "Intron" ~ "Introns",
    TRUE ~ name
  )
  
  # ----------------------------------------
  # A. LOCATION-AWARE STATISTICS GENERATION
  # ----------------------------------------
  if (name == "TSS") {
    stat_summary <- df_sub %>%
      filter(bin_id >= 38 & bin_id <= 43) %>% 
      group_by(species, clade) %>%
      summarise(avg_dens = mean(normalized_density, na.rm = TRUE), .groups = "drop")
  } else {
    stat_summary <- df_sub %>%
      group_by(species, clade) %>%
      summarise(avg_dens = mean(normalized_density, na.rm = TRUE), .groups = "drop")
  }
  
  mammal_vals <- stat_summary %>% filter(clade == "Mammals") %>% pull(avg_dens)
  
  stats_report <- map_df(c("Aves", "Reptiles", "Non-amniotes"), function(trg) {
    trg_vals <- stat_summary %>% filter(clade == trg) %>% pull(avg_dens)
    if(length(mammal_vals) > 0 && length(trg_vals) > 0) {
      wt <- wilcox.test(mammal_vals, trg_vals, exact = FALSE)
      data.frame(
        Feature = display_name,
        Comparison = trg,
        P_val = wt$p.value,
        Mean_Mammals = mean(mammal_vals, na.rm = TRUE),
        Mean_Target = mean(trg_vals, na.rm = TRUE)
      )
    }
  }) %>% mutate(P_adj = p.adjust(P_val, method = "BH"))
  
  write_tsv(stats_report, file.path(output_dir, paste0("Statistics_Output_Summary_", name, ".tsv")))
  
  p_val_aves <- format_plain_p(stats_report$P_adj[stats_report$Comparison == "Aves"])
  p_val_rept <- format_plain_p(stats_report$P_adj[stats_report$Comparison == "Reptiles"])
  p_val_amni <- format_plain_p(stats_report$P_adj[stats_report$Comparison == "Non-amniotes"])
  
  # ----------------------------------------
  # B. LOCAL AGGREGATION & GRAPH AXIS SETUP
  # ----------------------------------------
  df_plot <- df_sub %>% filter(species != "Ficedula_albicollis")
  bin_col <- if("bin_num" %in% colnames(df_plot)) "bin_num" else "bin_id"
  
  profile_agg <- df_plot %>%
    group_by(bin_id = .data[[bin_col]], clade) %>%
    summarise(mean_dens = mean(normalized_density, na.rm = TRUE),
              se_dens = sd(normalized_density, na.rm = TRUE) / sqrt(n()), .groups = "drop")
  
  max_y_ceiling <- max(profile_agg$mean_dens + profile_agg$se_dens, na.rm = TRUE) * 1.05
  
  text_size_val <- 5.8
  
  # CONDITIONALLY HANDLED: Specific x-coordinates and justification configurations per layout type
  if (name == "TSS") {
    grid_breaks   <- seq(1, 80, by = 2)
    center_line   <- 40.5
    scale_breaks  <- seq(1, 80, length.out = 9)
    scale_labels  <- c("-2.0kb", "-1.5", "-1.0", "-0.5", "TSS", "+0.5", "+1.0", "+1.5", "+2.0kb")
    highlight_min <- 38
    highlight_max <- 43
    x_axis_title  <- "Distance from TSS"
    
    text_x_pos    <- 0.1  # Set to left hand corner position
    text_hjust    <- 0    # Left-justified text flow
  } else {
    grid_breaks   <- seq(1, 100, by = 1)
    center_line   <- NA
    scale_breaks  <- c(0, 50, 100)
    scale_labels  <- c("Start", "", "End")
    highlight_min <- NA
    highlight_max <- NA
    x_axis_title  <- paste(display_name, "Bins (scaled)")
    
    text_x_pos    <- 98.0 # Preserved on the right-hand corner
    text_hjust    <- 1    # Right-justified text flow
  }
  
  df_html_annotations <- data.frame(
    x = rep(text_x_pos, 3),
    y = c(max_y_ceiling * 0.95, max_y_ceiling * 0.89, max_y_ceiling * 0.83),
    html_label = c(
      paste0("<span style='color:", clade_colors["Mammals"], "; font-weight:bold;'>Mammals</span> <span style='color:#000000;'>*vs*</span> <span style='color:", clade_colors["Aves"], "; font-weight:bold;'>Aves</span><span style='color:#000000;'>: *P*<sub>adj</sub> ", p_val_aves, "</span>"),
      paste0("<span style='color:", clade_colors["Mammals"], "; font-weight:bold;'>Mammals</span> <span style='color:#000000;'>*vs*</span> <span style='color:", clade_colors["Reptiles"], "; font-weight:bold;'>Reptiles</span><span style='color:#000000;'>: *P*<sub>adj</sub> ", p_val_rept, "</span>"),
      paste0("<span style='color:", clade_colors["Mammals"], "; font-weight:bold;'>Mammals</span> <span style='color:#000000;'>*vs*</span> <span style='color:", clade_colors["Non-amniotes"], "; font-weight:bold;'>Non-amniotes</span><span style='color:#000000;'>: *P*<sub>adj</sub> ", p_val_amni, "</span>")
    )
  )
  
  # ----------------------------------------
  # C. CORE GRAPHICS GENERATION ENGINE
  # ----------------------------------------
  p_base <- ggplot(profile_agg, aes(x = bin_id, y = mean_dens, color = clade, fill = clade)) +
    geom_vline(xintercept = grid_breaks, color = "grey95", linetype = "dashed", linewidth = 0.2)
  
  if (!is.na(center_line)) {
    p_base <- p_base + geom_vline(xintercept = center_line, linetype = "dashed", color = "blue", alpha = 0.5, linewidth = 1.2)
  }
  
  if (!is.na(highlight_min)) {
    p_base <- p_base + annotate("rect", xmin = highlight_min, xmax = highlight_max, ymin = -Inf, ymax = Inf, 
                                fill = "grey50", alpha = 0.10, color = NA)
  }
  
  p_final <- p_base +
    geom_rug(sides = "b", linewidth = 0.8, color = "grey30", alpha = 0.8, length = unit(0.15, "cm")) +
    geom_ribbon(aes(ymin = mean_dens - se_dens, ymax = mean_dens + se_dens), alpha = 0.15, color = NA) +
    geom_line(linewidth = 2.0) +
    
    # Utilizing conditional text_hjust settings to keep alignment clean against edges
    geom_richtext(data = df_html_annotations, 
                  aes(x = x, y = y, label = html_label), 
                  inherit.aes = FALSE, hjust = text_hjust, size = text_size_val,
                  fill = NA, label.color = NA, label.padding = unit(0, "lines")) +
    
    scale_color_manual(values = clade_colors) +
    scale_fill_manual(values = clade_colors) +
    scale_x_continuous(breaks = scale_breaks, labels = scale_labels) +
    theme_classic(base_size = 26) +
    labs(x = x_axis_title, y = "Mean normalized G4 density ± se") +
    theme(
      text = element_text(family = "sans", face = "plain", color = "black"),
      plot.title = element_blank(),
      axis.title.x = element_text(size = 28, margin = margin(t = 12)),
      axis.title.y = element_text(size = 28, margin = margin(r = 12)),
      axis.text = element_text(size = 22, color = "black"),
      axis.line = element_line(linewidth = 1.2, color = "black"),
      axis.ticks = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      legend.position = "none"
    )
  
  file_output_path <- file.path(output_dir, paste0("G4_Embedded_Profile_Master_", name, ".png"))
  ggsave(file_output_path, p_final, width = 11.5, height = 8.5, dpi = 300, bg = "white")
})

# ==========================================
# 5. CONSOLIDATED REPORTING METRICS
# ==========================================
cat("📦 Compiling global consolidated transcript reports (Excluding TSS)...\n")

landscape_elements <- list(df_5utr, df_exons, df_3utr) %>% compact()

if(length(landscape_elements) > 0) {
  box_data <- bind_rows(landscape_elements)
  
  box_summary <- box_data %>%
    group_by(species, clade, feature) %>%
    summarise(avg_dens = mean(normalized_density, na.rm = TRUE), .groups = "drop")
  
  global_landscape_stats <- map_df(levels(box_summary$feature), function(f) {
    feat_data <- box_summary %>% filter(feature == f)
    mammal_vals <- feat_data %>% filter(clade == "Mammals") %>% pull(avg_dens)
    
    map_df(c("Aves", "Reptiles", "Non-amniotes"), function(trg) {
      trg_vals <- feat_data %>% filter(clade == trg) %>% pull(avg_dens)
      if(length(mammal_vals) > 0 && length(trg_vals) > 0) {
        wt <- wilcox.test(mammal_vals, trg_vals, exact = FALSE)
        data.frame(Feature = f, Comparison = paste("Mammals vs", trg),
                   Mean_Mammal = mean(mammal_vals, na.rm=TRUE), Mean_Target = mean(trg_vals, na.rm=TRUE),
                   P_val = wt$p.value)
      }
    })
  }) %>% mutate(P_adj = p.adjust(P_val, method = "BH"))
  
  write_tsv(global_landscape_stats, file.path(output_dir, "Transcript_Landscape_Global_Combined_Statistics.tsv"))
}

cat(paste0("\n✨ SUCCESS! Conditional text layouts compiled natively.\nOutputs exported to: ", output_dir, "\n"))