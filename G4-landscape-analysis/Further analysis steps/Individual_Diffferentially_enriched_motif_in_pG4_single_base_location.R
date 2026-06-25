# 1. LOAD LIBRARIES
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel) 
  library(readxl) 
  library(cowplot) # For isolating and extracting standalone legends cleanly
})

# --- CONFIGURATION (TARGET DIRECTORY FOR FIGURE 3) ---
base_path     <- "Path//Checking"
meta_path     <- file.path(base_path, "Table S1.xlsx") 
clade_out_dir <- "Path//Checking/Final_for_fig_3"

if(!dir.exists(clade_out_dir)) dir.create(clade_out_dir, recursive = TRUE)

datasets <- list(
  Promoter = file.path(base_path, "Motif_Density_Difference_Promoters_G4_minus_Free.xlsx"),
  Genome   = file.path(base_path, "Motif_Density_Difference_Genomewide_G4_minus_Free.xlsx")
)

# Load and Group Metadata from Excel
metadata <- read_excel(meta_path)
meta_clean <- metadata %>%
  mutate(
    Species = gsub(" ", "_", .[[1]]),
    Clade = case_when(
      grepl("Non-", .[[2]], ignore.case = TRUE) ~ "Non-amniotes",
      .[[2]] == "Reptiles" ~ "Reptiles",
      .[[2]] == "Aves"     ~ "Aves",
      .[[2]] == "Mammals"  ~ "Mammals",
      TRUE                 ~ "Other"
    ),
    Thermal = case_when(
      Clade %in% c("Mammals", "Aves") ~ "Homeotherms",
      Clade %in% c("Reptiles", "Non-amniotes") ~ "Poikilotherms",
      TRUE ~ "Other"
    )
  )

# --- ANALYSIS FUNCTION ---
run_comparisons <- function(input_file, ds_name) {
  if(!file.exists(input_file)) {
    cat(sprintf("[WARNING] File not found: %s\n", input_file))
    return(NULL)
  }
  
  df <- read_excel(input_file)
  
  comp_list <- list(
    list("Homeotherms", "Poikilotherms") 
  )
  
  for (comp in comp_list) {
    target_grp <- comp[[1]] 
    ref_grp    <- comp[[2]] 
    
    s_ref <- meta_clean$Species[meta_clean$Clade == ref_grp | meta_clean$Thermal == ref_grp]
    s_tgt <- meta_clean$Species[meta_clean$Clade == target_grp | meta_clean$Thermal == target_grp]
    
    s_ref <- intersect(s_ref, colnames(df))
    s_tgt <- intersect(s_tgt, colnames(df))
    
    if(length(s_ref) < 2 | length(s_tgt) < 2) next
    
    results <- list()
    cat(sprintf("[STATS] Processing %s: %s / %s\n", ds_name, target_grp, ref_grp))
    
    for (i in 1:nrow(df)) {
      vals_ref <- as.numeric(df[i, s_ref]) 
      vals_tgt <- as.numeric(df[i, s_tgt]) 
      
      wt <- wilcox.test(vals_tgt, vals_ref)
      m_ref <- mean(vals_ref, na.rm=TRUE)
      m_tgt <- mean(vals_tgt, na.rm=TRUE)
      
      l2fc <- log2((m_tgt + 0.1) / (m_ref + 0.1))
      
      results[[i]] <- data.frame(
        Motif_ID = df$Motif_ID[i],
        Motif_Name = df$Motif_Name[i],
        Mean_Delta_Ref = m_ref,
        Mean_Delta_Target = m_tgt,
        Log2FC = l2fc,
        P_value = wt$p.value
      )
    }
    
    res_df <- do.call(rbind, results)
    res_df$Padj <- p.adjust(res_df$P_value, method = "BH")
    
    write_tsv(res_df, file.path(clade_out_dir, sprintf("Stats_%s_%s_vs_%s.tsv", ds_name, target_grp, ref_grp)))
    
    # --- VOLCANO PLOT MAP ---
    plot_df <- res_df %>%
      mutate(
        Significance = case_when(
          Padj < 0.05 & Log2FC > 1.0 ~ "Highly Enriched",
          Padj < 0.05 & Log2FC > 0.58 ~ "Moderately Enriched",
          Padj < 0.05 & Log2FC < -1.0 ~ "Highly Reduced",
          Padj < 0.05 & Log2FC < -0.58 ~ "Moderately Reduced",
          TRUE ~ "Not Significant"
        ),
        Search_Name = toupper(Motif_Name),
        Volcano_Score = abs(Log2FC) * (-log10(Padj))
      )
    
    # Limited label acquisition strictly to the top 10 elements
    label_df <- plot_df %>% 
      filter(abs(Log2FC) > 0.58 & Padj < 0.05) %>%
      filter(!Search_Name %in% c("NFIX", "CGGBP1")) %>%
      slice_max(order_by = Volcano_Score, n = 10, with_ties = FALSE)
    
    highlight_targets <- plot_df %>% 
      filter(Search_Name %in% c("NFIX", "CGGBP1"))
    
    max_y_val <- max(-log10(plot_df$Padj), na.rm = TRUE) + 1.0
    
    display_title <- if_else(ds_name == "Genome", 
                             "Genome-wide", 
                             "Promoter")
    
    # Build core graphics plot container
    p_base <- ggplot(plot_df, aes(x = Log2FC, y = -log10(Padj))) +
      geom_point(aes(color = Significance), alpha = 0.6, size = 8) +
      geom_point(
        data = highlight_targets,
        color = "black",       
        fill = NA,             
        shape = 21,            
        size = 12, 
        stroke = 1.5           
      ) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey70", linewidth = 0.4) +
      geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey40", linewidth = 0.4) +
      geom_vline(xintercept = c(-1, -0.58, 0.58, 1), linetype = "dashed", color = "grey70", linewidth = 0.3) +
      geom_text_repel(
        data = label_df, 
        aes(label = Motif_Name),
        size = 14,
        color = "black",
        max.overlaps = 100,
        fontface = "italic",
        segment.size = 0.3,
        segment.color = "grey60"
      ) +
      scale_color_manual(values = c(
        "Highly Enriched"     = "#ca0020", 
        "Moderately Enriched" = "#f4a582", 
        "Highly Reduced"      = "#0571b0", 
        "Moderately Reduced"  = "#92c5de", 
        "Not Significant"     = "grey90"
      )) +
      # MODIFICATION: Maintained physical limits at -7 to 7 but labeled ticks by 2 for a clean 1-point gap sequence
      scale_x_continuous(limits = c(-7, 7), breaks = seq(-6, 6, by = 2)) + 
      scale_y_continuous(limits = c(0, max_y_val)) +
      theme_classic(base_size = 15) + 
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey98", linewidth = 0.4),
        plot.title    = element_text(face = "plain", size = 46, hjust = 0.5), 
        plot.subtitle = element_text(face = "plain", size = 44, hjust = 0.5), 
        axis.title    = element_text(face = "plain", size = 40),
        axis.text     = element_text(color = "black", size = 40, face = "plain"),
        legend.title  = element_text(face = "plain", size = 24),
        legend.text   = element_text(face = "plain", size = 24),
        legend.position  = "bottom",          
        legend.direction = "horizontal",
        legend.box       = "horizontal"
      ) +
      labs(
        title = display_title,
        subtitle = bquote(plain("Δ Mean motif density (pG4 - pG4-free background)")),
        x = bquote(log[2] ~ "(" * plain("Δ") * ~ Homeotherms / plain("Δ") * ~ Poikilotherms * ")"),
        y = expression(paste("-log"[10], " (Adjusted P-value)"))
      )
    
    # =========================================================================
    # EXTRACT AND SAVE INDEPENDENT PLOT COMPONENTS
    # =========================================================================
    
    # 1. Capture the vertical legend object module explicitly 
    standalone_legend <- cowplot::get_legend(p_base)
    
    # 2. Strip the legend out of the main plot container completely
    p_clean_main <- p_base + theme(legend.position = "none")
    
    # 3. Export the clean Main Canvas (No Legend)
    ggsave(
      filename = file.path(clade_out_dir, sprintf("Volcano_%s_%s_vs_%s_NoLegend.png", ds_name, target_grp, ref_grp)), 
      plot     = p_clean_main, 
      width    = 18, 
      height   = 16, 
      dpi      = 300
    )
    
    # 4. Export the standalone Legend Asset Image
    ggsave(
      filename = file.path(clade_out_dir, sprintf("Legend_Asset_Vertical_%s_%s_vs_%s.png", ds_name, target_grp, ref_grp)), 
      plot     = standalone_legend, 
      width    = 20, 
      height   = 2, 
      dpi      = 300
    )
    
    cat(sprintf("[SUCCESS] Saved separated graphics for %s:\n  -> Volcano Main Canvas (X-axis ticks at 1-point intervals: -6, -4, -2, 0, 2, 4, 6)\n  -> Standalone Vertical Legend Asset\n", ds_name))
  }
}

# --- EXECUTE INDEPENDENT CONTEXTS ---
run_comparisons(datasets$Promoter, "Promoter")
run_comparisons(datasets$Genome, "Genome")

cat("\n[ALL DONE] Clear layout spacing script operation finished successfully!\n")