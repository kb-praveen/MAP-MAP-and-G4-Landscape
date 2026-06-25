# 1. LOAD LIBRARIES
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel) 
  library(readxl) 
  library(cowplot) # Required to extract and save the legend separately
})

# --- CONFIGURATION (TARGET DIRECTORY FOR FIGURE 3) ---
base_path     <- "Path/Checking"
meta_path     <- file.path(base_path, "Table S1.xlsx") 
clade_out_dir <- "Path/Checking/Final_for_fig_3"

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

# Global holder variable to safely catch the legend object
extracted_shared_legend <- NULL

# --- STANDALONE PLOT GENERATION PIPELINE ---
run_comparisons <- function(input_file, ds_name) {
  if(!file.exists(input_file)) {
    cat(sprintf("[WARNING] File not found: %s\n", input_file))
    return(NULL)
  }
  
  df <- read_excel(input_file)
  
  target_grp <- "Homeotherms"
  ref_grp    <- "Poikilotherms"
  
  s_ref <- meta_clean$Species[meta_clean$Clade == ref_grp | meta_clean$Thermal == ref_grp]
  s_tgt <- meta_clean$Species[meta_clean$Clade == target_grp | meta_clean$Thermal == target_grp]
  
  s_ref <- intersect(s_ref, colnames(df))
  s_tgt <- intersect(s_tgt, colnames(df))
  
  if(length(s_ref) < 2 | length(s_tgt) < 2) return(NULL)
  
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
  
  # Map points to clean discrete categories
  plot_df <- res_df %>%
    mutate(
      Significance = case_when(
        Padj < 0.05 & Log2FC > 1.0   ~ "Highly_Enriched",
        Padj < 0.05 & Log2FC > 0.58  ~ "Moderately_Enriched",
        Padj < 0.05 & Log2FC < -1.0  ~ "Highly_Reduced",
        Padj < 0.05 & Log2FC < -0.58 ~ "Moderately_Reduced",
        TRUE                         ~ "Not_Significant"
      ),
      Search_Name = toupper(Motif_Name),
      Volcano_Score = abs(Log2FC) * (-log10(Padj))
    )
  
  label_df <- plot_df %>% 
    filter(abs(Log2FC) > 0.58 & Padj < 0.05) %>%
    filter(!Search_Name %in% c("NFIX", "CGGBP1")) %>%
    slice_max(order_by = Volcano_Score, n = 10, with_ties = FALSE)
  
  highlight_targets <- plot_df %>% 
    filter(Search_Name %in% c("NFIX", "CGGBP1"))
  
  max_y_val <- max(-log10(plot_df$Padj), na.rm = TRUE) + 1.0
  display_subtitle <- if_else(ds_name == "Genome", "Genome-wide", "Promoters")
  
  p_base <- ggplot(plot_df, aes(x = Log2FC, y = -log10(Padj))) +
    geom_point(aes(color = Significance), alpha = 0.6, size = 9.5) +
    geom_point(
      data = highlight_targets,
      color = "black", fill = NA, shape = 21, size = 9.5, stroke = 1.5           
    ) +
    
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_vline(xintercept = c(-1, -0.58, 0.58, 1), linetype = "dashed", color = "black", linewidth = 0.3) +
    
    geom_text_repel(
      data = label_df, aes(label = toupper(Motif_Name)),
      size = 13, 
      color = "black", 
      max.overlaps = 100, 
      fontface = "italic",
      force = 3,                   
      box.padding = 0.6,           
      point.padding = 0.5,         
      segment.size = 0.5,          
      segment.color = "black",     
      min.segment.length = 0.1,    
      direction = "both"           
    ) +
    
    scale_color_manual(
      values = c(
        "Highly_Enriched"     = "#ca0020", 
        "Moderately_Enriched" = "#f4a582", 
        "Highly_Reduced"      = "#0571b0", 
        "Moderately_Reduced"  = "#92c5de", 
        "Not_Significant"     = "grey90"
      ),
      labels = c(
        "Highly_Enriched"     = "Highly Enriched\n(Log₂FC ≥ 1.0)",
        "Moderately_Enriched" = "Moderately Enriched\n(0.58 ≤ Log₂FC < 1.0)",
        "Highly_Reduced"      = "Highly Reduced\n(Log₂FC ≤ -1.0)",
        "Moderately_Reduced"  = "Moderately Reduced\n(-1.0 < Log₂FC ≤ -0.58)",
        "Not_Significant"     = "Not Significant\n(Baseline / ns)"
      )
    ) +
    
    guides(
      color = guide_legend(
        title.vjust = 1.0,
        label.vjust = 1.0,  
        label.theme = element_text(size = 32, face = "plain", lineheight = 0.9, vjust = 1.0)
      )
    ) +
    
    scale_x_continuous(limits = c(-7, 7), breaks = seq(-6, 6, by = 2)) + 
    scale_y_continuous(limits = c(0, max_y_val)) +
    theme_classic(base_size = 15) + 
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey98", linewidth = 0.4),
      
      plot.title    = element_text(face = "plain", size = 42, hjust = 0.5, margin = margin(b = 25)), 
      plot.subtitle = element_text(face = "plain", size = 42, hjust = 0.5, margin = margin(b = 25)), 
      axis.title.x  = element_text(face = "plain", size = 42, margin = margin(t = 24)),
      axis.title.y  = element_text(face = "plain", size = 42, margin = margin(r = 24)),
      axis.text.x   = element_text(color = "black", size = 42, face = "plain"),
      axis.text.y   = element_text(color = "black", size = 42, face = "plain"),
      
      legend.title  = element_text(face = "plain", size = 32, vjust = 1.0, lineheight = 1.1), 
      legend.position  = "bottom",          
      legend.direction = "horizontal",
      legend.box       = "horizontal",
      
      # Clean horizontal separation gaps
      legend.spacing.x     = unit(0.9, "cm"),  
      legend.key.spacing.x = unit(4.5, "cm"),  
      
      plot.margin      = margin(1, 1, 1, 1, "cm")
    ) +
    labs(
      title = bquote(plain("Δ Mean motif density (pG4 - pG4-free background)")),
      subtitle = display_subtitle,
      x = bquote(log[2] ~ "(" * plain("Δ") * ~ Homeotherms / plain("Δ") * ~ Poikilotherms * ")"),
      y = expression(paste("-log"[10], " (Adjusted P-value)")),
      color = "Enrichment Categories\n(All Tiers: Padj < 0.05)" 
    )
  
  # Harvest the legend asset on the first pass
  if (is.null(extracted_shared_legend)) {
    extracted_shared_legend <<- cowplot::get_legend(p_base)
  }
  
  # Strip internal legends completely before exporting individual volcano structures
  p_clean_main <- p_base + theme(legend.position = "none")
  
  ggsave(
    filename = file.path(clade_out_dir, sprintf("Volcano_%s_%s_vs_%s_NoLegend.png", ds_name, target_grp, ref_grp)), 
    plot     = p_clean_main, 
    width    = 18, 
    height   = 18.5, 
    dpi      = 300
  )
}

# --- 2. EXECUTE VISUAL GENERATION PASS ---
cat("[INFO] Running graphics pipeline pass...\n")
run_comparisons(datasets$Genome, "Genome")
run_comparisons(datasets$Promoter, "Promoter")

# --- 3. EXPORT EXCLUSIVELY SEPARATE SHARED LEGEND ASSET ---
if (!is.null(extracted_shared_legend)) {
  cat("[INFO] Exporting extra-wide global shared legend asset (Bypassing 50-inch safety boundary)...\n")
  ggsave(
    filename   = file.path(clade_out_dir, "Figure3_Shared_Legend_Asset.png"),
    plot       = extracted_shared_legend,
    width      = 54, 
    height     = 3.5, 
    dpi        = 300,
    limitsize  = FALSE # MODIFICATION: Fixed. Disables the 50-inch internal boundary gatekeeper
  )
}

cat("\n📈 SUCCESS! Over-sized horizontal vector asset written out safely to disk.\n")