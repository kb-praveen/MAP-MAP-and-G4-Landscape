suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

# --- CONFIGURATION ---
base_path <- "/scratch/kumarpraveen/G4_motif_analysis/Results"
meta_path <- "/scratch/kumarpraveen/G4_motif_analysis/Table_S1.tsv"

clade_out_dir <- file.path(base_path, "Clade_Comparisons")
if(!dir.exists(clade_out_dir)) dir.create(clade_out_dir, recursive = TRUE)

datasets <- list(
  Promoter = file.path(base_path, "Motif_Density_Difference_Promoters_G4_minus_Free.tsv"),
  Genome   = file.path(base_path, "Motif_Density_Difference_Genomewide_G4_minus_Free.tsv")
)

# 1. Load and Group Metadata
metadata <- read_tsv(meta_path, show_col_types = FALSE)
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
  if(!file.exists(input_file)) return(NULL)
  df <- read_tsv(input_file, show_col_types = FALSE)
  
  comp_list <- list(
    list("Aves", "Mammals"),
    list("Reptiles", "Mammals"),
    list("Non-amniotes", "Mammals"),
    list("Poikilotherms", "Homeotherms")
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
    
    # --- VOLCANO PLOT ---
    plot_df <- res_df %>%
      mutate(Significance = case_when(
        Padj < 0.05 & Log2FC > 1.0 ~ "Strongly Enriched",
        Padj < 0.05 & Log2FC > 0.58 ~ "Moderately Enriched",
        Padj < 0.05 & Log2FC < -1.0 ~ "Strongly Reduced",
        Padj < 0.05 & Log2FC < -0.58 ~ "Moderately Reduced",
        TRUE ~ "Not Significant"
      ))

    p <- ggplot(plot_df, aes(x = Log2FC, y = -log10(Padj), color = Significance)) +
      geom_point(alpha = 0.6, size = 1) +
      # Horizontal Significance Lines
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey70", size = 0.4) +
      geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey40", size = 0.4) +
      # Vertical Fold-Change Lines
      geom_vline(xintercept = c(-1, -0.58, 0.58, 1), linetype = "dashed", color = "grey70", size = 0.3) +
      scale_color_manual(values = c(
        "Strongly Enriched" = "#ca0020", 
        "Moderately Enriched" = "#f4a582", 
        "Strongly Reduced" = "#0571b0", 
        "Moderately Reduced" = "#92c5de", 
        "Not Significant" = "grey90"
      )) +
      # --- AXIS TICK LOGIC ---
      # Major ticks every 2 (labelled), Minor every 0.5 (unlabelled)
      scale_x_continuous(
        breaks = seq(-20, 20, by = 2), 
        minor_breaks = seq(-20, 20, by = 0.5)
      ) +
      scale_y_continuous(
        breaks = seq(0, 100, by = 2), 
        minor_breaks = seq(0, 100, by = 0.5)
      ) +
      theme_classic() +
      # Force minor grid lines to be visible if the theme hides them
      theme(
        panel.grid.minor = element_line(color = "grey95", size = 0.2),
        legend.position = "right", 
        plot.title = element_text(face="bold"),
        legend.title = element_text(face="bold")
      ) +
      labs(
        title = sprintf("%s: %s vs %s", ds_name, target_grp, ref_grp),
        subtitle = expression(paste("Wilcoxon Rank Sum Test on ", Delta, " Mean (G4 - G4-free background)")),
        x = bquote(log[2] ~ "(" * Delta * ~ .(target_grp) / Delta * ~ .(ref_grp) * ")"),
        y = expression(paste("-log"[10], " (Adjusted P-value)"))
      )

    ggsave(file.path(clade_out_dir, sprintf("Volcano_%s_%s_vs_%s.png", ds_name, target_grp, ref_grp)), 
           p, width = 8, height = 6, dpi = 300)
  }
}

# --- EXECUTE ---
run_comparisons(datasets$Promoter, "Promoter")
run_comparisons(datasets$Genome, "Genome")

cat("\n[SUCCESS] Completed 8 comparisons with refined axis ticks.\n")