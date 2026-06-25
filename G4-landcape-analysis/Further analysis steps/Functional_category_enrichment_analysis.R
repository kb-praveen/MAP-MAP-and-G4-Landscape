# 1. LOAD LIBRARIES
suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr) 
  library(ggvenn) # Ensure this is installed: install.packages("ggvenn")
})

# --- SETUP ---
master_file <- "/Users/praveenbishnoi/Desktop/G4_Final_Analysis_Complete/GC_bias_analysis_Final/Master_Retention_Analysis_Full.tsv"
new_dir     <- "/Users/praveenbishnoi/Desktop/G4_Final_Analysis_Complete/GC_bias_analysis_Final/New"
dir.create(new_dir, showWarnings = FALSE, recursive = TRUE)

# --- 2. LOAD & CATEGORIZE MOTIFS ---
df <- read_tsv(master_file, show_col_types = FALSE)

# Strict 4-column identity selection for Top 50
top_abs_df <- df %>% 
  arrange(desc(Absolute_G4_Shift)) %>% 
  head(50) %>% 
  mutate(Combined_ID = paste(Motif_ID, Motif_Name, Enrichment_Tier, Region, sep="_"))

top_norm_df <- df %>% 
  arrange(desc(Delta_Delta_Shift)) %>% 
  head(50) %>% 
  mutate(Combined_ID = paste(Motif_ID, Motif_Name, Enrichment_Tier, Region, sep="_"))

# --- 3. VENN DIAGRAM INTEGRATION ---
# Defining the list for the Venn
venn_list <- list(
  "Absolute Retention" = top_abs_df$Combined_ID,
  "Normalized Retention" = top_norm_df$Combined_ID
)

p_venn <- ggvenn(
  venn_list, 
  fill_color = c("#ca0020", "#0571b0"),
  stroke_size = 1.2, 
  set_name_size = 12, # Large for visibility
  text_size = 10      # Large for visibility
) +
  labs(title = "Top 50 Motif Overlap",
       subtitle = "Based on ID, Name, Tier, and Region") +
  theme(plot.title = element_text(size = 42, hjust = 0.5),
        plot.subtitle = element_text(size = 36, hjust = 0.5, color = "grey30"))

ggsave(file.path(new_dir, "Motif_Overlap_Venn_Final.png"), p_venn, width = 15, height = 15, dpi = 300)

# --- 4. PREP FOR GO ANALYSIS ---
motif_groups <- list(
  "Common" = intersect(top_abs_df$Combined_ID, top_norm_df$Combined_ID),
  "Absolute[unique]" = setdiff(top_abs_df$Combined_ID, top_norm_df$Combined_ID),
  "Normalized[unique]" = setdiff(top_norm_df$Combined_ID, top_abs_df$Combined_ID)
)

all_results <- list()
for (group_name in names(motif_groups)) {
  ids <- motif_groups[[group_name]]
  names_only <- sapply(strsplit(ids, "_"), `[`, 2) 
  genes <- unique(names_only)
  
  gene_conv <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  res <- enrichGO(gene = gene_conv$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF", 
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  
  if(!is.null(res) && nrow(as.data.frame(res)) > 0) {
    res_df <- as.data.frame(res) %>%
      mutate(Group = group_name,
             GeneRatio_Val = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)),
             BgRatio_Val   = as.numeric(sub("/.*", "", BgRatio)) / as.numeric(sub(".*/", "", BgRatio)),
             Fold_Enrichment = round(GeneRatio_Val / BgRatio_Val, 2))
    all_results[[group_name]] <- res_df
  }
}

master_stats <- bind_rows(all_results)
write_tsv(master_stats, file.path(new_dir, "Integrated_MF_Enrichment_Stats_New.tsv"))

# --- 5. DOT PLOT (ULTRA SCALE & WRAPPING) ---
plot_data <- master_stats %>%
  group_by(Group) %>%
  slice_min(order_by = p.adjust, n = 7) %>%
  ungroup() %>%
  mutate(Description = str_wrap(Description, width = 40))

x_labs <- c(
  "Common" = "Common",
  "Absolute[unique]" = expression(Absolute[unique]),
  "Normalized[unique]" = expression(Normalized[unique])
)

p_dot <- ggplot(plot_data, aes(x = Group, y = Description)) +
  geom_point(aes(size = Fold_Enrichment, color = p.adjust)) +
  scale_color_gradient(low = "#ca0020", high = "#0571b0", name = expression(italic(P)[adj])) +
  scale_size_continuous(
    name = "Fold Enrichment",
    range = c(12, 32), 
    breaks = c(50, 100, 150, 200),
    limits = c(0, 220)
  ) +
  scale_x_discrete(labels = x_labs) +
  theme_bw() + 
  labs(
    title = "Molecular Function Enrichment",
    subtitle = "Top GC-retaining TFBS in Homeotherms",
    x = "Motifs", 
    y = "GO: Molecular Function"
  ) +
  theme(
    plot.title = element_text(size = 42),
    plot.subtitle = element_text(size = 40, color = "black"),
    axis.text.y = element_text(size = 40, color = "black", lineheight = 0.8),
    axis.text.x = element_text(size = 40, color = "black"),
    axis.title = element_text(size = 40),
    legend.title = element_text(size = 40),
    legend.text = element_text(size = 35),
    legend.key.height = unit(3.5, "cm"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.7),
    plot.margin = margin(2, 2, 2, 2, "cm")
  )

ggsave(file.path(new_dir, "MF_Enrichment_Final_Publication.png"), p_dot, width = 30, height = 22, dpi = 300)

cat("\n[SUCCESS] Venn diagram and Dots saved in: ", new_dir)