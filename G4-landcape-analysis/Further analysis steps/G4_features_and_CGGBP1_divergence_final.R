# ==========================================
# 1. LOAD LIBRARIES
# ==========================================
suppressPackageStartupMessages({
  library(ape)
  library(geiger)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(phytools)
  library(openxlsx)
  library(ggrepel)
  library(cowplot)
})

# ==========================================
# 2. SETUP & PATHS
# ==========================================
tree_path    <- "/Users/praveenbishnoi/Desktop/G4/G4_analysis/Species.nwk"
protein_path <- "/Users/praveenbishnoi/Downloads/CGGBP1_alignment_aa_processed_rates.tsv"

fig_dir      <- "/Users/praveenbishnoi/Desktop/G4/G4_revision/Absolute_GC_analysis"
tab_dir      <- "/Users/praveenbishnoi/Desktop/G4/G4_revision/Absolute_GC_analysis"

if(!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if(!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE)

clade_order  <- c("Global", "Mammals", "Aves", "Reptiles")

colors_raw <- c(
  "Global"   = "#000000", 
  "Mammals"  = "#E41A1C", 
  "Aves"     = "#377EB8", 
  "Reptiles" = "#4DAF4A"  
)

# FIXED: Global color profile updated to black to maintain consistency across spaces
colors_pic <- c(
  "Global"   = "#000000", 
  "Mammals"  = "#E41A1C", 
  "Aves"     = "#377EB8", 
  "Reptiles" = "#4DAF4A"  
)

trait_tasks <- list(
  list(
    short = "Absolute_GC_Mean",
    g_col = "Absolute_G_plus_C_count_mean", p_col = "Absolute_G_plus_C_count_mean", 
    label = "Absolute G+C Count Mean"
  ),
  list(
    short = "GC_Content",
    g_col = "Genome_GC_content", p_col = "Promoters_GC_content", 
    label = "GC Content (%)"
  ),
  list(
    short = "G4_Density",
    g_col = "G4_per_Mb", p_col = "G4_per_Mb", 
    label = "G4 Density (per Mb)"
  ),
  list(
    short = "G4_Length",
    g_col = "Average_G4_Length", p_col = "Average_G4_Length", 
    label = "Average G4 Length (bp)"
  )
)

sub_clades <- c("Mammals", "Aves", "Reptiles")
target_y   <- "CGGBP1_divergence"

# ==========================================
# 3. DATA HARMONIZATION LAYER
# ==========================================
tree_master <- read.tree(tree_path)
df_protein  <- read.delim(protein_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

df_genomes_raw <- read.delim("/Users/praveenbishnoi/Desktop/G4/G4_analysis/Metadata/Genomes.tsv", header=TRUE, sep="\t", check.names=FALSE)
df_promots_raw <- read.delim("/Users/praveenbishnoi/Desktop/G4/G4_analysis/Metadata/Promoters.tsv", header=TRUE, sep="\t", check.names=FALSE)

clean_and_sync_dataset <- function(raw_df, tree_obj, protein_df) {
  df_mod <- raw_df %>%
    mutate(
      Species = gsub("Neovison_vison", "Neogale_vison", Species),
      Species = gsub("Stachyris_ruficeps", "Cyanoderma_ruficeps", Species),
      TreeNames = gsub(" ", "_", Species),
      Classification_1 = gsub("Non- Amniotes", "Non-amniotes", Classification_1)
    )
  rownames(df_mod) <- df_mod$TreeNames
  
  matched <- treedata(tree_obj, df_mod)
  df_res  <- as.data.frame(matched$data, stringsAsFactors = FALSE)
  df_res$Classification_1 <- df_mod[rownames(df_res), "Classification_1"]
  df_res$TreeNames <- rownames(df_res)
  
  df_res <- left_join(df_res, protein_df, by = "TreeNames")
  rownames(df_res) <- df_res$TreeNames
  
  df_clean <- df_res %>%
    filter(!is.na(CGGBP1_divergence)) %>%
    filter(Classification_1 != "Non-amniotes")
  return(df_clean)
}

df_genomes_master <- clean_and_sync_dataset(df_genomes_raw, tree_master, df_protein)
df_promots_master <- clean_and_sync_dataset(df_promots_raw, tree_master, df_protein)

tree_genomes_clean <- keep.tip(tree_master, rownames(df_genomes_master))
tree_promots_clean <- keep.tip(tree_master, rownames(df_promots_master))

get_p_string <- function(p) {
  if (is.nan(p) || is.na(p) || p > 0.05) return("= ns")
  if (p < 0.001) return("< 0.001")
  return(paste0("= ", round(p, 3)))
}

# ==========================================
# 4. INDIVIDUAL PLOT GENERATION LOOP
# ==========================================
wb_master_stats <- createWorkbook()

for (task in trait_tasks) {
  summary_stats <- data.frame()
  
  env_config <- list(
    list(mode = "Genomewide", df = df_genomes_master, tree = tree_genomes_clean, col = task$g_col),
    list(mode = "Promoters",  df = df_promots_master, tree = tree_promots_clean, col = task$p_col)
  )
  
  for (env in env_config) {
    df_working   <- env$df
    tree_working <- env$tree
    curr_x       <- env$col
    
    df_working[[curr_x]]   <- as.numeric(as.character(df_working[[curr_x]]))
    df_working[[target_y]] <- as.numeric(as.character(df_working[[target_y]]))
    
    # Statistical calculations
    raw_all  <- cor.test(df_working[[curr_x]], df_working[[target_y]], method = "spearman", exact = FALSE)
    pic_x_all  <- pic(setNames(df_working[[curr_x]], rownames(df_working)), tree_working)
    pic_y_all  <- pic(setNames(df_working[[target_y]], rownames(df_working)), tree_working)
    pic_lm_all <- lm(pic_y_all ~ pic_x_all + 0)
    
    summary_stats <- rbind(summary_stats, data.frame(
      Dataset = env$mode, Clade = "Global", 
      Raw_r = round(raw_all$estimate, 3), Raw_P = raw_all$p.value,
      PIC_AdjR2 = round(summary(pic_lm_all)$adj.r.squared, 3), PIC_P = summary(pic_lm_all)$coefficients[1, 4]
    ))
    
    clade_pic_trends <- data.frame()
    
    for (grp in sub_clades) {
      df_sub <- df_working[df_working$Classification_1 == grp, ]
      
      if (nrow(df_sub) > 3) {
        tree_sub <- keep.tip(tree_working, rownames(df_sub))
        raw_sub  <- cor.test(df_sub[[curr_x]], df_sub[[target_y]], method = "spearman", exact = FALSE)
        
        px_sub  <- pic(setNames(df_sub[[curr_x]], rownames(df_sub)), tree_sub)
        py_sub  <- pic(setNames(df_sub[[target_y]], rownames(df_sub)), tree_sub)
        pic_lm_sub <- lm(py_sub ~ px_sub + 0)
        
        clade_pic_trends <- rbind(clade_pic_trends, data.frame(px = px_sub, py = py_sub, Clade = grp))
        
        summary_stats <- rbind(summary_stats, data.frame(
          Dataset = env$mode, Clade = grp, 
          Raw_r = round(raw_sub$estimate, 3), Raw_P = raw_sub$p.value,
          PIC_AdjR2 = round(summary(pic_lm_sub)$adj.r.squared, 3), PIC_P = summary(pic_lm_sub)$coefficients[1, 4]
        ))
      }
    }
    
    # Hide visual outlier point safely from plot frames
    df_plot <- df_working %>% filter(TreeNames != "Ficedula_albicollis")
    df_focal_raw <- df_plot %>%
      filter(TreeNames %in% c("Anolis_carolinensis", "Gallus_gallus", "Homo_sapiens")) %>%
      mutate(Label = case_when(
        TreeNames == "Anolis_carolinensis" ~ "Ac",
        TreeNames == "Gallus_gallus"       ~ "Gg",
        TreeNames == "Homo_sapiens"         ~ "Hs"
      ))
    
    df_pic_all <- data.frame(px = pic_x_all, py = pic_y_all, NodeNames = names(pic_x_all), stringsAsFactors = FALSE)
    df_pic_plot <- df_pic_all %>% filter(NodeNames != "Ficedula_albicollis")
    df_focal_pic <- df_pic_plot %>%
      filter(NodeNames %in% c("Anolis_carolinensis", "Gallus_gallus", "Homo_sapiens")) %>%
      mutate(Label = case_when(
        NodeNames == "Anolis_carolinensis" ~ "Ac",
        NodeNames == "Gallus_gallus"       ~ "Gg",
        NodeNames == "Homo_sapiens"         ~ "Hs"
      ))
    
    theme_nature <- theme_classic(base_size = 26) + 
      theme(
        legend.position = "none", 
        plot.title = element_blank(), 
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
        axis.title = element_text(size = 28, face = "plain", colour = "black"), 
        axis.text = element_text(size = 22, colour = "black"),
        axis.line = element_line(linewidth = 1.0, colour = "black"),
        axis.ticks = element_line(linewidth = 1.0, colour = "black")
      )
    
    # Panel A: Raw Scatterplot (Global = Dashed, Clades = Solid)
    p_corr <- ggplot(df_plot, aes(x = .data[[curr_x]], y = .data[[target_y]])) +
      geom_point(color = "#E5E5E5", size = 4.5, alpha = 0.40, shape = 16) +
      geom_point(data = filter(df_plot, Classification_1 %in% sub_clades), 
                 aes(color = Classification_1), size = 4.5, alpha = 0.70, shape = 16) +
      geom_smooth(data = df_working, aes(x = .data[[curr_x]], y = .data[[target_y]]),
                  method = "lm", color = colors_raw["Global"], linetype = "dashed", linewidth = 2.5, se = FALSE) +
      geom_smooth(data = filter(df_working, Classification_1 %in% sub_clades), 
                  aes(color = Classification_1, group = Classification_1), 
                  method = "lm", linetype = "solid", linewidth = 2.8, se = FALSE) +
      geom_point(data = df_focal_raw, aes(fill = Classification_1), color = "#000000", size = 5.2, shape = 21, stroke = 1.0) +
      geom_text_repel(data = df_focal_raw, aes(label = Label), size = 6.5, color = "black", 
                      box.padding = 1.8, point.padding = 0.8, max.overlaps = Inf, fontface = "plain",
                      min.segment.length = 0, segment.color = "black", segment.linewidth = 0.5) +
      scale_fill_manual(values = colors_raw) + 
      scale_color_manual(values = colors_raw) + 
      theme_nature +
      labs(x = task$label, y = "Divergence from Hs CGGBP1")
    
    # Panel B: PIC Space Plot (Global = Solid, Clades = Dashed)
    p_pic <- ggplot(df_pic_plot, aes(x = px, y = py)) +
      geom_point(alpha = 0.30, color = "#B3B3B3", size = 4.0, shape = 16) +
      geom_smooth(data = data.frame(px = pic_x_all, py = pic_y_all), aes(x = px, y = py),
                  method = "lm", formula = y ~ x + 0, color = colors_pic["Global"], linetype = "solid", linewidth = 2.5, se = FALSE) +
      geom_smooth(data = clade_pic_trends, aes(x = px, y = py, color = Clade, group = Clade), 
                  method = "lm", formula = y ~ x + 0, linetype = "dashed", linewidth = 2.8, se = FALSE) +
      geom_point(data = df_focal_pic, aes(x = px, y = py), color = "#000000", fill = "#FFFFFF", size = 4.8, shape = 21, stroke = 1.0) +
      geom_text_repel(data = df_focal_pic, aes(x = px, y = py, label = Label), size = 6.5, color = "black",
                      box.padding = 1.8, point.padding = 0.8, max.overlaps = Inf, fontface = "plain",
                      min.segment.length = 0, segment.color = "black", segment.linewidth = 0.5) +
      scale_color_manual(values = colors_pic) +
      theme_nature + 
      labs(x = paste0("Contrast: ", task$label), y = "Contrast: Divergence from Hs CGGBP1")
    
    composite_side_by_side <- p_corr + p_pic + plot_layout(ncol = 2)
    
    file_output_grid <- file.path(fig_dir, paste0("Figure_Combined_Divergence_SideBySide_", env$mode, "_", task$short, ".png"))
    ggsave(file_output_grid, composite_side_by_side, width = 17.0, height = 8.0, dpi = 300, bg = "white")
  }
  addWorksheet(wb_master_stats, task$short)
  writeData(wb_master_stats, task$short, summary_stats)
}
saveWorkbook(wb_master_stats, file.path(tab_dir, "Summary_Calculated_Amniota_Stats.xlsx"), overwrite = TRUE)
cat("\n*** Divergence Master Pipeline Complete! Lines homogenized to black. ***\n")