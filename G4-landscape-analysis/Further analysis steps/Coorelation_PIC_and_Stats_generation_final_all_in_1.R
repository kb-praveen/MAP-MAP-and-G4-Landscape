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
  library(grid)
  library(gridExtra)
  library(openxlsx)
  library(cowplot)
})

# ==========================================
# 2. SETUP & PATHS
# ==========================================
tree_path <- "Path/G4_analysis/Species.nwk"
fig_dir   <- "Path/G4_revision/Figure_1_revision"
tab_dir   <- "Path/G4_revision/Figure_1_revision"

if(!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if(!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE)

clade_order <- c("Mammals", "Aves", "Reptiles", "Non-amniotes", "Others")
custom_colors <- c(
  "Mammals"      = "#E41A1C", 
  "Aves"         = "#377EB8", 
  "Reptiles"     = "#4DAF4A", 
  "Non-amniotes" = "#984EA3",
  "Others"       = "#D3D3D3"
)

datasets_config <- list(
  `Genomewide` = list(path = "Path/G4_analysis/Metadata/Genomes.tsv", 
                      gc_col = "Genome_GC_content"),
  `Promoters`  = list(path = "Path/G4_analysis/Metadata/Promoters.tsv", 
                      gc_col = "Promoters_GC_content")
)

analysis_tasks <- list(
  list(y_col = "G4_per_Mb", x_label = "GC-content (%)", y_label = "G4 density (per Mb)", short = "Density"),
  list(y_col = "G4_GC_Content_Percent", x_label = "GC-content (%)", y_label = "G4 GC content (%)", short = "Composition"),
  list(x_col_override = "Average_G4_Length", y_col = "G4_GC_Content_Percent", x_label = "Average G4 length (bp)", y_label = "G4 GC content (%)", short = "Length_GC")
)

sub_clades <- c("Mammals", "Aves", "Reptiles", "Non-amniotes")

# ==========================================
# 3. PROCESSING & GENERATION LOOP
# ==========================================
tree <- read.tree(tree_path)

for (data_name in names(datasets_config)) {
  conf <- datasets_config[[data_name]]
  all_panel_plots <- list()
  deviation_summary <- data.frame() 
  
  metadata <- read.delim(conf$path, header = TRUE, sep="\t", check.names = FALSE) %>%
    mutate(
      Species = gsub("Neovison_vison", "Neogale_vison", Species),
      Species = gsub("Stachyris_ruficeps", "Cyanoderma_ruficeps", Species),
      TreeNames = gsub(" ", "_", Species),
      Classification_1 = gsub("Non- Amniotes", "Non-amniotes", Classification_1)
    )
  rownames(metadata) <- metadata$TreeNames
  
  matched <- treedata(tree, metadata)
  tree_all <- matched$phy
  df_all <- as.data.frame(matched$data, stringsAsFactors = FALSE)
  
  numeric_cols <- c(conf$gc_col, "G4_per_Mb", "G4_GC_Content_Percent", "Average_G4_Length")
  df_all[numeric_cols] <- lapply(df_all[numeric_cols], function(x) as.numeric(as.character(x)))
  
  # Hide outlier bird visually from plots
  df_plot <- df_all %>% filter(TreeNames != "Ficedula_albicollis")
  
  wb_stats <- createWorkbook()
  
  theme_nature <- theme_classic(base_size = 26) + 
    theme(
      legend.position = "none", 
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
      axis.title = element_text(size = 28, face = "plain", colour = "black"), 
      axis.text = element_text(size = 22, colour = "black"),
      axis.line = element_line(linewidth = 1.0, colour = "black"),
      axis.ticks = element_line(linewidth = 1.0, colour = "black")
    )
  
  # --- PANEL COMPILATION ---
  for (t_idx in seq_along(analysis_tasks)) {
    task <- analysis_tasks[[t_idx]]
    curr_x <- if(is.null(task$x_col_override)) conf$gc_col else task$x_col_override
    
    # Model layers run on full dataset 
    raw_global <- cor.test(df_all[[curr_x]], df_all[[task$y_col]], method = "spearman", exact = FALSE)
    pic_x_all  <- pic(setNames(df_all[[curr_x]], rownames(df_all)), tree_all)
    pic_y_all  <- pic(setNames(df_all[[task$y_col]], rownames(df_all)), tree_all)
    
    df_pic_plot <- data.frame(px = pic_x_all, py = pic_y_all, NodeNames = names(pic_x_all)) %>%
      filter(NodeNames != "Ficedula_albicollis")
    
    deviation_summary <- rbind(deviation_summary, data.frame(
      Analysis = task$short, Clade = "Global", Spearman_Rho = round(raw_global$estimate, 3), 
      Spearman_P = raw_global$p.value, PIC_AdjR2 = round(summary(lm(pic_y_all ~ pic_x_all + 0))$adj.r.squared, 3), 
      Pagel_Lambda = round(phylosig(tree_all, setNames(df_all[[task$y_col]], rownames(df_all)), method="lambda")$lambda, 3)
    ))
    
    sheet_name <- substr(task$short, 1, 31)
    addWorksheet(wb_stats, sheet_name)
    writeData(wb_stats, sheet_name, data.frame(Clade="Global_Dataset", Rho=raw_global$estimate, P=raw_global$p.value), startRow=1)
    
    clade_pic_trends <- data.frame()
    
    for(grp in sub_clades) {
      df_sub <- df_all[df_all$Classification_1 == grp, ]
      if(nrow(df_sub) > 3) {
        tree_sub <- keep.tip(tree_all, rownames(df_sub))
        raw_sub  <- cor.test(df_sub[[curr_x]], df_sub[[task$y_col]], method = "spearman", exact = FALSE)
        px_sub   <- pic(setNames(df_sub[[curr_x]], rownames(df_sub)), tree_sub)
        py_sub   <- pic(setNames(df_sub[[task$y_col]], rownames(df_sub)), tree_sub)
        
        clade_pic_trends <- rbind(clade_pic_trends, data.frame(px = px_sub, py = py_sub, Clade = grp))
        
        deviation_summary <- rbind(deviation_summary, data.frame(
          Analysis = task$short, Clade = grp, Spearman_Rho = round(raw_sub$estimate, 3), 
          Spearman_P = raw_sub$p.value, PIC_AdjR2 = round(summary(lm(py_sub ~ px_sub + 0))$adj.r.squared, 3), 
          Pagel_Lambda = round(phylosig(tree_sub, setNames(df_sub[[task$y_col]], rownames(df_sub)), method="lambda")$lambda, 3)
        ))
      }
    }
    
    # 1. VISUAL LAYER: CORRELATION SPACE (UPDATED: Global is dashed, Clades are solid)
    p_corr <- ggplot(df_plot, aes_string(x = curr_x, y = task$y_col)) +
      geom_point(aes(color = Classification_1), size = 4.5, alpha = 0.40, stroke = 0) +
      # Global representation line changed to a dashed vector line
      geom_smooth(data = df_all, aes_string(x = curr_x, y = task$y_col), method = "lm", 
                  color = "#333333", linetype = "dashed", linewidth = 2.5, se = FALSE) +
      # Line paths for sub-clades changed to bold solid vector lines
      geom_smooth(data = filter(df_all, Classification_1 %in% sub_clades), 
                  aes_string(x = curr_x, y = task$y_col, color = "Classification_1", group = "Classification_1"), 
                  method = "lm", linetype = "solid", linewidth = 2.8, se = FALSE) +
      scale_color_manual(values = custom_colors) + 
      theme_nature +
      labs(x = task$x_label, y = task$y_label)
    
    # 2. VISUAL LAYER: PIC SPACE (Remains unchanged: Global solid, Clades dashed)
    p_pic <- ggplot(df_pic_plot, aes(x = px, y = py)) +
      geom_point(alpha = 0.30, color = "#555555", size = 4.0, stroke = 0) +
      geom_smooth(data = data.frame(px = pic_x_all, py = pic_y_all), aes(x = px, y = py),
                  method = "lm", formula = y ~ x + 0, color = "#333333", linetype = "solid", linewidth = 2.5, se = FALSE) +
      geom_smooth(data = clade_pic_trends, aes(x = px, y = py, color = Clade, group = Clade), 
                  method = "lm", formula = y ~ x + 0, linetype = "dashed", linewidth = 2.8, se = FALSE) +
      scale_color_manual(values = custom_colors) +
      theme_nature + 
      labs(x = paste0("Contrast: ", task$x_label), y = paste0("Contrast: ", task$y_label))
    
    all_panel_plots[[length(all_panel_plots)+1]] <- p_corr
    all_panel_plots[[length(all_panel_plots)+1]] <- p_pic
  }
  
  saveWorkbook(wb_stats, file.path(tab_dir, paste0("Stats_Calculations_", data_name, ".xlsx")), overwrite = TRUE)
  
  # --- 4. EXCLUSIVE EDITORIAL LEGEND (Updated to match Raw visual structures)
  dummy_plot <- ggplot(data.frame(Clade=factor(clade_order, levels=clade_order), x=1, y=1), aes(x, y, color=Clade)) +
    geom_smooth(data = data.frame(Clade=factor(c("Mammals", "Aves", "Reptiles", "Non-amniotes"), levels=clade_order), x=1, y=1),
                method="lm", linewidth=3, linetype="solid", key_glyph = draw_key_timeseries, se=FALSE) + 
    geom_smooth(data = data.frame(Clade=factor("Others", levels=clade_order), x=1, y=1),
                method="lm", linewidth=3, linetype="dashed", color="#333333", key_glyph = draw_key_timeseries, se=FALSE) +
    scale_color_manual(values=custom_colors) +
    theme_minimal() + 
    theme(
      legend.position="bottom", 
      legend.title=element_text(size=24, colour="black"),
      legend.text=element_text(size=22, colour="black"),
      legend.key.width = unit(2.8, "cm"),
      legend.spacing.x = unit(0.6, "cm")
    ) + 
    labs(color="Lineage Cohort:")
  m_legend <- get_legend(dummy_plot)
  
  # --- 5. GRID MATRIX CONFIGURATION & COHESION ---
  assembled_matrix <- wrap_plots(all_panel_plots, ncol = 6, nrow = 1)
  master_layout_canvas <- plot_grid(assembled_matrix, m_legend, ncol = 1, rel_heights = c(1, 0.08))
  
  file_base_path <- file.path(fig_dir, paste0("Figure_", data_name, "_Overlay_Master"))
  
  ggsave(paste0(file_base_path, ".pdf"), master_layout_canvas, width = 42, height = 8.5, device = "pdf", useDingbats = FALSE)
  ggsave(paste0(file_base_path, ".png"), master_layout_canvas, width = 42, height = 8.5, dpi = 300)
  
  write.xlsx(deviation_summary, file.path(tab_dir, paste0("Summary_Statistical_Output_", data_name, ".xlsx")))
}

cat("\n*** Success! Inverted line structures applied to Raw spaces, PIC spaces preserved safely. ***\n")