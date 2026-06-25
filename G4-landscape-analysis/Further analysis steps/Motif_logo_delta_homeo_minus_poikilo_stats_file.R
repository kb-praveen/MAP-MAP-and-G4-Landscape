# --- 0. LOAD LIBRARIES ---
suppressPackageStartupMessages({
  library(ggseqlogo)
  library(universalmotif)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(grid)
})

# --- 1. SETUP WORKING DIRECTORIES & DATA LOOKUPS ---
base_path   <- "Path/G4/G4_Final_Analysis_Complete/PWMS/MEME_PWMs"
test_dir    <- "Path/G4/G4_Final_Analysis_Complete/test"

out_bg_sub     <- "Path/G4/G4_Final_Analysis_Complete/Subtracted_PWMs/Background_subtracted"
out_homeo_poik <- "Path/G4/G4_Final_Analysis_Complete/Subtracted_PWMs/Homeo_minus_poikilo"

enrichments <- c("Highly_Enriched", "Highly_Reduced", "Moderately_Enriched", "Moderately_Reduced")
regions     <- c("Genome", "Promoter")

dir.create(test_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_bg_sub, showWarnings = FALSE, recursive = TRUE)
dir.create(out_homeo_poik, showWarnings = FALSE, recursive = TRUE)

stats_promoter_path <- "Path/Stats_Promoter_Homeotherms_vs_Poikilotherms.tsv"
stats_genome_path   <- "Path/Stats_Genome_Homeotherms_vs_Poikilotherms.tsv"

stats_promoter <- if(file.exists(stats_promoter_path)) read.delim(stats_promoter_path, sep="\t", stringsAsFactors=FALSE) else data.frame()
stats_genome   <- if(file.exists(stats_genome_path)) read.delim(stats_genome_path, sep="\t", stringsAsFactors=FALSE) else data.frame()

# --- REVERSE COMPLEMENT HELPER FOR C-RICH STRAND ENFORCEMENT ---
reverse_complement_matrix <- function(mat) {
  rownames_ordered <- c("A", "C", "G", "T")
  rc_mat <- mat[c("T", "G", "C", "A"), , drop = FALSE]
  rownames(rc_mat) <- rownames_ordered
  rc_mat <- rc_mat[, rev(seq_len(ncol(rc_mat))), drop = FALSE]
  return(rc_mat)
}

# --- CUSTOM ENGINE FOR TRUE PROPORTIONAL DIFFERENCE BAR LOGO ---
plot_true_difference_bars <- function(mat_h, mat_p, title_str, subtitle_str) {
  num_cols <- ncol(mat_h)
  color_map <- c("A" = "#109618", "C" = "#3366CC", "G" = "#FF9900", "T" = "#DC3912")
  
  diff_mat <- mat_h - mat_p
  colnames(diff_mat) <- paste0("PosColumn_", 1:num_cols)
  
  diff_df <- as.data.frame(diff_mat) %>%
    mutate(Base = rownames(diff_mat)) %>%
    pivot_longer(
      cols = starts_with("PosColumn_"), 
      names_to = "Position_Raw", 
      values_to = "Difference"
    ) %>%
    mutate(Position = as.numeric(gsub("PosColumn_", "", Position_Raw)))
  
  max_val <- max(abs(diff_df$Difference), na.rm = TRUE)
  y_lim <- max(0.05, ceiling(max_val * 20) / 20) 
  
  p_bars <- ggplot(diff_df, aes(x = factor(Position), y = Difference, fill = Base)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.8, color = NA) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
    scale_fill_manual(values = color_map) +
    # FIXED: Hardcoded exact coordinate axis bounds to achieve perfect vertical alignment with top panel
    scale_x_discrete(expand = c(0.05, 0.5)) +
    scale_y_continuous(limits = c(-y_lim, y_lim), expand = c(0, 0)) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "grey80", linewidth = 0.6, linetype = "dashed"),
      strip.text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10)
    ) +
    labs(
      x    = "Position", 
      y    = expression(Delta ~ "Frequency"), 
      fill = "Base"
    )
  
  return(p_bars)
}

# --- 2. MAIN VISUAL COMPILATION ENGINE LOOP ---
for (reg in regions) {
  cat("\n============================================\n")
  cat("[INFO] Commencing Pipeline Target Region:", reg, "\n")
  cat("============================================\n")
  
  bg_sub_reg_dir <- file.path(out_bg_sub, reg)
  hm_pk_reg_dir  <- file.path(out_homeo_poik, reg)
  
  dir.create(bg_sub_reg_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(hm_pk_reg_dir, showWarnings = FALSE, recursive = TRUE)
  
  current_meta_df <- if(reg == "Promoter") stats_promoter else stats_genome
  
  for (enr in enrichments) {
    ref_folder <- file.path(base_path, "Homeotherm", enr, paste0(reg, "_G4"))
    
    if (!dir.exists(ref_folder)) {
      cat("[SKIP] Folder missing:", ref_folder, "\n")
      next
    }
    
    enr_clean_text    <- gsub("_", " ", enr)
    full_title_header <- paste0(reg, " : ", enr_clean_text)
    
    motif_files <- list.files(ref_folder, pattern = "\\.meme$")
    
    for (f in motif_files) {
      motif_id_base <- gsub("\\.meme$", "", f)
      
      paths <- list(
        "Homeo_G4"      = file.path(base_path, "Homeotherm", enr, paste0(reg, "_G4"), f),
        "Homeo_Free"    = file.path(base_path, "Homeotherm", enr, paste0(reg, "_Free"), f),
        "Poikilo_G4"    = file.path(base_path, "Poikilotherm", enr, paste0(reg, "_G4"), f),
        "Poikilo_Free"  = file.path(base_path, "Poikilotherm", enr, paste0(reg, "_Free"), f)
      )
      
      if (!all(sapply(paths, file.exists))) next
      
      match_row <- current_meta_df[current_meta_df$Motif_ID == motif_id_base, ]
      meta <- if (nrow(match_row) > 0) {
        list(id = match_row$Motif_ID[1], name = match_row$Motif_Name[1])
      } else {
        list(id = motif_id_base, name = motif_id_base)
      }
      
      umotifs <- tryCatch({ 
        lapply(names(paths), function(n) { m <- read_meme(paths[[n]]); m@name <- n; m }) 
      }, error = function(e) return(NULL))
      
      if(is.null(umotifs)) next
      names(umotifs) <- names(paths)
      
      mat_h  <- as.matrix(umotifs[["Homeo_G4"]]@motif)
      mat_hf <- as.matrix(umotifs[["Homeo_Free"]]@motif)
      mat_p  <- as.matrix(umotifs[["Poikilo_G4"]]@motif)
      mat_pf <- as.matrix(umotifs[["Poikilo_Free"]]@motif)
      
      if (ncol(mat_h) != ncol(mat_hf) || ncol(mat_h) != ncol(mat_p) || ncol(mat_h) != ncol(mat_pf)) {
        cat("[WARN] Dimension mismatch for motif:", meta$name, "- skipping processing.\n")
        next
      }
      
      # --- STRAND RULE: ENFORCE C-RICH STRAND REPRESENTATION ---
      if (sum(mat_h["G", ]) > sum(mat_h["C", ])) {
        mat_h  <- reverse_complement_matrix(mat_h)
        mat_hf <- reverse_complement_matrix(mat_hf)
        mat_p  <- reverse_complement_matrix(mat_p)
        mat_pf <- reverse_complement_matrix(mat_pf)
      }
      
      # ------------------------------------------------------------------------
      # --- TRACK 1: BACKGROUND SUBTRACTED PLOTS (G4 - FREE BACKGROUND) ---
      # ------------------------------------------------------------------------
      bg_sub_homeo   <- pmax(mat_h - mat_hf, 0)
      bg_sub_poikilo <- pmax(mat_p - mat_pf, 0)
      
      if (sum(bg_sub_homeo) == 0 || sum(bg_sub_poikilo) == 0) {
        cat(sprintf("[SKIP BG LOGO] Subtracted matrix flat zero for %s (%s)\n", meta$name, meta$id))
      } else {
        list_bg_suburbed <- list(
          "Homeotherm (pG4-background)"   = bg_sub_homeo,
          "Poikilotherm (pG4-background)" = bg_sub_poikilo
        )
        
        p_bg_sub <- ggseqlogo(list_bg_suburbed, ncol = 1) +
          theme_bw() +
          theme(
            strip.text = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 12),
            axis.text  = element_text(size = 10),
            plot.title = element_text(size = 20, face = "bold"),
            plot.subtitle = element_text(size = 16, face = "italic")
          ) +
          labs(title = full_title_header, subtitle = paste(meta$name, "-", meta$id), y = "Bits")
        
        ggsave(
          filename = file.path(bg_sub_reg_dir, paste0(meta$id, "_", enr, "_BackgroundSubtracted.png")),
          plot     = p_bg_sub, width = 11, height = 7, dpi = 300
        )
      }
      
      # ------------------------------------------------------------------------
      # --- TRACK 2: HOMEOTHERM VS POIKILOTHERM G4 (TRUE PROPORTIONAL BARS) ---
      # ------------------------------------------------------------------------
      list_raws <- list("Homeotherm (pG4)" = mat_h, "Poikilotherm (pG4)" = mat_p)
      p_raw_logos <- ggseqlogo(list_raws, ncol = 1) + 
        theme_bw() + 
        theme(strip.text = element_text(size = 14, face = "bold"), 
              axis.title = element_text(size = 12), axis.text = element_text(size = 10)) + labs(y = "Bits")
      
      p_true_bars <- plot_true_difference_bars(mat_h, mat_p, full_title_header, paste(meta$name, "-", meta$id))
      p_true_bars <- p_true_bars + facet_wrap(~ "pG4 (homeotherm - poikilotherm)")
      
      p_assembled_thermal <- p_raw_logos / p_true_bars + 
        plot_layout(heights = c(2, 1)) +
        plot_annotation(
          title    = full_title_header,
          subtitle = paste(meta$name, "-", meta$id),
          theme    = theme(plot.title = element_text(size = 20, face = "bold"),
                           plot.subtitle = element_text(size = 16, face = "italic"))
        )
      
      p_raw_logos$labels$title <- NULL; p_raw_logos$labels$subtitle <- NULL
      p_true_bars$labels$title <- NULL; p_true_bars$labels$subtitle <- NULL
      
      ggsave(
        filename = file.path(hm_pk_reg_dir, paste0(meta$id, "_", enr, "_Homeo_minus_Poikilo.png")),
        plot     = p_assembled_thermal, width = 12, height = 10, dpi = 300
      )
      
      # --- POSITIONAL GC TRACKING METRICS GENERATION PASS ---
      for (i in 1:ncol(mat_h)) {
        h_gc  <- sum(mat_h[c("G", "C"), i])
        hf_gc <- sum(mat_hf[c("G", "C"), i])
        p_gc  <- sum(mat_p[c("G", "C"), i])
        pf_gc <- sum(mat_pf[c("G", "C"), i])
        
        global_locational_list[[length(global_locational_list)+1]] <- data.frame(
          Region=reg, Category=full_title_header, Motif_ID=meta$id, Motif_Name=meta$name, Position=i,
          Homeotherm_G4=h_gc, Homeotherm_G4_free=hf_gc, 
          Poikilotherm_G4=p_gc, Poikilotherm_G4_free=pf_gc, 
          stringsAsFactors=FALSE
        )
      }
      
      # --- SUMMARY BLOCK APPEND STEP ---
      global_motif_summary[[length(global_motif_summary)+1]] <- data.frame(
        Region=reg, Category=full_title_header, Motif_ID=meta$id, Motif_Name=meta$name,
        Homeotherm_Site_Count=umotifs[["Homeo_G4"]]@nsites, 
        Poikilotherm_Site_Count=umotifs[["Poikilo_G4"]]@nsites,
        Homeotherm_G4_Mean_GC=mean(colSums(mat_h[c("G", "C"),, drop=FALSE])),
        Poikilotherm_G4_Mean_GC=mean(colSums(mat_p[c("G", "C"),, drop=FALSE])),
        Homeotherm_G4_Free_Mean_GC=mean(colSums(mat_hf[c("G", "C"),, drop=FALSE])),
        Poikilotherm_G4_Free_Mean_GC=mean(colSums(mat_pf[c("G", "C"),, drop=FALSE])),
        stringsAsFactors=FALSE
      )
    }
  }
}

# --- 3. MASTER DATA METRICS EXPORT ---
if(length(global_motif_summary) > 0) {
  cat("\n[INFO] Saving compilation data out to disk...\n")
  motif_master <- bind_rows(global_motif_summary)
  loc_master   <- bind_rows(global_locational_list)
  write.table(motif_master, file.path(test_dir, "Overall_Motif_Enrichment_Summary.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
  write.table(loc_master, file.path(test_dir, "Locational_Subtracted_GC_Distribution.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
  cat("\n*** 📊 PIPELINE COMPLETE: HIGH-CONTRAST SEPARATION SEGREGATIONS GENERATED SUCCESSFULLY ***\n")
}