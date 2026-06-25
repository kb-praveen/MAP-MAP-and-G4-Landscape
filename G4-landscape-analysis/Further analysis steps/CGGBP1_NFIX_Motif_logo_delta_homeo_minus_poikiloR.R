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
jaspar_dir  <- "Path/G4/G4_Final_Analysis_Complete/PWMS/JASPAR2026meme"
output_dir  <- "Path/G4/G4_Final_Analysis_Complete/Subtracted_PWMs/Target_Filtered_Promoters"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

stats_promoter_path <- "Path/Stats_Promoter_Homeotherms_vs_Poikilotherms.tsv"
stats_promoter      <- if(file.exists(stats_promoter_path)) read.delim(stats_promoter_path, sep="\t", stringsAsFactors=FALSE) else data.frame()

target_motifs <- c("CGGBP1", "NFIX")
enrichments   <- c("Highly_Enriched", "Highly_Reduced", "Moderately_Enriched", "Moderately_Reduced")

# --- REVERSE COMPLEMENT HELPER FOR C-RICH STRAND ENFORCEMENT ---
reverse_complement_matrix <- function(mat) {
  rownames_ordered <- c("A", "C", "G", "T")
  rc_mat <- mat[c("T", "G", "C", "A"), , drop = FALSE]
  rownames(rc_mat) <- rownames_ordered
  rc_mat <- rc_mat[, rev(seq_len(ncol(rc_mat))), drop = FALSE]
  return(rc_mat)
}

# --- 2. MAIN VISUAL COMPILATION ENGINE LOOP ---
cat("\n============================================\n")
cat("[INFO] Running Uniform Font Size Customization Pipeline\n")
cat("============================================\n")

filtered_stats <- stats_promoter %>% 
  filter(toupper(Motif_Name) %in% target_motifs)

if(nrow(filtered_stats) == 0) {
  stop("[CRITICAL] Target motifs not found in promoter stats table metadata file.")
}

for (enr in enrichments) {
  ref_folder <- file.path(base_path, "Homeotherm", enr, "Promoter_G4")
  if (!dir.exists(ref_folder)) next
  
  for (i in 1:nrow(filtered_stats)) {
    target_id   <- filtered_stats$Motif_ID[i]
    target_name <- filtered_stats$Motif_Name[i]
    
    meme_file <- paste0(target_id, ".meme")
    homeo_g4_path   <- file.path(base_path, "Homeotherm", enr, "Promoter_G4", meme_file)
    poikilo_g4_path <- file.path(base_path, "Poikilotherm", enr, "Promoter_G4", meme_file)
    
    if (!file.exists(homeo_g4_path) || !file.exists(poikilo_g4_path)) next
    
    jaspar_match <- list.files(jaspar_dir, pattern = paste0("^", target_id, "(\\.[0-9]+)?\\.meme$"), full.names = TRUE)
    if (length(jaspar_match) == 0) next
    
    m_h <- read_meme(homeo_g4_path)
    m_p <- read_meme(poikilo_g4_path)
    m_j <- read_meme(jaspar_match[1])
    
    mat_h <- as.matrix(m_h@motif)
    mat_p <- as.matrix(m_p@motif)
    mat_j <- as.matrix(m_j@motif)
    
    if (ncol(mat_h) != ncol(mat_p) || ncol(mat_h) != ncol(mat_j)) next
    
    # --- STRAND ENFORCEMENT ENGINE RULE ---
    if (sum(mat_h["G", ]) > sum(mat_h["C", ])) {
      mat_h <- reverse_complement_matrix(mat_h)
      mat_p <- reverse_complement_matrix(mat_p)
      mat_j <- reverse_complement_matrix(mat_j)
    }
    
    num_cols   <- ncol(mat_h)
    pos_labels <- as.character(1:num_cols)
    color_map  <- c("A" = "#109618", "C" = "#3366CC", "G" = "#FF9900", "T" = "#DC3912")
    
    # ------------------------------------------------------------------------
    # --- VISUAL PANEL 1: STREAMLINED JASPAR LOGO (UNIFORM FONT) ---
    # ------------------------------------------------------------------------
    colnames(mat_j) <- pos_labels
    
    p_jaspar <- ggseqlogo(mat_j) + 
      scale_x_discrete(limits = pos_labels, expand = c(0.05, 0.5)) +
      theme_classic() + 
      theme(
        axis.line.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.x = element_blank(),
        # UNIFORM FONT SIZE CONFIGURATION (18pt Across Everything)
        axis.title.y = element_text(size = 22, face = "plain", color = "black", lineheight = 1.2),
        axis.text.y  = element_text(size = 22, color = "black"),
        plot.margin  = margin(b = 12, t = 10, r = 10, l = 10)
      ) + 
      labs(y = paste0(target_name, "\n(JASPAR: ", target_id, ")"))
    
    # ------------------------------------------------------------------------
    # --- VISUAL PANEL 2: STREAMLINED DELTA FREQUENCY BARS (UNIFORM FONT) ---
    # ------------------------------------------------------------------------
    diff_mat <- mat_h - mat_p
    colnames(diff_mat) <- paste0("Pos_", pos_labels)
    
    diff_df <- as.data.frame(diff_mat) %>%
      mutate(Base = rownames(diff_mat)) %>%
      pivot_longer(cols = starts_with("Pos_"), names_to = "Pos_Raw", values_to = "Difference") %>%
      mutate(Position = gsub("Pos_", "", Pos_Raw))
    
    max_val <- max(abs(diff_df$Difference), na.rm = TRUE)
    y_lim   <- max(0.05, ceiling(max_val * 20) / 20) 
    
    p_bars <- ggplot(diff_df, aes(x = factor(Position, levels = pos_labels), y = Difference, fill = Base)) +
      geom_col(position = position_dodge(width = 0.85), width = 0.8, color = NA) +
      geom_hline(yintercept = 0, color = "black", linewidth = 1.0) +
      scale_fill_manual(values = color_map) +
      scale_x_discrete(limits = pos_labels, expand = c(0.05, 0.5)) +
      scale_y_continuous(limits = c(-y_lim, y_lim), expand = c(0, 0)) +
      theme_classic() + 
      theme(
        # UNIFORM FONT SIZE CONFIGURATION (18pt Across Everything)
        axis.title.x = element_text(size = 22, face = "plain", color = "black", margin = margin(t = 12)),
        axis.text.x  = element_text(size = 22, face = "plain", color = "black"),
        axis.title.y = element_text(size = 22, face = "plain", color = "black", lineheight = 1.2),
        axis.text.y  = element_text(size = 22, color = "black"),
        legend.position = "right",
        legend.title = element_text(size = 22, face = "plain"),
        legend.text  = element_text(size = 22),
        plot.margin  = margin(t = 12, b = 10, r = 10, l = 10)
      ) +
      labs(
        x    = "Position", 
        y    = "Δ Frequency\n(Homeotherms - Poikilotherms)", 
        fill = "Base"
      )
    
    # ------------------------------------------------------------------------
    # --- VISUAL COMPILATION ASSEMBLY PASS VIA PATCHWORK ---
    # ------------------------------------------------------------------------
    p_assembled <- p_jaspar / plot_spacer() / p_bars + 
      plot_layout(heights = c(0.7, 0.06, 1.3))
    
    # Save the uniform high-resolution layout out to disk
    output_filename <- sprintf("%s_%s_Promoter_Streamlined_Alignment.png", target_name, enr)
    ggsave(
      filename = file.path(output_dir, output_filename),
      plot     = p_assembled, width = 11, height = 9, dpi = 300
    )
    
    cat(sprintf("[SUCCESS] Exported uniform-font chart for: %s (%s) [%s]\n", target_name, target_id, enr))
  }
}

cat("\n*** 📊 PIPELINE COMPLETE: ALL FONT SIZES EQUALIZED TO 18PT PERFECTLY ***\n")