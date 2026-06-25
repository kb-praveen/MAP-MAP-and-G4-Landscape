# --- 0. LOAD LIBRARIES ---
suppressPackageStartupMessages({
  library(ggseqlogo)
  library(universalmotif)
  library(ggplot2)
  library(ggdendro)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

# --- 1. SETUP ---
# Corrected local paths
base_path   <- "/Users/praveenbishnoi/Desktop/PWMS/MEME_PWMs"
jaspar_dir  <- "/Users/praveenbishnoi/Desktop/PWMS/JASPAR2026meme"
output_root <- "/Users/praveenbishnoi/Desktop/G4_Final_Analysis_Complete"
test_dir    <- file.path(output_root, "test")

if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)
if (!dir.exists(test_dir)) dir.create(test_dir, recursive = TRUE)

# Updated nomenclature matching the Python Script output
enrichments <- c("Highly_Enriched", "Highly_Reduced", 
                 "Moderately_Enriched", "Moderately_Reduced")
regions <- c("Genome", "Promoter")

global_locational_list <- list()
global_motif_summary   <- list()

get_jaspar_metadata <- function(filepath) {
  lines <- readLines(filepath, n = 20)
  motif_line <- lines[grep("^MOTIF", lines)]
  if(length(motif_line) > 0) {
    parts <- strsplit(motif_line, "\\s+")[[1]]
    return(list(id = parts[2], name = ifelse(length(parts) >= 3, parts[3], parts[2])))
  }
  return(list(id = "Unknown", name = "Unknown"))
}

# --- 2. MAIN PROCESSING LOOP ---
for (reg in regions) {
  cat("\nProcessing Region:", reg, "\n")
  for (enr in enrichments) {
    
    # Matching the structure: Thermal_Group / Enrichment_Tier / Category
    # Category is e.g., "Promoter_G4" or "Genome_G4"
    ref_folder <- file.path(base_path, "Homeotherm", enr, paste0(reg, "_G4"))
    
    if (!dir.exists(ref_folder)) {
      cat("Skipping missing folder:", ref_folder, "\n")
      next
    }
    
    enr_clean_text <- gsub("_", " ", enr)
    full_title_header <- paste0(reg, " : ", enr_clean_text, " in G4s")
    
    current_out_dir <- file.path(output_root, reg, paste0("Homeotherms_", tolower(enr)))
    if (!dir.exists(current_out_dir)) dir.create(current_out_dir, recursive = TRUE)
    
    motif_files <- list.files(ref_folder, pattern = "\\.meme$")
    
    for (f in motif_files) {
      motif_id_base <- gsub("\\.meme$", "", f)
      # Match JASPAR reference (allows for version numbers like .1, .2)
      jaspar_match <- list.files(jaspar_dir, pattern = paste0("^", motif_id_base, "(\\.[0-9]+)?\\.meme$"), full.names = TRUE)
      
      if (length(jaspar_match) == 0) next
      
      meta <- get_jaspar_metadata(jaspar_match[1])
      
      # Build absolute paths for all 5 required matrices
      paths <- list(
        "JASPAR Reference"      = jaspar_match[1],
        "Homeotherm (G4)"        = file.path(base_path, "Homeotherm", enr, paste0(reg, "_G4"), f),
        "Homeotherm (G4 free)"  = file.path(base_path, "Homeotherm", enr, paste0(reg, "_Free"), f),
        "Poikilotherm (G4)"     = file.path(base_path, "Poikilotherm", enr, paste0(reg, "_G4"), f),
        "Poikilotherm (G4 free)"= file.path(base_path, "Poikilotherm", enr, paste0(reg, "_Free"), f)
      )
      
      if (!all(sapply(paths, file.exists))) next
      
      umotifs <- tryCatch({ 
        lapply(names(paths), function(n) { m <- read_meme(paths[[n]]); m@name <- n; m }) 
      }, error = function(e) return(NULL))
      
      if(is.null(umotifs)) next
      names(umotifs) <- names(paths)
      
      m_h <- umotifs[["Homeotherm (G4)"]]; m_hf <- umotifs[["Homeotherm (G4 free)"]]
      m_p <- umotifs[["Poikilotherm (G4)"]]; m_pf <- umotifs[["Poikilotherm (G4 free)"]]
      m_r <- umotifs[["JASPAR Reference"]]
      
      # --- Positional GC Calculations ---
      motif_gc_diffs <- numeric(ncol(m_h@motif))
      for (i in 1:ncol(m_h@motif)) {
        h_gc  <- sum(m_h@motif[c("G", "C"), i])
        hf_gc <- sum(m_hf@motif[c("G", "C"), i])
        p_gc  <- sum(m_p@motif[c("G", "C"), i])
        pf_gc <- sum(m_pf@motif[c("G", "C"), i])
        r_gc  <- sum(m_r@motif[c("G", "C"), i])
        
        motif_gc_diffs[i] <- h_gc - r_gc # Retention metric
        
        global_locational_list[[length(global_locational_list)+1]] <- data.frame(
          Region=reg, Category=full_title_header, Motif_ID=meta$id, Motif_Name=meta$name, Position=i,
          Homeotherm_G4=h_gc, Homeotherm_G4_free=hf_gc, Poikilotherm_G4=p_gc, Poikilotherm_G4_free=pf_gc, 
          JASPAR_Reference=r_gc, Retention_Homeo_G4=h_gc-r_gc, stringsAsFactors=F)
      }
      
      # --- Summary Stats ---
      global_motif_summary[[length(global_motif_summary)+1]] <- data.frame(
        Region=reg, Category=full_title_header, Motif_ID=meta$id, Motif_Name=meta$name,
        Homeotherm_Site_Count=m_h@nsites, Poikilotherm_Site_Count=m_p@nsites,
        Homeotherm_G4_GC=mean(colSums(m_h@motif[c("G", "C"),, drop=FALSE])),
        Poikilotherm_G4_GC=mean(colSums(m_p@motif[c("G", "C"),, drop=FALSE])),
        Homeotherm_G4_Free_GC=mean(colSums(m_hf@motif[c("G", "C"),, drop=FALSE])),
        Poikilotherm_G4_Free_GC=mean(colSums(m_pf@motif[c("G", "C"),, drop=FALSE])),
        JASPAR_Reference_GC=mean(colSums(m_r@motif[c("G", "C"),, drop=FALSE])), stringsAsFactors=F
      )
      
      # --- 3. LOGO PLOTTING ---
      highlights <- which(motif_gc_diffs > 0.10) # Positions with >10% GC increase vs JASPAR
      sim_mat <- compare_motifs(umotifs, method = "PCC", tryRC = TRUE)
      hc <- hclust(as.dist(1 - sim_mat), method = "average")
      
      p_tree <- ggplot(segment(dendro_data(hc))) + 
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 1.2) +
        coord_flip() + scale_y_reverse() + theme_dendro()
      
      logo_list <- lapply(hc$labels[hc$order], function(n) as.matrix(umotifs[[n]]@motif))
      names(logo_list) <- hc$labels[hc$order]
      
      p_logos <- ggseqlogo(logo_list, ncol = 1) + 
        annotate('rect', xmin = highlights - 0.5, xmax = highlights + 0.5, ymin = 0, ymax = 2, alpha = .15, fill = "red") + 
        theme_bw() + theme(strip.text = element_text(size = 14, face = "bold"),
                           axis.title = element_text(size = 12),
                           axis.text = element_text(size = 10)) + labs(y = "Bits")
      
      p_final <- p_tree + p_logos + 
        plot_layout(widths = c(1, 4)) +
        plot_annotation(title = full_title_header,
                        subtitle = paste(meta$id, "-", meta$name),
                        theme = theme(plot.title = element_text(size = 20, face = "bold"),
                                      plot.subtitle = element_text(size = 16, face = "italic")))
      
      ggsave(file.path(current_out_dir, paste0(meta$id, "_Clustered_Grid.png")), p_final, width = 14, height = 12, dpi = 300)
    }
  }
}

# --- 4. MASTER FILE GENERATION ---
if(length(global_motif_summary) > 0) {
  motif_master <- bind_rows(global_motif_summary)
  loc_master   <- bind_rows(global_locational_list)
  
  # Proportion test helper
  calc_pval <- function(gc_obs, n_sites, gc_ref) {
    if(is.na(n_sites) || n_sites < 10 || is.na(gc_obs) || is.na(gc_ref)) return(NA)
    tryCatch({ prop.test(x = round(gc_obs * n_sites), n = n_sites, p = gc_ref)$p.value }, error = function(e) return(NA))
  }
  
  master_sig <- motif_master %>%
    mutate(
      H_Ret_Pct = (Homeotherm_G4_GC - JASPAR_Reference_GC) * 100,
      P_Ret_Pct = (Poikilotherm_G4_GC - JASPAR_Reference_GC) * 100,
      Abs_Divergence = abs(H_Ret_Pct - P_Ret_Pct),
      Target_GC = ifelse(abs(H_Ret_Pct) >= abs(P_Ret_Pct), Homeotherm_G4_GC, Poikilotherm_G4_GC),
      Target_N  = ifelse(abs(H_Ret_Pct) >= abs(P_Ret_Pct), Homeotherm_Site_Count, Poikilotherm_Site_Count)
    ) %>%
    rowwise() %>%
    mutate(Pval = calc_pval(Target_GC, Target_N, JASPAR_Reference_GC)) %>%
    ungroup() %>%
    mutate(Adj_P_val_FDR = p.adjust(Pval, method = "fdr")) %>%
    filter(Adj_P_val_FDR < 0.05)
  
  # Save results
  write.table(master_sig, file.path(test_dir, "Master_Significant_Motifs_Full_Data.tsv"), sep="\t", row.names=F, quote=F)
  write.table(loc_master, file.path(output_root, "Master_Locational_GC_Retention.tsv"), sep="\t", row.names=F, quote=F)
  write.table(motif_master, file.path(output_root, "Overall_Motif_Enrichment_Summary.tsv"), sep="\t", row.names=F, quote=F)
  
  cat("\n*** PIPELINE COMPLETE ***\n")
} else {
  cat("\n*** ERROR: No motifs processed. Check your directory structure and MEME file names. ***\n")
}