# 1. LOAD REQUIRED LIBRARIES
suppressPackageStartupMessages({
  library(readxl)         # For reading the multi-sheet summary workbook
  library(ComplexHeatmap) # For advanced publication heatmaps
  library(circlize)       # For dynamic color scales
  library(dendextend)     # For branch width line thickness customizations
  library(grid)           # For graphical metric tracking parameters
  library(tidyverse)
})

# --- SETUP WORKING DIRECTORIES ---
promoter_dir <- "Path/Promoter"
input_file   <- file.path(promoter_dir, "Final_all_mutations_summary_promoter.xlsx")

# Destinations for the optimized horizontal plots
output_master <- file.path(promoter_dir, "CtoT_horizontal_heatmap_1_all_motifs.png")
output_subset <- file.path(promoter_dir, "CtoT_horizontal_heatmap_2_selected_subset.png")
output_rest   <- file.path(promoter_dir, "CtoT_horizontal_heatmap_3_remaining_motifs.png")

# Separate standalone unique legend output files
legend_master_out <- file.path(promoter_dir, "Legend_Heatmap_1_All_Motifs.png")
legend_subset_out <- file.path(promoter_dir, "Legend_Heatmap_2_Filtered_Subset.png")
legend_rest_out   <- file.path(promoter_dir, "Legend_Heatmap_3_Remaining_Background.png")

# --- ADJUSTABLE PARAMETERS ---
dend_lwd       <- 4.0      # Thickened dendrogram lines for the sample rows and columns
sample_columns <- c("Lc", "Ac", "Gg", "Hs", "Ev", "Input")

# List of 22 explicit target motifs to isolate into Heatmap 2
target_motifs <- c(
  "CEBPE", "ZFP28", "ZBTB32", "PITX2", "OTX2", "LEUTX", "SNAI2", "SNAI1", 
  "SNAI3", "MITF", "USF2", "TFE3", "ELK1", "ETV4", "ETV5", "FOSL2", "PAX9", 
  "GCM1", "ZNF20", "HES1", "HES5", "Rarb"
)
# ----------------------------------

# --- 2. IMPORT DATA & CONSERVE MOTIF NAMES ---
if (!file.exists(input_file)) {
  stop(sprintf("[ERROR] Could not find required summary workbook at:\n%s", input_file))
}

raw_df <- read_excel(input_file, sheet = "C_to_T")
raw_df <- as.data.frame(raw_df)

rownames(raw_df) <- make.unique(if_else(is.na(raw_df$`Motif Name`) | raw_df$`Motif Name` == "", raw_df[[1]], raw_df$`Motif Name`))

# --- 3. CONVERT MATRIX NUMERIC VALUES ---
sub_df <- raw_df[, sample_columns, drop = FALSE]
sub_df[] <- lapply(sub_df, as.numeric)
matrix_master_raw <- data.matrix(sub_df)

# --- 4. SUBSET MATRICES GENERATION & TRANSPOSITION ---
all_row_names <- rownames(matrix_master_raw)
target_pattern <- paste0("^(", paste(target_motifs, collapse = "|"), ")(\\.|_|$)")
subset_indices <- grep(target_pattern, all_row_names, ignore.case = TRUE)

if (length(subset_indices) == 0) {
  stop("[ERROR] Could not find matching signatures inside the matrix for your specified subset list.")
}

matrix_master <- t(matrix_master_raw)
matrix_subset <- t(matrix_master_raw[subset_indices, , drop = FALSE])
matrix_rest   <- t(matrix_master_raw[-subset_indices, , drop = FALSE])

colnames(matrix_master) <- toupper(colnames(matrix_master))
colnames(matrix_subset) <- toupper(colnames(matrix_subset))
colnames(matrix_rest)   <- toupper(colnames(matrix_rest))


# --- 5. DEFINE DISTINCT HIGH-CONTRAST COLOR SCALES ---
max_val <- max(abs(matrix_master), na.rm = TRUE)
col_fun_master <- colorRamp2(
  c(-max_val, -max_val*0.33, 0, max_val*0.33, max_val), 
  c("#1A4475", "#6BAED6", "white", "#F4A582", "#B81D1D")
)

col_fun_subset <- colorRamp2(
  c(-0.4, -0.12, 0, 0.12, 0.4), 
  c("#1A4475", "#6BAED6", "white", "#F4A582", "#B81D1D")
)

col_fun_rest <- colorRamp2(
  c(-0.3, -0.09, 0, 0.09, 0.3), 
  c("#1A4475", "#6BAED6", "white", "#F4A582", "#B81D1D")
)

# --- CREATION OF THE 40MM ANNOTATION SPACER ---
spacer_annotation <- rowAnnotation(
  spacer = anno_empty(width = unit(40, "mm"), border = FALSE)
)


# --- 6. GENERATE HEATMAP 1: MASTER (ALL MOTIFS) ---
cat("[INFO] Rendering Heatmap 1 (Master - Sample and Motif Trees Active)...\n")
row_dend_master <- hclust(dist(matrix_master), method = "ward.D2")
row_dend_master <- as.dendrogram(row_dend_master)
row_dend_master <- dendextend::set(row_dend_master, "branches_lwd", dend_lwd)

col_dend_master <- hclust(dist(t(matrix_master)), method = "ward.D2")
col_dend_master <- as.dendrogram(col_dend_master)
col_dend_master <- dendextend::set(col_dend_master, "branches_lwd", dend_lwd)

png(filename = output_master, res = 300, width = 10500, height = 2600) # Increased height slightly for column dendrogram
p1 <- Heatmap(
  matrix_master,
  col = col_fun_master,
  cluster_rows = row_dend_master, 
  cluster_columns = col_dend_master,        
  show_row_dend = TRUE,          
  show_column_dend = TRUE,               # MODIFICATION: Enabled column dendrogram tracking tree for motifs
  row_dend_side = "left",
  column_dend_side = "top",
  row_dend_width = unit(55, "mm"),      
  column_dend_height = unit(40, "mm"),   # MODIFICATION: Set controlled vertical tree height spacing buffer
  left_annotation = spacer_annotation,  
  column_names_rot = 90,          
  row_names_side = "left",       
  row_names_gp = gpar(fontsize = 30, fontface = "plain"),     
  column_names_gp = gpar(fontsize = 18, fontface = "plain"),  
  show_heatmap_legend = FALSE      
)
draw(p1, padding = unit(c(2, 95, 2, 2), "mm")) 
dev.off()


# --- 7. GENERATE HEATMAP 2: PRE-SELECTED SUBSET (22 MOTIFS) ---
cat("[INFO] Rendering Heatmap 2 (22 Filtered Motifs - Sample and Motif Trees Active)...\n")
row_dend_subset <- hclust(dist(matrix_subset), method = "ward.D2")
row_dend_subset <- as.dendrogram(row_dend_subset)
row_dend_subset <- dendextend::set(row_dend_subset, "branches_lwd", dend_lwd)

col_dend_subset <- hclust(dist(t(matrix_subset)), method = "ward.D2")
col_dend_subset <- as.dendrogram(col_dend_subset)
col_dend_subset <- dendextend::set(col_dend_subset, "branches_lwd", dend_lwd)

png(filename = output_subset, res = 300, width = 5700, height = 2600) # Increased height slightly for column dendrogram
p2 <- Heatmap(
  matrix_subset,
  col = col_fun_subset,            
  cluster_rows = row_dend_subset,    
  cluster_columns = col_dend_subset, 
  column_dend_reorder = TRUE,
  show_row_dend = TRUE,          
  show_column_dend = TRUE,               # MODIFICATION: Enabled column dendrogram tracking tree for motifs
  row_dend_side = "left",
  column_dend_side = "top",
  row_dend_width = unit(55, "mm"),      
  column_dend_height = unit(40, "mm"),   # MODIFICATION: Set controlled vertical tree height spacing buffer
  left_annotation = spacer_annotation,  
  column_names_rot = 90,
  row_names_side = "left",       
  row_names_gp = gpar(fontsize = 30, fontface = "plain"),     
  column_names_gp = gpar(fontsize = 18, fontface = "plain"),  
  show_heatmap_legend = FALSE      
)
draw(p2, padding = unit(c(2, 95, 2, 2), "mm"))
dev.off()


# --- 8. GENERATE HEATMAP 3: REMAINING ELEMENTS (75 MOTIFS) ---
cat("[INFO] Rendering Heatmap 3 (75 Background Motifs - Sample and Motif Trees Active)...\n")
row_dend_rest <- hclust(dist(matrix_rest), method = "ward.D2")
row_dend_rest <- as.dendrogram(row_dend_rest)
row_dend_rest <- dendextend::set(row_dend_rest, "branches_lwd", dend_lwd)

col_dend_rest <- hclust(dist(t(matrix_rest)), method = "ward.D2")
col_dend_rest <- as.dendrogram(col_dend_rest)
col_dend_rest <- dendextend::set(col_dend_rest, "branches_lwd", dend_lwd)

png(filename = output_rest, res = 300, width = 10500, height = 2600) # Increased height slightly for column dendrogram
p3 <- Heatmap(
  matrix_rest,
  col = col_fun_rest,              
  cluster_rows = row_dend_rest,      
  cluster_columns = col_dend_rest, 
  column_dend_reorder = TRUE,
  show_row_dend = TRUE,          
  show_column_dend = TRUE,               # MODIFICATION: Enabled column dendrogram tracking tree for motifs
  row_dend_side = "left",
  column_dend_side = "top",
  row_dend_width = unit(55, "mm"),      
  column_dend_height = unit(40, "mm"),   # MODIFICATION: Set controlled vertical tree height spacing buffer
  left_annotation = spacer_annotation,  
  column_names_rot = 90,
  row_names_side = "left",       
  row_names_gp = gpar(fontsize = 30, fontface = "plain"),     
  column_names_gp = gpar(fontsize = 18, fontface = "plain"),  
  show_heatmap_legend = FALSE      
)
draw(p3, padding = unit(c(2, 95, 2, 2), "mm"))
dev.off()


# --- 9. EXPORT SEPARATE INDIVIDUAL COLOR MAP LEGEND ASSETS ---
cat("[INFO] Saving isolated, custom-scaled legends to individual files...\n")

# A. Standalone Legend for Heatmap 1 (All Motifs Range)
legend_1 <- Legend(
  col_fun = col_fun_master, 
  title = "Heatmap 1: C→T Fraction\n(Global Data Range)",
  title_gp = gpar(fontsize = 26, fontface = "plain", rot = 90),
  labels_gp = gpar(fontsize = 24),
  direction = "vertical",
  grid_width = unit(28, "mm"),       
  legend_height = unit(180, "mm"),   
  title_gap = unit(25, "mm"),        
  title_position = "topleft"
)
png(filename = legend_master_out, res = 300, width = 1800, height = 3400)
grid.newpage()
draw(legend_1, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "center")
dev.off()

# B. Standalone Legend for Heatmap 2 (22 Filtered Motifs Range: -0.4 to +0.4)
legend_2 <- Legend(
  col_fun = col_fun_subset, 
  title = "Heatmap 2: C→T Fraction\n(Filtered Scale: -0.4 to +0.4)",
  title_gp = gpar(fontsize = 26, fontface = "plain", rot = 90),
  labels_gp = gpar(fontsize = 24),
  direction = "vertical",
  grid_width = unit(28, "mm"),       
  legend_height = unit(180, "mm"),   
  title_gap = unit(25, "mm"),        
  title_position = "topleft"
)
png(filename = legend_subset_out, res = 300, width = 1800, height = 3400)
grid.newpage()
draw(legend_2, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "center")
dev.off()

# C. Standalone Legend for Heatmap 3 (75 Background Motifs Range: -0.3 to +0.3)
legend_3 <- Legend(
  col_fun = col_fun_rest, 
  title = "Heatmap 3: C→T Fraction\n(Background Scale: -0.3 to +0.3)",
  title_gp = gpar(fontsize = 26, fontface = "plain", rot = 90),
  labels_gp = gpar(fontsize = 24),
  direction = "vertical",
  grid_width = unit(28, "mm"),       
  legend_height = unit(180, "mm"),   
  title_gap = unit(25, "mm"),        
  title_position = "topleft"
)
png(filename = legend_rest_out, res = 300, width = 1800, height = 3400)
grid.newpage()
draw(legend_3, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "center")
dev.off()

cat(sprintf("\n🚀 SUCCESS! Motif column dendrograms are active. High-resolution images exported cleanly to:\n%s\n", promoter_dir))