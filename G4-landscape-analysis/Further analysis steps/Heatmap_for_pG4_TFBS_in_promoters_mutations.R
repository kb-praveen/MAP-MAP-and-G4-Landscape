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
output_file  <- file.path(promoter_dir, "CtoT_fraction_heatmap_final_tall_dend.png")

# --- ADJUSTABLE PARAMETERS ---
dend_lwd       <- 4     # Thickened dendrogram lines
sample_columns <- c("Lc", "Ac", "Hs", "Gg", "Ev", "Input")
# ----------------------------------

# --- 2. IMPORT DATA & CONSERVE MOTIF NAMES ---
if (!file.exists(input_file)) {
  stop(sprintf("[ERROR] Could not find required summary workbook at:\n%s", input_file))
}

# Read the specific "C_to_T" sheet directly from the summary spreadsheet workbook
raw_df <- read_excel(input_file, sheet = "C_to_T")
raw_df <- as.data.frame(raw_df)

# Assign the 'Motif Name' column as row names (falls back to ID column if name is blank)
rownames(raw_df) <- make.unique(if_else(is.na(raw_df$`Motif Name`) | raw_df$`Motif Name` == "", raw_df[[1]], raw_df$`Motif Name`))

# --- 3. CONVERT MATRIX NUMERIC VALUES ---
sub_df <- raw_df[, sample_columns, drop = FALSE]
sub_df[] <- lapply(sub_df, as.numeric)
matrix_data <- data.matrix(sub_df)

# MODIFICATION: Enforce a rigid color scale mapping bounded strictly between -1 and 1
col_fun <- colorRamp2(
  c(-0.5, 0, 0.5), 
  c("#2B61A1", "white", "#B81D1D") # Blue (Low/Negative), White (Zero), Red (High/Positive)
)

# --- 4. CALCULATE DENDROGRAM CLUSTERING STRUCTURES VIA WARD.D2 ---
row_dend <- hclust(dist(matrix_data), method = "ward.D2")
row_dend <- as.dendrogram(row_dend)
row_dend <- dendextend::set(row_dend, "branches_lwd", dend_lwd)

col_dend <- hclust(dist(t(matrix_data)), method = "ward.D2")
col_dend <- as.dendrogram(col_dend)
col_dend <- dendextend::set(col_dend, "branches_lwd", dend_lwd)

# --- 5. EXPORT MASTER HEATMAP (FIXED COLOR SCALE, MAXIMUM FONTS) ---
png(filename = output_file, res = 300, width = 3200, height = 8800)

p <- Heatmap(
  matrix_data,
  name = "C→T Fraction\n(Mean Subtracted)",
  col = col_fun,
  cluster_rows = row_dend,
  
  # Natural data-driven tree sorting enabled
  cluster_columns = col_dend,      
  column_dend_reorder = TRUE,      
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  column_dend_side = "top",        
  
  # Tall dendrogram line dimensions
  row_dend_width = unit(50, "mm"),
  column_dend_height = unit(30, "mm"),
  
  # Title configurations
  column_title = "C to T transition", 
  column_title_gp = gpar(fontsize = 30, fontface = "plain"), 
  
  # Sample labels configuration (Horizontal and Centered)
  column_names_rot = 0,
  column_names_centered = TRUE,  
  
  # Global font configurations
  row_names_gp = gpar(fontsize = 18, fontface = "plain"),     
  column_names_gp = gpar(fontsize = 28, fontface = "plain"),  
  
  # Layout metrics positioning the title text under the color bar
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 24, fontface = "plain"), 
    labels_gp = gpar(fontsize = 24),                     
    direction = "horizontal",
    legend_width = unit(180, "mm"),
    title_position = "topcenter"   # Centers title string text beneath horizontal color scales
  )
)

# Render plot to canvas with horizontal legend anchored below
draw(p, heatmap_legend_side = "bottom")
dev.off()

cat(sprintf("\n✅ SUCCESS! Fixed color scale [-1 to 1] applied. Heatmap updated smoothly at:\n%s\n", output_file))