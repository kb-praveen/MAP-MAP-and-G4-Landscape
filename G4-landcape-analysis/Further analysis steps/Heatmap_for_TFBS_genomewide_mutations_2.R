# --- Load Required Libraries ---
library(readxl)         # For reading Excel files (.xlsx)
library(ComplexHeatmap) # For advanced heatmap plotting
library(circlize)       # For color mapping (used by ComplexHeatmap)
library(dendextend)     # For customizing dendrograms (thick branches, etc.)
library(grid)           # For graphical parameters (e.g., units, gpar)

# --- USER INPUT: File Paths ---
input_file <- "/Users/praveenbishnoi/Desktop/GC_rich_TFBS_mutation_profile_genome_wide/final_CtoT_fraction_subtracted.xlsx"  # Path to input data file
output_file_with_col_dend <- "/Users/praveenbishnoi/Desktop/GC_rich_TFBS_mutation_profile_genome_wide/CtoT_fraction_heatmap_with_col_dend.png"  # Output file: heatmap with both dendrograms
output_file_without_col_dend <- "/Users/praveenbishnoi/Desktop/GC_rich_TFBS_mutation_profile_genome_wide/CtoT_fraction_heatmap_without_col_dend.png"  # Output file: heatmap with only row dendrogram

# --- USER ADJUSTABLE PARAMETERS ---
column_width_mm <- 6      # Width of each column tile (cell) in millimeters
row_height_mm <- 1.5      # Height of each row tile (cell) in millimeters
column_title_angle <- 0   # Angle for the main column title (0 = horizontal, 90 = vertical)
row_title_angle <- 0      # Angle for the main row title (0 = horizontal, 90 = vertical)
dend_lwd <- 2             # Dendrogram line width (thickness)
# ----------------------------------

# --- Data Import and Preparation ---
df <- readxl::read_excel(input_file) # Read Excel file into a data frame
df <- as.data.frame(df)              # Ensure it's a data frame
rownames(df) <- df[[1]]              # Set first column as row names (motif names)
df <- df[ , -1]                      # Remove first column (now redundant)
df[] <- lapply(df, as.numeric)       # Convert all columns to numeric

# --- Color Scale Definition ---
col_fun <- colorRamp2(
  c(min(df, na.rm = TRUE), 0, max(df, na.rm = TRUE)),  # Range: min to 0 to max
  c("blue", "white", "red")                            # Colors: blue (low), white (zero), red (high)
)

# --- Prepare Dendrograms with Custom Thickness ---
row_dend <- hclust(dist(df), method = "ward.D2")       # Hierarchical clustering for rows
row_dend <- as.dendrogram(row_dend)                    # Convert to dendrogram object
row_dend <- dendextend::set(row_dend, "branches_lwd", dend_lwd) # Set row dendrogram line width

col_dend <- hclust(dist(t(df)), method = "ward.D2")    # Hierarchical clustering for columns
col_dend <- as.dendrogram(col_dend)                    # Convert to dendrogram object
col_dend <- dendextend::set(col_dend, "branches_lwd", dend_lwd) # Set column dendrogram line width

# --- 1. Heatmap WITH Column Dendrogram ---
png(filename = output_file_with_col_dend, width = 2000, height = 7000, res = 300)  # Open PNG device (image file)
Heatmap(
  as.matrix(df),                       # Data matrix for heatmap
  name = "C→T",                        # Legend title
  col = col_fun,                       # Color mapping function
  cluster_rows = row_dend,             # Use custom row dendrogram
  cluster_columns = col_dend,          # Use custom column dendrogram
  show_row_dend = TRUE,                # Show row dendrogram
  show_column_dend = TRUE,             # Show column dendrogram
  column_title = "C→T Transition in GC rich TFBS", # Main title above columns
  column_title_gp = gpar(fontsize = 16, col = "black"), # Main column title style
  column_title_rot = column_title_angle,              # Main column title angle
  column_names_rot = 90,                # Rotate column labels (sample names) vertically
  row_title = NULL,                     # No main row title
  row_title_gp = gpar(fontsize = 8, col = "black"),   # (If row title used) style
  row_title_rot = row_title_angle,                     # (If row title used) angle
  row_names_gp = gpar(fontsize = 4),                  # Row label (motif) font size
  column_names_gp = gpar(fontsize = 8),               # Column label (sample) font size
  width = unit(ncol(df) * column_width_mm, "mm"),     # Total heatmap width (columns × tile width)
  height = unit(nrow(df) * row_height_mm, "mm")       # Total heatmap height (rows × tile height)
  # No heatmap_legend_param: legend uses default settings
)
dev.off()  # Close PNG device

# --- 2. Heatmap WITHOUT Column Dendrogram ---
png(filename = output_file_without_col_dend, width = 2000, height = 7000, res = 300)  # Open PNG device
Heatmap(
  as.matrix(df),
  name = "C→T",
  col = col_fun,
  cluster_rows = row_dend,
  cluster_columns = FALSE,              # Do not cluster columns (keep input order)
  show_row_dend = TRUE,
  show_column_dend = FALSE,             # Do not show column dendrogram
  column_title = "C→T Transition in GC rich TFBS",
  column_title_gp = gpar(fontsize = 16, col = "black"),
  column_title_rot = column_title_angle,
  column_names_rot = 90,                # Rotate column labels vertically
  row_title = NULL,
  row_title_gp = gpar(fontsize = 8, col = "black"),
  row_title_rot = row_title_angle,
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 8),
  width = unit(ncol(df) * column_width_mm, "mm"),
  height = unit(nrow(df) * row_height_mm, "mm")
  # No heatmap_legend_param: legend uses default settings
)
dev.off()  # Close PNG device
