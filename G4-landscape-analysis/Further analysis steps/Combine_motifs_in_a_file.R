# ========================================================================
# AUTOMATED LARGE-SCALE MOTIF GRID COMPILATION ENGINE
# ========================================================================
# Objective: Stitch high volumes of individual motif logo charts 
#            into high-resolution, aligned composite grids per category.
# ========================================================================

# --- 0. LOAD REQUIRED LIBRARIES ---
suppressPackageStartupMessages({
  library(magick)
  library(stringr)
  library(purrr)
})

# --- 1. SETUP FILE PATH STRUCTURES ---
base_input_dir  <- "Path/G4/G4_Final_Analysis_Complete"
output_grid_dir <- file.path(base_input_dir, "Combined_Master_Grids")

# Safely establish the destination grid repository directory
dir.create(output_grid_dir, showWarnings = FALSE, recursive = TRUE)

# Compartment routing vectors matching your local folder structures
contexts    <- c("Genome", "Promoter")
enrichments <- c("Homeotherms_highly_enriched", "Homeotherms_moderately_enriched",
                 "Homeotherms_highly_reduced", "Homeotherms_moderately_reduced")

# --- 2. HIGH-VOLUME GRID GRAPHICS ENGINE FUNCTION ---
compile_large_motif_folder_to_grid <- function(input_path, output_filepath, context_label, tier_label) {
  
  # Scan the current directory pass for targeted clustered grid plots
  img_files <- list.files(input_path, pattern = "_Clustered_Grid\\.png$", full.names = TRUE)
  
  if (length(img_files) == 0) {
    cat(sprintf("[SKIP] No clustered grid images found in: %s/%s\n", context_label, tier_label))
    return(NULL)
  }
  
  total_plots <- length(img_files)
  cat(sprintf("[PROCESSING] Stitched matrix pass for %s - %s (%d files found)...\n", 
              context_label, tier_label, total_plots))
  
  # Load the entire sequence of raw target image paths into an image stack vector
  loaded_images <- tryCatch({
    image_read(img_files)
  }, error = function(e) {
    cat(sprintf("[ERROR] Failed to read structural image objects in: %s\n", input_path))
    return(NULL)
  })
  
  if (is.null(loaded_images)) return(NULL)
  
  # --- OPTIMIZED HIGH-VOLUME COLUMN MATRIX ALLOCATION ---
  # Dynamically widens matrix rows based on the extensive file arrays (91-97 files total)
  if (total_plots <= 4) {
    grid_cols <- total_plots
  } else if (total_plots <= 10) {
    grid_cols <- 3
  } else if (total_plots <= 18) {
    grid_cols <- 4
  } else if (total_plots <= 28) {
    grid_cols <- 5
  } else {
    grid_cols <- 6 # 6 columns stretches higher loads (e.g. 30+ files) into a readable, symmetric layout
  }
  
  # Clean up formatting for the master category page banner title
  clean_tier_name <- str_to_title(gsub("_", " ", tier_label))
  header_text     <- sprintf("%s: %s Motif Alignments (Total: %d)", context_label, clean_tier_name, total_plots)
  
  # --- GRAPHICS ASSEMBLY LINE EXECUTION PASS ---
  cat(sprintf("  -> Arranging sub-panels into an optimized %d-column canvas...\n", grid_cols))
  
  # Build row-ordered image montage panel layout grids
  grid_canvas <- image_montage(
    loaded_images, 
    tile     = as.character(grid_cols), 
    geometry = "x1200+30+30", # Locks vertical alignment height to 1200px with a 30px clear safety margin gutter
    bg       = "white"
  )
  
  # Build a separate independent master text header block canvas matching the exact width profile
  header_canvas <- image_blank(width = image_info(grid_canvas)$width[1], height = 180, color = "white") %>%
    image_annotate(
      text    = header_text, 
      gravity = "center", 
      size    = 65, 
      color   = "black", 
      font    = "sans"
    )
  
  # FIXED: Concatenate the components as a combined image vector array pass to achieve vertical stacking
  final_composite <- image_append(c(header_canvas, grid_canvas), stack = TRUE)
  
  # Export the consolidated master chart to disk at complete lossless texture rendering specifications
  image_write(final_composite, path = output_filepath, format = "png", quality = 100)
  cat(sprintf("[SUCCESS] Master composite canvas written smoothly to:\n  %s\n\n", output_filepath))
  
  # Free system memory buffers
  rm(loaded_images, grid_canvas, header_canvas, final_composite)
  gc()
}

# --- 3. MAIN RUNTIME ITERATION CONTROLLER LOOP ---
cat("========================================================================\n")
cat("[INFO] Launching Large-Scale Automated Motif Grid Compilation Engine\n")
cat("========================================================================\n")

for (ctx in contexts) {
  for (enr in enrichments) {
    
    target_path <- file.path(base_input_dir, ctx, enr)
    
    if (!dir.exists(target_path)) {
      cat(sprintf("[WARN] Target subfolder path absent on system disk: %s\n", target_path))
      next
    }
    
    # Generate clean output grid file identifiers
    out_filename  <- sprintf("Master_Grid_%s_%s.png", ctx, enr)
    full_out_path <- file.path(output_grid_dir, out_filename)
    
    # Execute folder compilation cycle pass
    compile_large_motif_folder_to_grid(
      input_path      = target_path, 
      output_filepath = full_out_path, 
      context_label   = ctx, 
      tier_label      = enr
    )
  }
}

cat("========================================================================\n")
cat("*** 📊 PIPELINE COMPLETE: COMPOSITE MATRICES STITCHED SUCCESSFULLY ***\n")
cat(sprintf("Your publication-ready composite graphics maps are saved inside:\n%s\n", output_grid_dir))
cat("========================================================================\n")