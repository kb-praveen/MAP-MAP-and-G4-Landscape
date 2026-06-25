# ========================================================================
# HOMEO-MINUS-POIKILO SUBTRACTION GRAPHICS GRID COMPILATION ENGINE
# ========================================================================
# Objective: Process flat 'Genome' and 'Promoter' folders, dynamically
#            group files by enrichment tier from names, and output grids.
# ========================================================================

# --- 0. LOAD REQUIRED LIBRARIES ---
suppressPackageStartupMessages({
  library(magick)
  library(stringr)
  library(dplyr)
  library(purrr)
})

# --- 1. SETUP FILE PATH STRUCTURES ---
base_input_dir  <- "Path/G4/G4_Final_Analysis_Complete/Subtracted_PWMs/Homeo_minus_poikilo"
output_grid_dir <- file.path(base_input_dir, "Combined_Subtraction_Grids")

# Safely establish the destination grid repository directory
dir.create(output_grid_dir, showWarnings = FALSE, recursive = TRUE)

# Main context compartments
contexts <- c("Genome", "Promoter")

# --- 2. GRID ASSEMBLY GENERATION WORKER PASS ---
assemble_tier_grid <- function(img_paths, context_label, tier_label, output_filepath) {
  total_plots <- length(img_paths)
  
  # Load the raw target image paths into an image stack vector array
  loaded_images <- tryCatch({
    image_read(img_paths)
  }, error = function(e) {
    cat(sprintf("  [ERROR] Failed to read image vectors for tier: %s\n", tier_label))
    return(NULL)
  })
  
  if (is.null(loaded_images)) return(NULL)
  
  # --- OPTIMIZED COLUMN MATRIX ALLOCATION MATRIX ---
  if (total_plots <= 4) {
    grid_cols <- total_plots
  } else if (total_plots <= 10) {
    grid_cols <- 3
  } else if (total_plots <= 18) {
    grid_cols <- 4
  } else if (total_plots <= 28) {
    grid_cols <- 5
  } else {
    grid_cols <- 6
  }
  
  # Format clean master string labels for the top page layout banner 
  clean_tier_name <- str_to_title(gsub("_", " ", tier_label))
  header_text     <- sprintf("%s: Δ Frequency (%s) Matrix [Total: %d]", 
                             context_label, clean_tier_name, total_plots)
  
  cat(sprintf("  -> Arranging %d plots into a clean %d-column collage grid...\n", total_plots, grid_cols))
  
  # Build row-ordered image montage canvas panels
  grid_canvas <- image_montage(
    loaded_images, 
    tile     = as.character(grid_cols), 
    geometry = "x1200+30+30", # Locks vertical height alignment safely to 1200px with a 30px clear margin
    bg       = "white"
  )
  
  # Build a matching text title element header block canvas
  header_canvas <- image_blank(width = image_info(grid_canvas)$width[1], height = 180, color = "white") %>%
    image_annotate(
      text    = header_text, 
      gravity = "center", 
      size    = 60, 
      color   = "black", 
      font    = "sans"
    )
  
  # Concatenate the title layer and panel vector grid seamlessly together
  final_composite <- image_append(c(header_canvas, grid_canvas), stack = TRUE)
  
  # Export the consolidated graphic map out to disk
  image_write(final_composite, path = output_filepath, format = "png", quality = 100)
  cat(sprintf("  [SUCCESS] Grid output written smoothly to:\n    %s\n\n", output_filepath))
  
  # Free system memory buffers
  rm(loaded_images, grid_canvas, header_canvas, final_composite)
  gc()
}

# --- 3. MAIN RUNTIME ITERATION CONTROLLER LOOP ---
cat("========================================================================\n")
cat("[INFO] Launching Flat-Directory Subtraction Grid Assembly Engine\n")
cat("========================================================================\n")

for (ctx in contexts) {
  target_path <- file.path(base_input_dir, ctx)
  
  if (!dir.exists(target_path)) {
    cat(sprintf("[WARN] Context folder absent on system disk: %s\n", target_path))
    next
  }
  
  cat(sprintf("\n[PROCESSING CONTEXT] Scanning Flat Directory Pass: %s\n", ctx))
  
  # Find all Homeo-minus-Poikilo subtraction image plots inside the flat context folder
  all_files <- list.files(target_path, pattern = "_Homeo_minus_Poikilo\\.png$", full.names = TRUE)
  
  if (length(all_files) == 0) {
    cat(sprintf("  [SKIP] No compatible subtraction plots found inside: %s\n", ctx))
    next
  }
  
  # --- DYNAMIC FILENAME PARSING & TIER GROUPING SYSTEM ---
  # Instantly map and isolate string metadata patterns out from file signatures
  file_metadata_df <- tibble(file_path = all_files) %>%
    mutate(
      file_name = basename(file_path),
      # Drops the motif prefix ID and trailing suffix string to extract exact raw category tiers
      tier_extracted = str_match(file_name, "^[A-Za-z0-9\\.]+(?:_Part[0-9]+)?_(Highly_Enriched|Moderately_Enriched|Highly_Reduced|Moderately_Reduced)_")[,2]
    ) %>%
    filter(!is.na(tier_extracted)) # Drops un-parseable file artifacts safely
  
  if (nrow(file_metadata_df) == 0) {
    cat("  [WARN] Filenames could not be categorized. Check alignment naming formats.\n")
    next
  }
  
  # Split paths based on our newly extracted enrichment tier groups
  tier_groups <- split(file_metadata_df$file_path, file_metadata_df$tier_extracted)
  
  for (tier_name in names(tier_groups)) {
    paths_to_process <- tier_groups[[tier_name]]
    
    # Generate clean output naming conventions
    out_filename  <- sprintf("Master_Grid_Subtraction_%s_%s.png", ctx, tier_name)
    full_out_path <- file.path(output_grid_dir, out_filename)
    
    # Run the compiled graphics generation process for the isolated subset array
    assemble_tier_grid(
      img_paths       = paths_to_process, 
      context_label   = ctx, 
      tier_label      = tier_name, 
      output_filepath = full_out_path
    )
  }
}

cat("========================================================================\n")
cat("*** 📊 PIPELINE COMPLETE: FLAT COMPARTMENT SUBTRACTIONS MATRIX COMPILED ***\n")
cat(sprintf("Your clean subtraction composite grid maps are located inside:\n%s\n", output_grid_dir))
cat("========================================================================\n")