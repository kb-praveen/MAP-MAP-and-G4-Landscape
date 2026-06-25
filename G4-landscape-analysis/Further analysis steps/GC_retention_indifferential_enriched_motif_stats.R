# 1. LOAD LIBRARIES
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

# --- SETUP PATHS ---
raw_species_file <- "Path/Motif_GC_content_summary.xlsx"
name_map_file   <- "Path/Checking/Motif_ID_and_name.xlsx"
output_dir      <- "Path/Checking/GC_retention_analysis"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 2. LOAD DATA & CLEAN STRUCTURE ---
if (!file.exists(name_map_file) | !file.exists(raw_species_file)) {
  stop("[ERROR] One or more input files are missing from your Desktop directory.")
}

names_df   <- read_excel(name_map_file)
id_to_name <- deframe(names_df %>% dplyr::select(Motif_ID, Motif_Name))
raw_data   <- read_excel(raw_species_file)

# Standardize names and fix sorting factors
gc_clean <- raw_data %>%
  mutate(
    Motif_Name = id_to_name[Motif_ID],
    Enrichment_Tier = gsub("_", " ", Enrichment_Tier),
    Enrichment_Tier = factor(Enrichment_Tier, levels = c("Highly Enriched", "Moderately Enriched", 
                                                         "Moderately Reduced", "Highly Reduced")),
    Thermal_Group = factor(Thermal_Group, levels = c("Poikilotherm", "Homeotherm"))
  ) %>%
  # Focus exclusively on significant enrichment/reduction tiers (excludes "Not Significant")
  filter(Enrichment_Tier %in% c("Highly Enriched", "Moderately Enriched", "Moderately Reduced", "Highly Reduced"))

# --- 3. ALL-MOTIF UNPAIRED SPECIES PROCESSING FUNCTION ---
process_all_motifs_species_wise <- function(df, context_keyword, filename_base) {
  
  # Filter context layer (Genome or Promoter)
  context_df <- df %>% 
    filter(str_detect(Category, regex(context_keyword, ignore_case = TRUE))) %>%
    mutate(Region_Type = if_else(str_detect(Category, "_G4"), "G4-forming regions", "G4-free regions"),
           Region_Type = factor(Region_Type, levels = c("G4-forming regions", "G4-free regions")))
  
  if (nrow(context_df) == 0) return(NULL)
  
  # --- STEP A: RESHAPE SPECIES REPLICATES INTO A WIDE FORMAT PER MOTIF ---
  # This isolates species values to calculate clean descriptive parameters
  motif_wide_summary <- context_df %>%
    group_by(Motif_ID, Motif_Name, Enrichment_Tier, Region_Type) %>%
    summarise(
      Mean_GC_Homeotherm    = mean(Mean_GC[Thermal_Group == "Homeotherm"], na.rm = TRUE),
      Mean_GC_Poikilotherm  = mean(Mean_GC[Thermal_Group == "Poikilotherm"], na.rm = TRUE),
      Hit_Count_Homeotherm  = sum(Hit_Count[Thermal_Group == "Homeotherm"], na.rm = TRUE),
      Hit_Count_Poikilotherm = sum(Hit_Count[Thermal_Group == "Poikilotherm"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Individual_GC_Delta = Mean_GC_Homeotherm - Mean_GC_Poikilotherm
    )
  
  # --- STEP B: RUN INDEPENDENT UNPAIRED TESTS PER MOTIF ---
  # We loop through each row to collect the actual unpaired statistical values across species vectors
  test_results <- context_df %>%
    group_by(Motif_ID, Motif_Name, Enrichment_Tier, Region_Type) %>%
    do({
      homeo_vector <- .$Mean_GC[.$Thermal_Group == "Homeotherm"]
      poiki_vector <- .$Mean_GC[.$Thermal_Group == "Poikilotherm"]
      
      # Safety validation check to guarantee statistical variance exists in both clades
      if (length(na.omit(homeo_vector)) >= 2 & length(na.omit(poiki_vector)) >= 2) {
        u_test <- wilcox.test(homeo_vector, poiki_vector, paired = FALSE)
        p_val  <- u_test$p.value
      } else {
        p_val  <- NA_real_
      }
      tibble(Motif_Individual_Wilcox_P = p_val)
    }) %>%
    ungroup()
  
  # --- STEP C: INTEGRATE METRICS & CALCULATE BH ADJUSTMENT ---
  final_matrix <- motif_wide_summary %>%
    left_join(test_results, by = c("Motif_ID", "Motif_Name", "Enrichment_Tier", "Region_Type")) %>%
    group_by(Region_Type) %>%
    # Multi-test correction applied across the total test pool within each architectural panel
    mutate(Motif_Individual_Wilcox_Padj_BH = p.adjust(Motif_Individual_Wilcox_P, method = "BH")) %>%
    ungroup() %>%
    arrange(Enrichment_Tier, desc(Individual_GC_Delta)) %>%
    # Restructure into your exact clean column requirements
    dplyr::select(
      Motif_ID, 
      Motif_Name, 
      Enrichment_Tier, 
      Region_Type, 
      Mean_GC_Homeotherm, 
      Mean_GC_Poikilotherm, 
      Hit_Count_Homeotherm, 
      Hit_Count_Poikilotherm, 
      Individual_GC_Delta, 
      Motif_Individual_Wilcox_P, 
      Motif_Individual_Wilcox_Padj_BH
    )
  
  # Output file generation
  write_tsv(final_matrix, file.path(output_dir, paste0("Species_Wise_Individual_Motifs_Stats_", filename_base, ".tsv")))
}

# --- 4. RUN COMPREHENSIVE BATCH ANALYSIS ---
process_all_motifs_species_wise(gc_clean, "Genome", "Genomewide")
process_all_motifs_species_wise(gc_clean, "Promoter", "Promoters")

cat(sprintf("\n[SUCCESS] Independent species-level pipeline complete! All individual motif p-values successfully output into:\n%s\n", output_dir))