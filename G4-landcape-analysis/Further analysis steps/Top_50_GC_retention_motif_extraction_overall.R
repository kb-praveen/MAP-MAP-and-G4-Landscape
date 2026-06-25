# 1. LOAD LIBRARIES
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

# --- CONFIGURATION ---
input_xlsx  <- "/Users/praveenbishnoi/Desktop/G4_submission/Supplemantary_tables_final/Table_S19_.xlsx"
output_dir  <- "/Users/praveenbishnoi/Desktop/G4_Final_Analysis_Complete/GC_bias_analysis_Final/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 2. LOAD DATA FROM MULTIPLE SHEETS
message("Step 1: Loading Raw Data (Sheet 2) and Motif Names (Sheet 4)...")
raw_data     <- read_excel(input_xlsx, sheet = 2)
motif_names  <- read_excel(input_xlsx, sheet = 4) %>% 
  select(Motif_ID, Motif_Name) %>% 
  distinct()

# 3. RESHAPE AND JOIN
message("Step 2: Reshaping and Joining Data...")
refined_data <- raw_data %>%
  separate(Category, into = c("Region", "Type"), sep = "_") %>%
  mutate(Group = paste(Thermal_Group, Type, sep = "_")) %>%
  select(Motif_ID, Enrichment_Tier, Region, Group, Mean_GC, Hit_Count) %>%
  pivot_wider(
    id_cols = c(Motif_ID, Enrichment_Tier, Region),
    names_from = Group,
    values_from = c(Mean_GC, Hit_Count)
  ) %>%
  left_join(motif_names, by = "Motif_ID")

# 4. STATISTICAL CALCULATIONS (Chi-Square)
message("Step 3: Calculating Statistics...")
final_analysis <- refined_data %>%
  mutate(
    Statistical_Test = "Chi-Square Test of Proportions",
    Absolute_G4_Shift = Mean_GC_Homeotherm_G4 - Mean_GC_Poikilotherm_G4,
    Homeo_Delta   = Mean_GC_Homeotherm_G4 - Mean_GC_Homeotherm_Free,
    Poikilo_Delta = Mean_GC_Poikilotherm_G4 - Mean_GC_Poikilotherm_Free,
    Delta_Delta_Shift = Homeo_Delta - Poikilo_Delta,
    Dominant_Clade = ifelse(Absolute_G4_Shift > 0, "Homeotherm", "Poikilotherm"),
    
    # Chi-Square Math (GC vs AT balance)
    H_N = Hit_Count_Homeotherm_G4,
    P_N = Hit_Count_Poikilotherm_G4,
    H_GC_Count = H_N * (Mean_GC_Homeotherm_G4 / 100),
    P_GC_Count = P_N * (Mean_GC_Poikilotherm_G4 / 100),
    Total_N = H_N + P_N,
    Total_GC = H_GC_Count + P_GC_Count,
    Prob_GC = Total_GC / Total_N,
    Exp_H_GC = H_N * Prob_GC,
    Exp_P_GC = P_N * Prob_GC,
    
    Chisq_Stat = ((H_GC_Count - Exp_H_GC)^2 / Exp_H_GC) + 
      ((P_GC_Count - Exp_P_GC)^2 / Exp_P_GC) +
      (((H_N - H_GC_Count) - (H_N - Exp_H_GC))^2 / (H_N - Exp_H_GC)) +
      (((P_N - P_GC_Count) - (P_N - Exp_P_GC))^2 / (P_N - Exp_P_GC)),
    
    P_Value = pchisq(Chisq_Stat, df = 1, lower.tail = FALSE)
  ) %>%
  mutate(P_Adj = p.adjust(P_Value, method = "BH")) %>%
  filter(P_Adj < 0.05)

# 5. EXPORT TOP 50 TSV FILES
message("Step 4: Exporting TSV Results...")

export_tsv <- function(df, metric_col, method_name) {
  # Overall
  top_ov <- df %>% arrange(desc(!!sym(metric_col)), P_Adj) %>% head(50) %>%
    select(Motif_ID, Motif_Name, Enrichment_Tier, Region, !!sym(metric_col), 
           Dominant_Clade, P_Adj, Statistical_Test, everything())
  write_tsv(top_ov, paste0(output_dir, "Top50_", method_name, "_Overall.tsv"))
  
  # By Region
  top_reg <- df %>% group_by(Region) %>% 
    slice_max(order_by = tibble(!!sym(metric_col), -P_Adj), n = 50) %>% ungroup() %>%
    select(Motif_ID, Motif_Name, Enrichment_Tier, Region, !!sym(metric_col), 
           Dominant_Clade, P_Adj, Statistical_Test, everything())
  write_tsv(top_reg, paste0(output_dir, "Top50_", method_name, "_By_Region.tsv"))
}

export_tsv(final_analysis, "Absolute_G4_Shift", "Absolute_Retention")
export_tsv(final_analysis, "Delta_Delta_Shift", "Background_Corrected_Retention")

# Master Stats TSV
write_tsv(final_analysis, paste0(output_dir, "Master_Retention_Analysis_Full.tsv"))

message("DONE! Files saved as TSV with Motif Names and Test metadata.")