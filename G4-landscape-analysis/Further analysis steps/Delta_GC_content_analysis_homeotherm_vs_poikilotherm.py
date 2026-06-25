import os
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import shutil
from multiprocessing import Pool

# --- PATH CONFIGURATION ---
BASE = "Path/G4_motif_analysis"
INPUT_DIR = f"{BASE}/Results/Homeotherm_vs_Poikilotherm"
NEW_DIR = f"{INPUT_DIR}/GC_content_analysis"

# Create new directory
os.makedirs(NEW_DIR, exist_ok=True)

def run_analysis():
    # 1. Copy the source data to the new directory to preserve original
    source_file = f"{INPUT_DIR}/Detailed_Thermal_GC_Stats.tsv"
    local_copy = f"{NEW_DIR}/Detailed_Thermal_GC_Stats_Copy.tsv"
    shutil.copy2(source_file, local_copy)
    
    # 2. Load and Pivot Data
    df = pd.read_csv(local_copy, sep='\t')
    
    # Pivot to get Categories as columns so we can calculate Delta
    # We use pivot_table to handle any potential duplicate entries safely
    df_pivot = df.pivot_table(
        index=['Motif_ID', 'Enrichment_Tier', 'Thermal_Group', 'Clade', 'Species'],
        columns='Category',
        values='Mean_GC'
    ).reset_index()
    
    # 3. Calculate Deltas
    # Delta = (G4 Region GC) - (Background GC)
    df_pivot['Delta_Promoter'] = df_pivot['Promoter_G4'] - df_pivot['Promoter_Free']
    df_pivot['Delta_Genome'] = df_pivot['Genome_G4'] - df_pivot['Genome_Free']
    
    # Save the expanded Delta table
    df_pivot.to_csv(f"{NEW_DIR}/Motif_Delta_GC_Per_Species.tsv", sep='\t', index=False)
    
    # 4. Statistical Testing (Mann-Whitney U / Wilcoxon Rank Sum)
    stats_results = []
    tiers = [
        "Highly_Enriched_Homeotherm", "Moderately_Enriched_Homeotherm",
        "Highly_Reduced_Homeotherm", "Moderately_Reduced_Homeotherm"
    ]
    regions = ['Delta_Promoter', 'Delta_Genome']
    
    for region in regions:
        for tier in tiers:
            subset = df_pivot[df_pivot['Enrichment_Tier'] == tier].dropna(subset=[region])
            
            homeo = subset[subset['Thermal_Group'] == 'Homeotherm'][region]
            poikilo = subset[subset['Thermal_Group'] == 'Poikilotherm'][region]
            
            if len(homeo) > 0 and len(poikilo) > 0:
                stat, p = mannwhitneyu(homeo, poikilo, alternative='two-sided')
                
                stats_results.append({
                    'Region': region.replace('Delta_', ''),
                    'Enrichment_Tier': tier,
                    'P_Value': p,
                    'Homeo_Mean_Delta': homeo.mean(),
                    'Poikilo_Mean_Delta': poikilo.mean(),
                    'Homeo_N': len(homeo),
                    'Poikilo_N': len(poikilo)
                })

    stats_df = pd.DataFrame(stats_results)
    stats_df.to_csv(f"{NEW_DIR}/Delta_GC_Statistical_Significance.tsv", sep='\t', index=False)
    
    print(f"Analysis complete. All files saved in: {NEW_DIR}")

if __name__ == "__main__":
    # While the data processing is mostly pandas-based, 
    # the workflow is structured to respect your 48-core environment for future scaling.
    run_analysis()