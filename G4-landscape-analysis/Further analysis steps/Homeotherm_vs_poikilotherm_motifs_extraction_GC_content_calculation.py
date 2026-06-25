import os
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from multiprocessing import Pool
import gzip
from collections import defaultdict

# --- PATH CONFIGURATION ---
BASE = "Path/G4_motif_analysis"
RESULTS_DIR = f"{BASE}/Results/Homeotherm_vs_Poikilotherm"
STATS_PROM = f"{BASE}/Results/Clade_Comparisons/Stats_Promoter_Poikilotherms_vs_Homeotherms.tsv"
STATS_GENO = f"{BASE}/Results/Clade_Comparisons/Stats_Genome_Poikilotherms_vs_Homeotherms.tsv"
OUT_FASTA = f"{RESULTS_DIR}/Motif_fasta"

os.makedirs(OUT_FASTA, exist_ok=True)
CORES = 48

# --- 1. LOAD & CLASSIFY MOTIFS BY REGION ---
def get_tiers(file_path):
    df = pd.read_csv(file_path, sep='\t')
    # Filter out NAs in Log2FC or Padj
    df = df.dropna(subset=['Log2FC', 'Padj'])
    
    def assign_tier(row):
        lfc = row['Log2FC']
        padj = row['Padj']
        if padj < 0.05 and lfc > 1: return "Highly_Enriched_Homeotherm"
        if padj < 0.05 and lfc > 0: return "Moderately_Enriched_Homeotherm"
        if padj < 0.05 and lfc < -1: return "Highly_Reduced_Homeotherm"
        if padj < 0.05 and lfc < 0: return "Moderately_Reduced_Homeotherm"
        return "Stable"
    
    df['Tier'] = df.apply(assign_tier, axis=1)
    return dict(zip(df['Motif_ID'].astype(str), df['Tier']))

# region_tiers holds two dictionaries: one for promoter stats, one for genome stats
region_tiers = {
    "Promoter": get_tiers(STATS_PROM),
    "Genome": get_tiers(STATS_GENO)
}

CAT_MAP = {
    "Promoter_G4": ("G4_in_promoters", "Promoter"),
    "Promoter_Free": ("Fimo_on_G4_free_promoters", "Promoter"),
    "Genome_G4": ("Fimo_in_G4s_genomewide", "Genome"),
    "Genome_Free": ("Fimo_G4_free_background", "Genome")
}

def process_task(args):
    cat_label, dir_name, stat_type, clade, species, thermal_group = args
    fimo_file = f"{BASE}/{dir_name}/{clade}/{species}/fimo.tsv.gz"
    
    if not os.path.exists(fimo_file):
        return []

    results = []
    motif_groups = defaultdict(list)
    tier_map = region_tiers[stat_type]

    try:
        with gzip.open(fimo_file, 'rt') as f:
            next(f) # Skip header
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) < 10: continue
                m_id = cols[0].strip()
                
                if m_id in tier_map:
                    tier = tier_map[m_id]
                    if tier == "Stable": continue
                    
                    out_path = f"{OUT_FASTA}/{cat_label}/{thermal_group}/{tier}/{species}"
                    os.makedirs(out_path, exist_ok=True)
                    
                    record = SeqRecord(Seq(cols[9]), id=f"{cols[2]}:{cols[3]}", description=m_id)
                    motif_groups[(m_id, out_path, tier)].append(record)

        for (m_id, out_path, tier), records in motif_groups.items():
            m_file = f"{out_path}/{m_id}.fasta.gz"
            with gzip.open(m_file, "wt") as f_out:
                SeqIO.write(records, f_out, "fasta")
            
            gc_vals = [gc_fraction(r.seq) * 100 for r in records]
            results.append({
                "Motif_ID": m_id,
                "Enrichment_Tier": tier,
                "Category": cat_label,
                "Thermal_Group": thermal_group,
                "Clade": clade,
                "Species": species,
                "Mean_GC": sum(gc_vals) / len(gc_vals),
                "Hit_Count": len(records)
            })
        return results
    except Exception:
        return []

if __name__ == "__main__":
    metadata = pd.read_csv(f"{BASE}/Table_S1.tsv", sep='\t')
    homeotherm_clades = ["Mammals", "Birds", "Aves"] 

    tasks = []
    for label, (d_name, s_type) in CAT_MAP.items():
        for _, row in metadata.iterrows():
            clade = str(row['Classification_1'])
            thermal = "Homeotherm" if any(h in clade for h in homeotherm_clades) else "Poikilotherm"
            spp = str(row['Species']).replace(" ", "_")
            tasks.append((label, d_name, s_type, clade, spp, thermal))

    print(f"Submitting {len(tasks)} tasks to {CORES} cores...")
    with Pool(CORES) as pool:
        all_results = pool.map(process_task, tasks)

    flat_data = [item for sublist in all_results if sublist for item in sublist]
    
    if flat_data:
        df = pd.DataFrame(flat_data)
        # Detailed output
        df.to_csv(f"{RESULTS_DIR}/Detailed_Thermal_GC_Stats.tsv", sep='\t', index=False)
        
        # Summary Perspective (Mean of GC per Motif per Thermal Group)
        summary = df.groupby(['Motif_ID', 'Enrichment_Tier', 'Category', 'Thermal_Group']).agg({
            'Mean_GC': 'mean',
            'Hit_Count': 'sum'
        }).reset_index()
        summary.to_csv(f"{RESULTS_DIR}/Summary_Thermal_GC_Averages.tsv", sep='\t', index=False)
        
        print(f"Success! Data saved to {RESULTS_DIR}")
    else:
        print("No motifs processed. Verify Motif IDs in Stats files.")