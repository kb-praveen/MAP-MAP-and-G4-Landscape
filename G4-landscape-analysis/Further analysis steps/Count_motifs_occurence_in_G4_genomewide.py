import csv
import gzip
import os
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool

# --- CONFIGURATION ---
CORES = 48
ID_NAME_FILE = Path('Path/G4_motif_analysis/Background_G4_free/JASPAR_vertebrate_Non_redundant_2026_Ids_names.txt')
G4_OUTPUT_BASE = Path('Path/G4_motif_analysis/Fimo_in_G4s_genomewide')
CLADE_ORDER = ['Pisces', 'Amphibians', 'Reptiles', 'Aves', 'Mammals']

def count_motifs_in_species(args):
    """Worker function to count motifs in a single species fimo file."""
    clade, s_name, fimo_path = args
    counts = defaultdict(int)
    try:
        with gzip.open(fimo_path, 'rt') as f:
            # Skip comments and use DictReader
            lines = (line for line in f if not line.startswith('#'))
            reader = csv.DictReader(lines, delimiter='\t')
            for row in reader:
                counts[row['motif_id']] += 1
        return s_name, counts
    except Exception as e:
        print(f"Error processing {s_name}: {e}")
        return s_name, {}

def main():
    print("="*60 + "\nPARALLEL AGGREGATION (48 CORES): G4 MOTIF COUNTS\n" + "="*60)
    
    # 1. Load Motif Metadata
    motif_ids, motif_info = [], {}
    with open(ID_NAME_FILE, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            mid = row['Motif_ID']
            motif_ids.append(mid)
            motif_info[mid] = row['Motif_Name']

    # 2. Prepare Tasks
    tasks = []
    for clade in CLADE_ORDER:
        clade_dir = G4_OUTPUT_BASE / clade
        if not clade_dir.exists(): continue
        
        species_dirs = sorted([d for d in clade_dir.iterdir() if d.is_dir()])
        for s_dir in species_dirs:
            target_fimo = s_dir / "fimo.tsv.gz"
            if target_fimo.exists() and target_fimo.stat().st_size > 100:
                tasks.append((clade, s_dir.name, target_fimo))

    # 3. Execute in Parallel
    print(f"[INFO] Processing {len(tasks)} species across {CORES} cores...")
    with Pool(processes=CORES) as pool:
        results = pool.map(count_motifs_in_species, tasks)

    # 4. Organize Results
    species_data = {s_name: counts for s_name, counts in results}
    # Maintain phylogenetic order for the header
    ordered_species_list = [t[1] for t in tasks]

    # 5. Write Final TSV
    output_file = G4_OUTPUT_BASE / "Motif_counts_in_G4_regions_genomewide_in_vertebrates.tsv"
    with open(output_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Motif_ID', 'Motif_Name'] + ordered_species_list)
        for mid in motif_ids:
            row = [mid, motif_info[mid]]
            for s_name in ordered_species_list:
                row.append(species_data[s_name].get(mid, 0))
            writer.writerow(row)
            
    print(f"\n[SUCCESS] G4 Count table generated at: {output_file}")

if __name__ == "__main__":
    main()