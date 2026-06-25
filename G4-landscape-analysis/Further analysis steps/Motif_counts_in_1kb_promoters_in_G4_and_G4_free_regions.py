import csv
import gzip
import os
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool

# --- GLOBAL CONFIGURATION ---
CORES = 48
ID_NAME_FILE = Path('Path/G4_motif_analysis/Background_G4_free/JASPAR_vertebrate_Non_redundant_2026_Ids_names.txt')
CLADE_ORDER = ['Pisces', 'Amphibians', 'Reptiles', 'Aves', 'Mammals']

# Define the two tasks
JOBS = [
    {
        "input_dir": Path('Path/G4_motif_analysis/Fimo_on_G4_free_promoters'),
        "output_name": "Motif_counts_in_G4_free_promoters.tsv"
    },
    {
        "input_dir": Path('Path/G4_motif_analysis/G4_in_promoters'),
        "output_name": "Motif_counts_in_G4_containing_promoters.tsv"
    }
]

def count_motifs_in_species(args):
    """Worker function to count motifs in a single species fimo file."""
    s_name, fimo_path = args
    counts = defaultdict(int)
    try:
        # Check if file exists and is not empty
        if not fimo_path.exists() or fimo_path.stat().st_size < 100:
            return s_name, {}
            
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

def run_aggregation(input_base, output_filename, motif_ids, motif_info):
    print(f"\n[STARTING] Processing directory: {input_base}")
    
    # 1. Prepare Tasks
    tasks = []
    for clade in CLADE_ORDER:
        clade_dir = input_base / clade
        if not clade_dir.exists(): continue
        
        # Collect species directories
        species_dirs = sorted([d for d in clade_dir.iterdir() if d.is_dir()])
        for s_dir in species_dirs:
            target_fimo = s_dir / "fimo.tsv.gz"
            if target_fimo.exists():
                tasks.append((s_dir.name, target_fimo))

    # 2. Execute in Parallel
    print(f"[INFO] Found {len(tasks)} species. Running on {CORES} cores...")
    with Pool(processes=CORES) as pool:
        results = pool.map(count_motifs_in_species, tasks)

    # 3. Organize Results
    species_data = {s_name: counts for s_name, counts in results}
    ordered_species_list = [t[0] for t in tasks]

    # 4. Write TSV
    output_path = input_base / output_filename
    with open(output_path, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Motif_ID', 'Motif_Name'] + ordered_species_list)
        for mid in motif_ids:
            row = [mid, motif_info[mid]]
            for s_name in ordered_species_list:
                row.append(species_data[s_name].get(mid, 0))
            writer.writerow(row)
            
    print(f"[SUCCESS] Table generated: {output_path}")

def main():
    # Load Motif Metadata once
    motif_ids, motif_info = [], {}
    with open(ID_NAME_FILE, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            mid = row['Motif_ID']
            motif_ids.append(mid)
            motif_info[mid] = row['Motif_Name']

    # Run both jobs sequentially
    for job in JOBS:
        run_aggregation(job['input_dir'], job['output_name'], motif_ids, motif_info)

if __name__ == "__main__":
    main()