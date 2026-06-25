import os
import subprocess
import shutil
import csv
from pathlib import Path
from multiprocessing import Pool
from collections import defaultdict

# --- HPC CONFIGURATION ---
CORES = 48
MEME_FILE = Path('/scratch/kumarpraveen/G4_motif_analysis/Background_G4_free/JASPAR_vertebrate_MEME_Non_redundant_2026.txt')
ID_NAME_FILE = Path('/scratch/kumarpraveen/G4_motif_analysis/Background_G4_free/JASPAR_vertebrate_Non_redundant_2026_Ids_names.txt')
INPUT_BASE = Path('/scratch/kumarpraveen/G4_motif_analysis/Background_G4_free')
OUTPUT_BASE = Path('/scratch/kumarpraveen/G4_motif_analysis/Fimo_G4_free_background')

# Taxon order for the final table header
CLADE_ORDER = ['Pisces', 'Amphibians', 'Reptiles', 'Aves', 'Mammals']

def split_fasta(input_fasta, temp_dir, n_chunks):
    chunk_paths = [temp_dir / f"chunk_{i}.fa" for i in range(n_chunks)]
    chunk_handles = [open(p, 'w') for p in chunk_paths]
    with open(input_fasta, 'r') as f:
        chunk_idx = 0
        current_entry = []
        for line in f:
            if line.startswith('>'):
                if current_entry:
                    chunk_handles[chunk_idx % n_chunks].write("".join(current_entry))
                    chunk_idx += 1
                current_entry = [line]
            else:
                current_entry.append(line)
        if current_entry:
            chunk_handles[chunk_idx % n_chunks].write("".join(current_entry))
    for h in chunk_handles: h.close()
    return chunk_paths

def run_fimo_worker(args):
    chunk_path, chunk_out_dir = args
    cmd = ["fimo", "--o", str(chunk_out_dir), str(MEME_FILE), str(chunk_path)]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return chunk_out_dir / "fimo.tsv"
    except:
        return None

def process_species_background(clade, fasta_path):
    species_name = fasta_path.name.replace('_background.fa', '').replace('_background.fasta', '')
    species_dir = OUTPUT_BASE / clade / species_name
    final_tsv = species_dir / "fimo.tsv"

    if final_tsv.exists() and final_tsv.stat().st_size > 0:
        return species_name

    temp_dir = species_dir / "temp_work"
    if temp_dir.exists(): shutil.rmtree(temp_dir)
    species_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n[RUNNING] {species_name}...")
    chunk_files = split_fasta(fasta_path, temp_dir, CORES)
    tasks = [(cf, temp_dir / f"out_{i}") for i, cf in enumerate(chunk_files)]
    with Pool(processes=CORES) as pool:
        tsv_results = pool.map(run_fimo_worker, tasks)
    
    header_written = False
    with open(final_tsv, 'w') as fout:
        for tsv in tsv_results:
            if tsv and tsv.exists():
                with open(tsv, 'r') as fin:
                    for line in fin:
                        if line.startswith(('#', 'motif_id')) and header_written: continue
                        if line.startswith('motif_id'): header_written = True
                        fout.write(line)
    shutil.rmtree(temp_dir)
    return species_name

def aggregate_counts_no_pandas():
    print("\n" + "="*60 + "\nAGGREGATING MOTIF COUNTS (Pure Python Mode)\n" + "="*60)
    
    # 1. Load Motif ID to Name mapping
    motif_ids = []
    motif_info = {} # ID -> Name
    with open(ID_NAME_FILE, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            mid = row['Motif_ID']
            motif_ids.append(mid)
            motif_info[mid] = row['Motif_Name']

    # 2. Prepare Data Storage
    # Structure: species_data[species_name][motif_id] = count
    species_data = {}
    ordered_species_list = []

    for clade in CLADE_ORDER:
        clade_dir = OUTPUT_BASE / clade
        if not clade_dir.exists(): continue
        
        species_dirs = sorted([d for d in clade_dir.iterdir() if d.is_dir()])
        for s_dir in species_dirs:
            s_name = s_dir.name
            fimo_file = s_dir / "fimo.tsv"
            
            if fimo_file.exists():
                print(f"Counting: {clade} | {s_name}")
                ordered_species_list.append(s_name)
                counts = defaultdict(int)
                
                with open(fimo_file, 'r') as f:
                    # Skip comment lines and use DictReader for the rest
                    lines = (line for line in f if not line.startswith('#'))
                    reader = csv.DictReader(lines, delimiter='\t')
                    for row in reader:
                        counts[row['motif_id']] += 1
                
                species_data[s_name] = counts

    # 3. Write Output TSV
    output_file = OUTPUT_BASE / "Motif_counts_in_G4_free_regions_equivalent_genomewide_in_vertebrates.tsv"
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        
        # Header: Motif_ID, Motif_Name, Species1, Species2...
        header = ['Motif_ID', 'Motif_Name'] + ordered_species_list
        writer.writerow(header)
        
        # Rows
        for mid in motif_ids:
            row = [mid, motif_info[mid]]
            for s_name in ordered_species_list:
                row.append(species_data[s_name].get(mid, 0))
            writer.writerow(row)

    print(f"\n[SUCCESS] Master counts file written to: {output_file}")

if __name__ == "__main__":
    if not OUTPUT_BASE.exists():
        OUTPUT_BASE.mkdir(parents=True, exist_ok=True)

    # 1. Run FIMO
    for clade in CLADE_ORDER:
        # Note: Check if your path uses 'Amphibians' (script) or 'Amphibia' (your text)
        # Using a list comprehension to catch both
        clade_input_dir = INPUT_BASE / clade / 'fasta'
        if not clade_input_dir.exists():
            continue
        
        print(f"\n{'#'*60}\nCLADE: {clade.upper()}\n{'#'*60}")
        fasta_files = list(clade_input_dir.glob("*_background.fa")) + list(clade_input_dir.glob("*_background.fasta"))
        for f in sorted(fasta_files, key=lambda x: os.path.getsize(x)):
            process_species_background(clade, f)

    # 2. Aggregate
    aggregate_counts_no_pandas()