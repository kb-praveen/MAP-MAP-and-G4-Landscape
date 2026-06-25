import os
import subprocess
import shutil
import csv
import gzip
from pathlib import Path
from multiprocessing import Pool

# --- HPC CONFIGURATION ---
CORES = 48
FIMO_THRESH = "1e-4" 
MEME_FILE = Path('Path/G4_motif_analysis/Background_G4_free/JASPAR_vertebrate_MEME_Non_redundant_2026.txt')

# Base Directories
FASTA_BASE = Path('Path/G4_motif_analysis/Entire_promoter_fasta')
G4_IN_PROMOTER_BASE = Path('Path/G4_motif_analysis/G4_in_promoters')
FIMO_ENTIRE_BASE = Path('Path/G4_motif_analysis/Fimo_on_entire_promoters')
FIMO_G4_FREE_BASE = Path('Path/G4_motif_analysis/Fimo_on_G4_free_promoters')

CLADES = ['Pisces', 'Amphibians', 'Reptiles', 'Aves', 'Mammals']

# --- COORDINATE & FILE HELPERS ---

def parse_g4_info(name, motif_start, motif_stop):
    """Translates G4-fragment coordinates to absolute genomic coordinates."""
    if not name: return None, 0, 0
    parts = name.split(':')
    if len(parts) >= 3:
        chr_name = parts[0]
        try:
            frag_start = int(parts[1].split('-')[0])
            g4_offset = int(parts[2].split('-')[0])
            abs_start = str(frag_start + g4_offset + int(motif_start) - 1)
            abs_stop = str(frag_start + g4_offset + int(motif_stop) - 1)
            return chr_name, abs_start, abs_stop
        except: return chr_name, motif_start, motif_stop
    return name, motif_start, motif_stop

def split_fasta(input_fasta, temp_dir, n_chunks):
    chunk_paths = [temp_dir / f"chunk_{i}.fa" for i in range(n_chunks)]
    chunk_handles = [open(p, 'w') for p in chunk_paths]
    with open(input_fasta, 'r') as f:
        idx = 0
        entry = []
        for line in f:
            if line.startswith('>'):
                if entry:
                    chunk_handles[idx % n_chunks].write("".join(entry))
                    idx += 1
                entry = [line]
            else: entry.append(line)
        if entry: chunk_handles[idx % n_chunks].write("".join(entry))
    for h in chunk_handles: h.close()
    return chunk_paths

def run_fimo_worker(args):
    chunk_path, chunk_out_gz = args
    cmd = ["fimo", "--thresh", FIMO_THRESH, "--text", str(MEME_FILE), str(chunk_path)]
    try:
        with gzip.open(chunk_out_gz, 'wt') as f_gz:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True)
            for line in proc.stdout:
                f_gz.write(line)
            proc.wait()
        if chunk_path.exists(): os.remove(chunk_path)
        return chunk_out_gz
    except: return None

# --- CORE PROCESSING LOGIC ---

def process_species(clade, fasta_path):
    # Extract species name: e.g. Erpetoichthys_calabaricus
    species_name = fasta_path.name.split('.')[0]
    
    # Setup paths
    species_fimo_dir = FIMO_ENTIRE_BASE / clade / species_name
    entire_fimo_gz = species_fimo_dir / "fimo.tsv.gz"
    g4_fimo_file = G4_IN_PROMOTER_BASE / clade / species_name / "fimo.tsv.gz"
    output_dir = FIMO_G4_FREE_BASE / clade / species_name
    final_gz = output_dir / "fimo.tsv.gz"

    # Skip if final output already exists
    if final_gz.exists() and final_gz.stat().st_size > 1000:
        print(f"  [SKIP] {species_name} already completed.")
        return

    print(f"\n{'#'*40}\nProcessing: {species_name} ({clade})\n{'#'*40}")

    # PHASE 1: PARALLEL FIMO
    if species_fimo_dir.exists(): shutil.rmtree(species_fimo_dir)
    species_fimo_dir.mkdir(parents=True, exist_ok=True)
    
    chunks = split_fasta(fasta_path, species_fimo_dir, CORES)
    tasks = [(c, species_fimo_dir / f"c_{i}.gz") for i, c in enumerate(chunks)]
    
    with Pool(processes=CORES) as pool:
        results = pool.map(run_fimo_worker, tasks)

    header_written = False
    with gzip.open(entire_fimo_gz, 'wt') as fout:
        for r_gz in results:
            if r_gz and r_gz.exists():
                with gzip.open(r_gz, 'rt') as fin:
                    for line in fin:
                        if line.startswith('motif_id'):
                            if not header_written:
                                fout.write(line); header_written = True
                            continue
                        fout.write(line)
                os.remove(r_gz)

    # PHASE 2: PRUNING
    if not g4_fimo_file.exists():
        print(f"  [WARN] G4 FIMO file missing for {species_name}. Skipping pruning.")
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    g4_keys = set()
    with gzip.open(g4_fimo_file, 'rt') as f:
        header = None
        for line in f:
            if line.startswith('motif_id'):
                header = line.strip().split('\t'); break
        if not header: return
        reader = csv.DictReader(f, fieldnames=header, delimiter='\t')
        for row in reader:
            if not row['motif_id'] or row['motif_id'].startswith('#'): continue
            c_name, a_start, a_stop = parse_g4_info(row['sequence_name'], row['start'], row['stop'])
            g4_keys.add((row['motif_id'], c_name, a_start, a_stop, row['strand']))

    total, removed = 0, 0
    with gzip.open(entire_fimo_gz, 'rt') as f_in, gzip.open(final_gz, 'wt') as f_out:
        header = None
        for line in f_in:
            if line.startswith('motif_id'):
                header = line.strip().split('\t'); f_out.write(line); break
        
        if header:
            reader = csv.DictReader(f_in, fieldnames=header, delimiter='\t')
            writer = csv.DictWriter(f_out, fieldnames=header, delimiter='\t')
            for row in reader:
                total += 1
                if (row['motif_id'], row['sequence_name'], row['start'], row['stop'], row['strand']) in g4_keys:
                    removed += 1
                else: writer.writerow(row)

    print(f"  [DONE] {species_name}: Total {total:,} | Removed {removed:,}")

# --- MAIN LOOP ---

if __name__ == "__main__":
    for clade in CLADES:
        clade_fasta_dir = FASTA_BASE / clade
        if not clade_fasta_dir.exists(): continue
        
        print(f"\n{'='*60}\nCLADE: {clade.upper()}\n{'='*60}")
        
        # Get all FASTA files in clade
        fastas = sorted(list(clade_fasta_dir.glob("*.fa")) + list(clade_fasta_dir.glob("*.fasta")))
        
        for f in fastas:
            process_species(clade, f)

    print("\n[COMPLETE] All clades processed successfully.")