import os
import subprocess
import shutil
import csv
import gzip
from pathlib import Path
from multiprocessing import Pool
from collections import defaultdict

# --- HPC CONFIGURATION ---
CORES = 48
FIMO_THRESH = "1e-4" 

MEME_FILE = Path('Path/G4_motif_analysis/Background_G4_free/JASPAR_vertebrate_MEME_Non_redundant_2026.txt')
ID_NAME_FILE = Path('Path/G4_motif_analysis/Background_G4_free/JASPAR_vertebrate_Non_redundant_2026_Ids_names.txt')
INPUT_BASE = Path('Path/G4_motif_analysis/Background_G4_free')
OUTPUT_BASE = Path('Path/G4_motif_analysis/Fimo_G4_free_background')

CLADE_ORDER = ['Pisces', 'Amphibians', 'Reptiles', 'Aves', 'Mammals']

def split_fasta(input_fasta, temp_dir, n_chunks):
    chunk_paths = [temp_dir / "chunk_{0}.fa".format(i) for i in range(n_chunks)]
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
    chunk_path, chunk_out_gz = args
    cmd = ["fimo", "--thresh", FIMO_THRESH, "--text", str(MEME_FILE), str(chunk_path)]
    try:
        with gzip.open(chunk_out_gz, 'wt') as f_gz:
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True)
            for line in process.stdout:
                f_gz.write(line)
            process.wait()
        if chunk_path.exists(): os.remove(chunk_path)
        return chunk_out_gz
    except Exception as e:
        print("Error in worker: {0}".format(e))
        return None

def process_species_background(clade, fasta_path):
    species_name = fasta_path.name.replace('_background.fa', '').replace('_background.fasta', '')
    species_dir = OUTPUT_BASE / clade / species_name
    final_gz = species_dir / "fimo.tsv.gz"

    # We assume if the file exists and is > 100 bytes (header size), it is done
    if final_gz.exists() and final_gz.stat().st_size > 100:
        return species_name

    temp_dir = species_dir / "temp_work"
    if temp_dir.exists(): shutil.rmtree(temp_dir)
    species_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    print("\n[RUNNING] {0} (Streaming to Gzip)...".format(species_name))
    chunk_files = split_fasta(fasta_path, temp_dir, CORES)
    tasks = [(cf, temp_dir / "chunk_{0}.tsv.gz".format(i)) for i, cf in enumerate(chunk_files)]
    
    with Pool(processes=CORES) as pool:
        chunk_gz_results = pool.map(run_fimo_worker, tasks)
    
    header_written = False
    with gzip.open(final_gz, 'wt') as fout:
        for gz_file in chunk_gz_results:
            if gz_file and gz_file.exists():
                with gzip.open(gz_file, 'rt') as fin:
                    for line in fin:
                        if line.startswith('motif_id'):
                            if not header_written:
                                fout.write(line)
                                header_written = True
                            continue
                        fout.write(line)
                os.remove(gz_file)

    shutil.rmtree(temp_dir)
    print("[SUCCESS] {0} -> fimo.tsv.gz".format(species_name))
    return species_name

def aggregate_counts_no_pandas():
    print("\n" + "="*60 + "\nAGGREGATING MOTIF COUNTS\n" + "="*60)
    motif_ids, motif_info = [], {}
    with open(ID_NAME_FILE, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            mid = row['Motif_ID']
            motif_ids.append(mid)
            motif_info[mid] = row['Motif_Name']

    species_data, ordered_species_list = {}, []
    for clade in CLADE_ORDER:
        clade_dir = OUTPUT_BASE / clade
        if not clade_dir.exists(): continue
        species_dirs = sorted([d for d in clade_dir.iterdir() if d.is_dir()])
        for s_dir in species_dirs:
            target_fimo = s_dir / "fimo.tsv.gz"
            if target_fimo.exists() and target_fimo.stat().st_size > 100:
                s_name = s_dir.name
                print("Counting: {0} | {1}".format(clade, s_name))
                ordered_species_list.append(s_name)
                counts = defaultdict(int)
                with gzip.open(target_fimo, 'rt') as f:
                    lines = (line for line in f if not line.startswith('#'))
                    reader = csv.DictReader(lines, delimiter='\t')
                    for row in reader: counts[row['motif_id']] += 1
                species_data[s_name] = counts

    output_file = OUTPUT_BASE / "Motif_counts_in_G4_free_regions_equivalent_genomewide_in_vertebrates.tsv"
    with open(output_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Motif_ID', 'Motif_Name'] + ordered_species_list)
        for mid in motif_ids:
            row = [mid, motif_info[mid]]
            for s_name in ordered_species_list:
                row.append(species_data[s_name].get(mid, 0))
            writer.writerow(row)
    print("\n[DONE] Final table: {0}".format(output_file))

if __name__ == "__main__":
    if not OUTPUT_BASE.exists(): OUTPUT_BASE.mkdir(parents=True, exist_ok=True)

    for clade in CLADE_ORDER:
        clade_input_dir = INPUT_BASE / clade / 'fasta'
        if not clade_input_dir.exists(): continue
        
        # --- RECONCILIATION PHASE ---
        fasta_files = list(clade_input_dir.glob("*_background.fa")) + list(clade_input_dir.glob("*_background.fasta"))
        to_process = []
        already_done = []

        for f in fasta_files:
            species_name = f.name.replace('_background.fa', '').replace('_background.fasta', '')
            check_file = OUTPUT_BASE / clade / species_name / "fimo.tsv.gz"
            
            # Check if file exists AND is larger than a standard header (~100 bytes)
            if check_file.exists() and check_file.stat().st_size > 100:
                already_done.append(species_name)
            else:
                to_process.append(f)

        print("\n{0}\nCLADE: {1}\n{0}".format('#'*60, clade.upper()))
        print("Status: {0} Finished, {1} Pending".format(len(already_done), len(to_process)))
        
        if len(already_done) > 0:
            print("Already Done: {0}...".format(", ".join(already_done[:5])))

        # 1. Run FIMO only on pending species
        for f in sorted(to_process, key=lambda x: os.path.getsize(x)):
            process_species_background(clade, f)

    # 2. Final Aggregation
    aggregate_counts_no_pandas()