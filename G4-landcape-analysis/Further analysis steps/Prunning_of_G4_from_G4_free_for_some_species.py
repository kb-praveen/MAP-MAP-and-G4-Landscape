import csv
import gzip
import os
import shutil
from pathlib import Path

# --- CONFIGURATION ---
G4_IN_PROMOTER_BASE = Path('/scratch/kumarpraveen/G4_motif_analysis/G4_in_promoters')
FIMO_ENTIRE_BASE = Path('/scratch/kumarpraveen/G4_motif_analysis/Fimo_on_entire_promoters')
FIMO_G4_FREE_BASE = Path('/scratch/kumarpraveen/G4_motif_analysis/Fimo_on_G4_free_promoters')

# The mapping of (Clade, Long_Name_in_Entire) -> Short_Name_in_G4
MISSING_SPECIES = [
    ("Aves", "Anas_platyrhynchos_platyrhynchos", "Anas_platyrhynchos"),
    ("Aves", "Aquila_chrysaetos_chrysaetos", "Aquila_chrysaetos"),
    ("Reptiles", "Pseudonaja_textilis", "Pseudonaja_textlis") # Note the typo 'textlis' in your G4 folder
]

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

def prune_only(clade, entire_name, g4_name):
    print(f"\n[PRUNING] Clade: {clade} | {entire_name} using G4 data from {g4_name}")
    
    # Paths
    entire_fimo_gz = FIMO_ENTIRE_BASE / clade / entire_name / "fimo.tsv.gz"
    g4_fimo_file = G4_IN_PROMOTER_BASE / clade / g4_name / "fimo.tsv.gz"
    output_dir = FIMO_G4_FREE_BASE / clade / entire_name
    final_gz = output_dir / "fimo.tsv.gz"

    if not entire_fimo_gz.exists():
        print(f"  [ERROR] Entire FIMO file not found at {entire_fimo_gz}")
        return
    if not g4_fimo_file.exists():
        print(f"  [ERROR] G4 FIMO file not found at {g4_fimo_file}")
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load G4 Keys
    g4_keys = set()
    with gzip.open(g4_fimo_file, 'rt') as f:
        lines = (line for line in f if not line.startswith('#'))
        reader = csv.DictReader(lines, delimiter='\t')
        for row in reader:
            c_name, a_start, a_stop = parse_g4_info(row['sequence_name'], row['start'], row['stop'])
            g4_keys.add((row['motif_id'], c_name, a_start, a_stop, row['strand']))

    # 2. Prune and Write
    total, removed = 0, 0
    with gzip.open(entire_fimo_gz, 'rt') as f_in, gzip.open(final_gz, 'wt') as f_out:
        # Pass through comments
        line = f_in.readline()
        while line.startswith('#'):
            f_out.write(line)
            line = f_in.readline()
        
        # Header
        f_out.write(line)
        header = line.strip().split('\t')
        
        reader = csv.DictReader(f_in, fieldnames=header, delimiter='\t')
        writer = csv.DictWriter(f_out, fieldnames=header, delimiter='\t')
        
        for row in reader:
            total += 1
            # Comparison key
            if (row['motif_id'], row['sequence_name'], row['start'], row['stop'], row['strand']) in g4_keys:
                removed += 1
            else:
                writer.writerow(row)

    print(f"  [DONE] Total: {total:,} | Removed (G4-overlaps): {removed:,} | Saved: {total-removed:,}")

if __name__ == "__main__":
    for clade, entire, g4 in MISSING_SPECIES:
        prune_only(clade, entire, g4)
    print("\n[COMPLETE] Targeted pruning finished.")