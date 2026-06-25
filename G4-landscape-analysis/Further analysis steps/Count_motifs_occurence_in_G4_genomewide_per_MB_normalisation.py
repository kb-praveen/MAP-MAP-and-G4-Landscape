import csv
import os
from pathlib import Path

# --- CONFIGURATION ---
METADATA_FILE = Path('Path/G4_motif_analysis/Table_S1.tsv')

# Input paths
G4_COUNTS_FILE = Path('Path/G4_motif_analysis/Fimo_in_G4s_genomewide/Motif_counts_in_G4_regions_genomewide_in_vertebrates.tsv')
BG_COUNTS_FILE = Path('Path/G4_motif_analysis/Fimo_G4_free_background/Motif_counts_in_G4_free_regions_equivalent_genomewide_in_vertebrates.tsv')

# Output directories
RESULTS_BASE = Path('Path/G4_motif_analysis/Results')
G4_OUT_DIR = RESULTS_BASE / "Genomewide_comparison"
BG_OUT_DIR = RESULTS_BASE / "Genomewide_G4_free"

# Create directories
G4_OUT_DIR.mkdir(parents=True, exist_ok=True)
BG_OUT_DIR.mkdir(parents=True, exist_ok=True)

def normalize_matrix_per_mb(input_path, output_path, g4_totals):
    """Function to read a matrix, apply normalization per Million, and save result."""
    if not input_path.exists():
        print(f"[ERROR] Input file not found: {input_path}")
        return

    print(f"[VERBOSE] Normalizing: {input_path.name} to Per Mb Density")
    
    with open(input_path, 'r') as f_in, open(output_path, 'w', newline='') as f_out:
        reader = csv.reader(f_in, delimiter='\t')
        writer = csv.writer(f_out, delimiter='\t')
        
        header = next(reader)
        writer.writerow(header)
        species_list = header[2:]
        
        for row in reader:
            m_id, m_name = row[0], row[1]
            new_row = [m_id, m_name]
            
            raw_counts = row[2:]
            for i, count in enumerate(raw_counts):
                species = species_list[i]
                n_g4s = g4_totals.get(species, 0)
                
                if n_g4s > 0:
                    # NORMALIZATION: (Count * 1,000,000) / G4_Predicted
                    density_per_mb = (float(count) * 1000000.0) / n_g4s
                    new_row.append(round(density_per_mb, 4))
                else:
                    new_row.append(0.0)
            writer.writerow(new_row)

def main():
    print("="*60 + "\nGENOME-WIDE DENSITY NORMALIZATION (PER 10^6 SEQUENCES)\n" + "="*60)

    # 1. Load Metadata (Species -> G4_Predicted)
    g4_totals = {}
    print(f"[VERBOSE] Loading G4 predicted counts from {METADATA_FILE}...")
    with open(METADATA_FILE, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # Clean species name to match underscores in count files
            species_name = row['Species'].replace(' ', '_')
            try:
                g4_totals[species_name] = float(row['G4_Predicted'])
            except ValueError:
                continue

    # 2. Process G4 Regions (Density per Million G4s)
    g4_output = G4_OUT_DIR / "Motif_Density_per_Mb_G4s_Matrix.tsv"
    normalize_matrix_per_mb(G4_COUNTS_FILE, g4_output, g4_totals)

    # 3. Process G4-free Background Regions (Density per Million Background Seqs)
    bg_output = BG_OUT_DIR / "Motif_Density_per_Mb_G4free_Matrix.tsv"
    normalize_matrix_per_mb(BG_COUNTS_FILE, bg_output, g4_totals)

    print("\n" + "="*60)
    print(f"[SUCCESS] G4 results saved: {g4_output}")
    print(f"[SUCCESS] BG results saved: {bg_output}")
    print("="*60)

if __name__ == "__main__":
    main()