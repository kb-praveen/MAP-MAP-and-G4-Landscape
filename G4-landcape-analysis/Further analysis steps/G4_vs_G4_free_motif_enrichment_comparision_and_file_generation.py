import os
import csv
import math
import sys
from pathlib import Path
from collections import defaultdict

# --- CONFIGURATION ---
CORES = 48
G4_RAW_COUNTS = Path('/scratch/kumarpraveen/G4_motif_analysis/Fimo_in_G4s_genomewide/Motif_counts_in_G4_regions_genomewide_in_vertebrates.tsv')
BG_RAW_COUNTS = Path('/scratch/kumarpraveen/G4_motif_analysis/Fimo_G4_free_background/Motif_counts_in_G4_free_regions_equivalent_genomewide_in_vertebrates.tsv')

METADATA_FILE = Path('/scratch/kumarpraveen/G4_motif_analysis/Table_S1.tsv')
OUTPUT_DIR = Path('/scratch/kumarpraveen/G4_motif_analysis/Statistical_Analysis')

def calculate_bh_correction(p_values):
    """Multiple testing correction using Benjamini-Hochberg."""
    n = len(p_values)
    if n == 0: return []
    sorted_p = sorted(enumerate(p_values), key=lambda x: x[1])
    adj_p = [0] * n
    prev_bh = 1.0
    for i in range(n - 1, -1, -1):
        idx, p = sorted_p[i]
        rank = i + 1
        bh = p * n / rank
        bh = min(bh, prev_bh)
        adj_p[idx] = bh
        prev_bh = bh
    return adj_p

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # 1. Load Metadata (Species -> G4_Predicted as N)
    print(f"[VERBOSE] Loading total counts (N) from {METADATA_FILE}...")
    species_N = {}
    with open(METADATA_FILE, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            species_N[row['Species'].replace(' ', '_')] = int(row['G4_Predicted'])

    # 2. Load Raw Counts
    def load_matrix(path):
        data = defaultdict(dict)
        motifs = []
        with open(path, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)
            species_list = header[2:]
            for row in reader:
                mid, mname = row[0], row[1]
                motifs.append((mid, mname))
                for i, count in enumerate(row[2:]):
                    data[species_list[i]][mid] = int(count)
        return data, motifs, species_list

    g4_data, motifs, processed_species = load_matrix(G4_RAW_COUNTS)
    bg_data, _, _ = load_matrix(BG_RAW_COUNTS)

    # 3. Statistical Phase
    total_tests = len(motifs) * len(processed_species)
    print(f"[VERBOSE] Running {total_tests:,} statistical tests...")
    
    all_p_values = []
    results_flat = []
    
    # Pre-initialize matrices for grid output
    # matrix[motif_id][species_name]
    fdr_matrix = defaultdict(dict)
    fold_matrix = defaultdict(dict)

    count = 0
    for mid, mname in motifs:
        for s_name in processed_species:
            N = species_N.get(s_name)
            if not N: continue
            
            a = g4_data[s_name].get(mid, 0)
            b = bg_data[s_name].get(mid, 0)
            c = max(0, N - a)
            d = max(0, N - b)

            # Odds Ratio with pseudocount
            ap, bp, cp, dp = a+0.5, b+0.5, c+0.5, d+0.5
            odds_ratio = (ap / cp) / (bp / dp)
            
            try:
                se = math.sqrt(1/ap + 1/bp + 1/cp + 1/dp)
                z = math.log(odds_ratio) / se
                p = math.erfc(abs(z) / math.sqrt(2))
            except: 
                p = 1.0

            all_p_values.append(p)
            results_flat.append({
                's': s_name, 'mid': mid, 'mname': mname, 
                'a': a, 'b': b, 'c': c, 'd': d, 'N': N,
                'fold': odds_ratio, 'p': p
            })
            
            count += 1
            if count % 10000 == 0:
                sys.stdout.write(f"\r[PROGRESS] Calculating Stats: {(count/total_tests)*100:.1f}%")
                sys.stdout.flush()

    sys.stdout.write("\n")
    print(f"[VERBOSE] Applying FDR correction...")
    adj_p_values = calculate_bh_correction(all_p_values)
    
    # 4. Generate Output Files
    summary_path = OUTPUT_DIR / "Significant_Findings_Summary_FDR.tsv"
    fdr_path = OUTPUT_DIR / "FDR_Adjusted_P_Matrix.tsv"
    fold_path = OUTPUT_DIR / "Fold_Enrichment_Matrix.tsv"

    print(f"[VERBOSE] Writing final matrices and summary...")
    
    with open(summary_path, 'w') as f_sum:
        f_sum.write("Species\tMotif_ID\tMotif_Name\tCounts_G4(a)\tCounts_BG(b)\tTotal(N)\tG4_NoMotif(c)\tBG_NoMotif(d)\tOdds_Ratio\tP_Value\tFDR_Adj_P\tResult\n")
        
        for i, res in enumerate(results_flat):
            fdr = adj_p_values[i]
            # Store in matrix dictionary for Step 5
            fdr_matrix[res['mid']][res['s']] = fdr
            fold_matrix[res['mid']][res['s']] = res['fold']
            
            if fdr < 0.05:
                res_type = "G4_Enriched" if res['fold'] > 1 else "BG_Enriched"
                f_sum.write(f"{res['s']}\t{res['mid']}\t{res['mname']}\t{res['a']}\t{res['b']}\t{res['N']}\t{res['c']}\t{res['d']}\t{res['fold']:.4f}\t{res['p']:.2e}\t{fdr:.2e}\t{res_type}\n")

    # 5. Write the Square Matrices
    header = "Motif_ID\tMotif_Name\t" + "\t".join(processed_species) + "\n"
    
    with open(fdr_path, 'w') as f_fdr, open(fold_path, 'w') as f_fold:
        f_fdr.write(header)
        f_fold.write(header)
        
        for mid, mname in motifs:
            fdr_vals = [str(fdr_matrix[mid][s]) for s in processed_species]
            fold_vals = [f"{fold_matrix[mid][s]:.4f}" for s in processed_species]
            
            f_fdr.write(f"{mid}\t{mname}\t" + "\t".join(fdr_vals) + "\n")
            f_fold.write(f"{mid}\t{mname}\t" + "\t".join(fold_vals) + "\n")

    print(f"\n[SUCCESS] Files generated in {OUTPUT_DIR}:")
    print(f"1. Significant_Findings_Summary_FDR.tsv (Filtered)")
    print(f"2. FDR_Adjusted_P_Matrix.tsv (Full Grid)")
    print(f"3. Fold_Enrichment_Matrix.tsv (Full Grid)")

if __name__ == "__main__":
    main()