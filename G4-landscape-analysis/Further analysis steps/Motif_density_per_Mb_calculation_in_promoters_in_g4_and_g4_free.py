import csv
import gzip
import os
import subprocess
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool

# --- CONFIGURATION ---
CORES = 48
ID_NAME_FILE = Path('Path/G4_motif_analysis/Background_G4_free/JASPAR_vertebrate_Non_redundant_2026_Ids_names.txt')
CLADE_ORDER = ['Pisces', 'Amphibians', 'Reptiles', 'Aves', 'Mammals']
RESULTS_DIR = Path('Path/G4_motif_analysis/Results')
ENTIRE_PROM_FASTA = Path('Path/G4_motif_analysis/Entire_promoter_fasta')
G4_PROM_FASTA = Path('Path/G4_motif_analysis/Promoter_G4_fasta')

JOBS = [
    {
        "name": "G4_free_promoter",
        "input_dir": Path('Path/G4_motif_analysis/Fimo_on_G4_free_promoters'),
        "count_file": "Motif_counts_in_G4_free_promoters_raw.tsv"
    },
    {
        "name": "G4_in_promoter",
        "input_dir": Path('Path/G4_motif_analysis/G4_in_promoters'),
        "count_file": "Motif_counts_in_G4_containing_promoters_raw.tsv"
    }
]

def get_fasta_length(fasta_path):
    total_len = 0
    try:
        f = gzip.open(fasta_path, 'rt') if fasta_path.suffix == '.gz' else open(fasta_path, 'r')
        for line in f:
            if not line.startswith('>'):
                total_len += len(line.strip())
        f.close()
    except Exception as e:
        print(f"Error reading fasta {fasta_path}: {e}")
    return total_len

def process_fimo_counts(args):
    s_name, fimo_path = args
    counts = defaultdict(int)
    try:
        if not fimo_path.exists() or fimo_path.stat().st_size < 100:
            return s_name, {}
        with gzip.open(fimo_path, 'rt') as f:
            lines = (line for line in f if not line.startswith('#'))
            reader = csv.DictReader(lines, delimiter='\t')
            for row in reader:
                counts[row['motif_id']] += 1
        return s_name, counts
    except Exception as e:
        return s_name, {}

def main():
    print("="*70 + "\nSTEP 1: CALCULATING SEARCH SPACE FROM FASTA\n" + "="*70)
    final_lengths = defaultdict(lambda: {'Total': 0, 'G4': 0, 'Free': 0})

    for clade in CLADE_ORDER:
        for f in (ENTIRE_PROM_FASTA / clade).glob("*.fa*"):
            s_name = f.name.split('.')[0]
            final_lengths[s_name]['Total'] = get_fasta_length(f)
        for f in (G4_PROM_FASTA / clade).glob("*.fa*"):
            s_name = f.name.split('.')[0]
            final_lengths[s_name]['G4'] = get_fasta_length(f)

    for s in final_lengths:
        final_lengths[s]['Free'] = final_lengths[s]['Total'] - final_lengths[s]['G4']

    print("\n" + "="*70 + "\nSTEP 2: AGGREGATING MOTIF COUNTS\n" + "="*70)
    motif_ids, motif_info = [], {}
    with open(ID_NAME_FILE, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            mid = row['Motif_ID']; motif_ids.append(mid)
            motif_info[mid] = row['Motif_Name']

    for job in JOBS:
        tasks = []
        for clade in CLADE_ORDER:
            clade_dir = job['input_dir'] / clade
            if not clade_dir.exists(): continue
            for s_dir in sorted([d for d in clade_dir.iterdir() if d.is_dir()]):
                if (s_dir / "fimo.tsv.gz").exists():
                    tasks.append((s_dir.name, s_dir / "fimo.tsv.gz"))

        with Pool(processes=CORES) as pool:
            results = pool.map(process_fimo_counts, tasks)

        species_counts = {s: counts for s, counts in results}
        ordered_species = [t[0] for t in tasks]

        output_count_path = job['input_dir'] / job['count_file']
        with open(output_count_path, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['Motif_ID', 'Motif_Name'] + ordered_species)
            for mid in motif_ids:
                row = [mid, motif_info[mid]] + [species_counts[s].get(mid, 0) for s in ordered_species]
                writer.writerow(row)

    summary_path = RESULTS_DIR / 'Promoter_Search_Space_Summary.tsv'
    with open(summary_path, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Species', 'G4_BP', 'Free_BP', 'Total_BP', 'G4_Mb', 'Free_Mb'])
        for s, lts in final_lengths.items():
            writer.writerow([s, lts['G4'], lts['Free'], lts['Total'], lts['G4']/1e6, lts['Free']/1e6])

    print("\n" + "="*70 + "\nSTEP 3: R - ROBUST NORMALIZATION\n" + "="*70)
    
    r_script = """
library(readr)
library(dplyr)
library(stringr)

base <- "Path/G4_motif_analysis"
lens <- read_tsv(file.path(base, "Results/Promoter_Search_Space_Summary.tsv"), show_col_types = FALSE)

# Robust matching helper
get_mb <- function(target_species, mb_col) {
  # Try exact match first
  idx <- match(target_species, lens$Species)
  # Fallback: match first two components (Genus_species)
  if(is.na(idx)) {
    short_target <- word(target_species, 1, 2, sep="_")
    short_lens <- word(lens$Species, 1, 2, sep="_")
    idx <- match(short_target, short_lens)
  }
  return(lens[[mb_col]][idx])
}

norm_mb <- function(in_f, mb_col, out_dir, out_f) {
  if(!file.exists(in_f)) return(NULL)
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  d <- read_tsv(in_f, show_col_types = FALSE)
  motif_meta <- d[,1:2]
  mat <- as.matrix(d[,-c(1,2)])
  
  # Calculate divisors using robust matching
  divs <- sapply(colnames(mat), get_mb, mb_col = mb_col)
  divs[divs == 0] <- NA
  
  dens <- sweep(mat, 2, divs, "/")
  write_tsv(cbind(motif_meta, as.data.frame(dens)), file.path(out_dir, out_f))
  cat(paste("[SUCCESS] Saved:", out_f, "\\n"))
}

norm_mb(file.path(base, "G4_in_promoters/Motif_counts_in_G4_containing_promoters_raw.tsv"), 
        "G4_Mb", file.path(base, "Results/G4_in_promoter"), "Motif_Density_per_Mb_G4_Promoters.tsv")
norm_mb(file.path(base, "Fimo_on_G4_free_promoters/Motif_counts_in_G4_free_promoters_raw.tsv"), 
        "Free_Mb", file.path(base, "Results/G4_free_promoter"), "Motif_Density_per_Mb_G4_free_Promoters.tsv")
"""
    with open("robust_norm.R", "w") as f: f.write(r_script)
    subprocess.run(["Rscript", "robust_norm.R"])
    os.remove("robust_norm.R")

if __name__ == "__main__":
    main()