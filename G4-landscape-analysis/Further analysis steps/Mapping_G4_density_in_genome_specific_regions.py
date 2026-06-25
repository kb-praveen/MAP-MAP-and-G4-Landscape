--------------------------------------------------------------------------------------------
#Run this for UTRs - normalised by Genome GC-content
--------------------------------------------------------------------------------------------
import os
import subprocess
import pandas as pd
import pathlib
import tempfile
import sys

# --- 1. Configuration ---
BASE_DIR = "path to the directory /Genome_specific_regions/Genome_specific_coordinates"
G4_DIR = "path to the directory /Genome_specific_regions/Predicted_G4_bed_genomes"
METADATA_PATH = "path to the directory /Genome_specific_regions/Metadata/Genomes.tsv"
MASTER_OUTPUT = "path to the directory /Genome_specific_regions/Output"

FEATURES_PRIORITY = [
    ("5UTR", "5UTR_regions.bed"),
    ("3UTR", "3UTR_regions.bed")
]

# --- 2. Load Metadata (Column 3 Check) ---
meta = pd.read_csv(METADATA_PATH, sep='\t')
# Column 0: Species, Column 2: Genome_GC_content (Index 2 is the 3rd column)
gc_map = dict(zip(meta.iloc[:, 0], meta.iloc[:, 2]))

print(f"\n{'='*80}")
print("🚀 G4 METAGENE PIPELINE: UTR SEQUENTIAL MODE")
print(f"Normalizing by Genome GC Content (Column 3: {meta.columns[2]})")
print(f"{'='*80}")

# --- 3. Main Loop ---
for feat_name, feat_file in FEATURES_PRIORITY:
    print(f"\n💎 CURRENT FEATURE: {feat_name.upper()}")
    
    feat_dir = os.path.join(MASTER_OUTPUT, feat_name)
    os.makedirs(feat_dir, exist_ok=True)
    out_path = os.path.join(feat_dir, f"All_Vertebrates_{feat_name}_Scaled_G4.tsv")
    
    # Reset/Initialize File with Header
    with open(out_path, 'w', buffering=1) as f_out: # buffering=1 forces line-by-line writing
        f_out.write("bin_num\tcount\tnormalized_density\tspecies\tclade\n")
    
    taxa_list = ["Amphibians", "Aves", "Reptiles", "Pisces", "Mammals"]
    
    for taxon in taxa_list:
        taxon_path = pathlib.Path(BASE_DIR) / taxon
        g4_path = pathlib.Path(G4_DIR) / taxon
        if not taxon_path.exists(): continue
            
        species_dirs = [d for d in taxon_path.iterdir() if d.is_dir() and not d.name.startswith('.')]
        print(f"\n📂 Clade: {taxon} ({len(species_dirs)} species)")
        
        for idx, s_dir in enumerate(species_dirs, 1):
            species_name = s_dir.name
            standard_name = species_name.replace(" ", "_")
            g4_file = g4_path / f"{standard_name}_G4.bed"
            coord_file = s_dir / feat_file
            
            if not g4_file.exists() or not coord_file.exists():
                print(f"    → [{idx}/{len(species_dirs)}] {species_name}: SKIP (Missing Bed)")
                continue
            
            try:
                # Get line count for progress
                total_lines = int(subprocess.check_output(f"wc -l < '{coord_file}'", shell=True))
                if total_lines == 0: continue

                with tempfile.TemporaryDirectory() as tmpdirname:
                    tmp_clean = os.path.join(tmpdirname, "clean.bed")
                    tmp_bins = os.path.join(tmpdirname, "bins.bed")
                    tmp_counts = os.path.join(tmpdirname, "counts.bed")
                    
                    # 
                    # Filter and scale
                    subprocess.run(f"awk '($3-$2) >= 100 && $2 >= 0' '{coord_file}' > {tmp_clean}", shell=True)
                    if os.path.getsize(tmp_clean) == 0:
                        print(f"    → {species_name}: SKIP (Too short)")
                        continue

                    subprocess.run(f"bedtools makewindows -b {tmp_clean} -n 100 -i srcwinnum > {tmp_bins}", shell=True, capture_output=True)
                    subprocess.run(f"bedtools intersect -c -a {tmp_bins} -b '{g4_file}' > {tmp_counts}", shell=True, capture_output=True)
                    
                    counts = pd.read_csv(tmp_counts, sep='\t', header=None, names=['chr', 'start', 'end', 'id_info', 'count'], dtype={0: str})
                    counts['bin_num'] = counts['id_info'].str.extract(r'_(\d+)$').astype(int)
                    
                    # Strand alignment
                    strands = pd.read_csv(tmp_clean, sep='\t', header=None, usecols=[3, 5], names=['feat_id', 'strand'], dtype={0: str})
                    counts['feat_id'] = counts['id_info'].str.replace(r'_\d+$', '', regex=True)
                    counts = counts.merge(strands, on='feat_id')
                    
                    # 
                    counts.loc[counts['strand'] == '-', 'bin_num'] = 101 - counts['bin_num']
                    
                    # Normalization: Space/Underscore alignment
                    lookup_name = species_name.replace("_", " ")
                    species_gc_val = gc_map.get(lookup_name, 50.0)
                    species_gc_ratio = float(species_gc_val) / 100.0

                    summary = counts.groupby('bin_num')['count'].mean().reset_index()
                    summary['normalized_density'] = summary['count'] / species_gc_ratio
                    summary['species'] = species_name
                    summary['clade'] = taxon
                    
                    # Append and hard-sync to disk
                    with open(out_path, 'a') as f:
                        summary.to_csv(f, sep='\t', header=False, index=False)
                        f.flush()
                        os.fsync(f.fileno())
                    
                    # Print exact Norm ratio for verification
                    sys.stdout.write(f"\r    → {species_name}: 100.0% ✅ DONE (Norm: {species_gc_ratio:.4f})")
                    sys.stdout.flush()
                    print()

            except Exception as e:
                print(f"    → {species_name}: ❌ ERROR: {e}")

print(f"\n{'='*80}\n🎉 UTR PROCESSING COMPLETE\n{'='*80}")


--------------------------------------------------------------------------------------------
##Run this for Introns and Entire Gene  - normalised by Genome GC-content
--------------------------------------------------------------------------------------------

import os
import subprocess
import pandas as pd
import pathlib
import tempfile
import sys

# --- 1. Configuration ---
BASE_DIR = "path to the directory /Genome_specific_regions/Genome_specific_coordinates"
G4_DIR = "path to the directory /Genome_specific_regions/Predicted_G4_bed_genomes"
METADATA_PATH = "path to the directory /Genome_specific_regions/Metadata/Genomes.tsv"
MASTER_OUTPUT = "path to the directory /Genome_specific_regions/Output"

FEATURES_PRIORITY = [
    ("Intron", "Intron_coordinates.bed"),
    ("Gene", "Gene_Body_coordinates.bed")
]

BATCH_SIZE = 20000 

# --- 2. Load Metadata (Standardized to Column 3) ---
meta = pd.read_csv(METADATA_PATH, sep='\t')
# Standardizing: Species is Column 0, Genome_GC_content is Column 2 (Index 2)
gc_map = dict(zip(meta.iloc[:, 0], meta.iloc[:, 2]))

print(f"\n{'='*80}")
print("🚀 G4 HEAVY FEATURE PIPELINE: INTRON & GENE MODE")
print(f"Normalizing by Genome GC Content (Column 3 of metadata)")
print(f"{'='*80}")

# --- 3. Main Loop ---
for feat_name, feat_file in FEATURES_PRIORITY:
    print(f"\n💎 FEATURE: {feat_name.upper()}")
    
    feat_dir = os.path.join(MASTER_OUTPUT, feat_name)
    os.makedirs(feat_dir, exist_ok=True)
    out_path = os.path.join(feat_dir, f"All_Vertebrates_{feat_name}_Scaled_G4.tsv")
    
    with open(out_path, 'w') as f_out:
        f_out.write("bin_num\tcount\tnormalized_density\tspecies\tclade\n")
    
    taxa_list = ["Amphibians", "Aves", "Reptiles", "Pisces", "Mammals"]
    
    for taxon in taxa_list:
        taxon_path = pathlib.Path(BASE_DIR) / taxon
        g4_path = pathlib.Path(G4_DIR) / taxon
        if not taxon_path.exists(): continue
            
        species_dirs = [d for d in taxon_path.iterdir() if d.is_dir() and not d.name.startswith('.')]
        print(f"\n📂 CLADE: {taxon}")
        
        for s_idx, s_dir in enumerate(species_dirs, 1):
            species_name = s_dir.name
            standard_name = species_name.replace(" ", "_")
            g4_file = g4_path / f"{standard_name}_G4.bed"
            coord_file = s_dir / feat_file
            
            if not g4_file.exists() or not coord_file.exists():
                print(f"  → [{s_idx}/{len(species_dirs)}] {species_name}: SKIP (Missing File)")
                continue

            try:
                total_lines = int(subprocess.check_output(f"wc -l < '{coord_file}'", shell=True))
                if total_lines == 0: continue

                bin_totals = {i: [0.0, 0] for i in range(1, 101)}
                processed_regions = 0
                
                for chunk in pd.read_csv(coord_file, sep='\t', header=None, chunksize=BATCH_SIZE, comment='#'):
                    chunk.columns = ['chrom', 'start', 'end', 'id', 'score', 'strand'] + list(chunk.columns[6:])
                    original_size = len(chunk)

                    chunk = chunk[(chunk['end'] - chunk['start']) >= 100].copy()
                    if chunk.empty: 
                        processed_regions += original_size
                        continue
                    
                    with tempfile.TemporaryDirectory() as tmpdir:
                        c_bed, c_bins, c_counts = os.path.join(tmpdir, "ch.bed"), os.path.join(tmpdir, "bi.bed"), os.path.join(tmpdir, "co.bed")
                        chunk.to_csv(c_bed, sep='\t', index=False, header=False)
                        
                        subprocess.run(f"bedtools makewindows -b '{c_bed}' -n 100 -i srcwinnum > '{c_bins}'", shell=True, capture_output=True)
                        subprocess.run(f"bedtools intersect -c -a '{c_bins}' -b '{g4_file}' > '{c_counts}'", shell=True, capture_output=True)
                        
                        if not os.path.exists(c_counts) or os.path.getsize(c_counts) == 0: 
                            processed_regions += original_size
                            continue

                        counts = pd.read_csv(c_counts, sep='\t', header=None, names=['chr', 's', 'e', 'name_bin', 'cnt'], dtype={0: str})
                        counts[['feat_id', 'bin_num']] = counts['name_bin'].str.rsplit('_', n=1, expand=True)
                        counts['bin_num'] = counts['bin_num'].astype(int)
                        
                        strand_map = chunk[['id', 'strand']].drop_duplicates()
                        counts = counts.merge(strand_map, left_on='feat_id', right_on='id')
                        
                        # 
                        counts.loc[counts['strand'] == '-', 'bin_num'] = 101 - counts['bin_num']
                        
                        batch_sum = counts.groupby('bin_num')['cnt'].sum()
                        for b_num, b_sum in batch_sum.items():
                            bin_totals[b_num][0] += b_sum
                            bin_totals[b_num][1] += len(chunk)
                        
                        processed_regions += original_size

                    percent = (processed_regions / total_lines) * 100
                    sys.stdout.write(f"\r  → {species_name}: {percent:.1f}% ({processed_regions}/{total_lines})")
                    sys.stdout.flush()

                if processed_regions > 0:
                    # --- FIX: Align species_name with metadata spaces ---
                    lookup_name = species_name.replace("_", " ")
                    species_gc_val = gc_map.get(lookup_name, 50)
                    species_gc_ratio = species_gc_val / 100.0
                    
                    final_data = []
                    for b in range(1, 101):
                        avg_cnt = bin_totals[b][0] / bin_totals[b][1] if bin_totals[b][1] > 0 else 0
                        norm_dens = avg_cnt / species_gc_ratio
                        final_data.append([b, avg_cnt, norm_dens, species_name, taxon])
                    
                    with open(out_path, 'a') as f:
                        pd.DataFrame(final_data).to_csv(f, sep='\t', header=False, index=False)
                        f.flush()
                        os.fsync(f.fileno())
                        
                    sys.stdout.write(f" ✅ DONE (Norm: {species_gc_ratio:.4f})\n")
                else:
                    sys.stdout.write(" ⚠️ SKIPPED (No valid regions)\n")

            except Exception as e:
                sys.stdout.write(f" ❌ ERROR in {species_name}: {e}\n")

print(f"\n{'='*80}\n🎉 PIPELINE FINISHED: Intron and Gene results updated.\n{'='*80}")

--------------------------------------------------------------------------------------------
#Run this for CDS (Coding DNA Sequence) - normalised by Genome GC-content
--------------------------------------------------------------------------------------------

import os
import subprocess
import pandas as pd
import pathlib
import tempfile
import sys

# --- 1. Configuration ---
BASE_DIR = "path to the directory /Genome_specific_regions/Genome_specific_coordinates"
G4_DIR = "path to the directory /Genome_specific_regions/Predicted_G4_bed_genomes"
METADATA_PATH = "path to the directory /Genome_specific_regions/Metadata/Genomes.tsv"
MASTER_OUTPUT = "path to the directory /Genome_specific_regions/Output"

# Syncing logic with your Intron feature processing
FEATURES_PRIORITY = [("CDS", "CDS_coordinates.bed")]
BATCH_SIZE = 20000 

# --- 2. Load Metadata (Standardized to Column 3) ---
meta = pd.read_csv(METADATA_PATH, sep='\t')
# Standardizing: Species is Column 0, Genome_GC_content is Column 2 (Index 2)
gc_map = dict(zip(meta.iloc[:, 0], meta.iloc[:, 2]))

print(f"\n{'='*80}")
print("🚀 G4 HEAVY FEATURE PIPELINE: CDS MODE (SYNCHRONIZED)")
print(f"Normalizing by Genome GC Content (Column 3 of metadata)")
print(f"{'='*80}")

# --- 3. Main Loop ---
for feat_name, feat_file in FEATURES_PRIORITY:
    print(f"\n💎 FEATURE: {feat_name.upper()}")
    
    feat_dir = os.path.join(MASTER_OUTPUT, feat_name)
    os.makedirs(feat_dir, exist_ok=True)
    out_path = os.path.join(feat_dir, f"All_Vertebrates_{feat_name}_Scaled_G4.tsv")
    
    # Initialize Master Output with Header
    with open(out_path, 'w') as f_out:
        f_out.write("bin_num\tcount\tnormalized_density\tspecies\tclade\n")
    
    taxa_list = ["Amphibians", "Aves", "Reptiles", "Pisces", "Mammals"]
    
    for taxon in taxa_list:
        taxon_path = pathlib.Path(BASE_DIR) / taxon
        g4_path = pathlib.Path(G4_DIR) / taxon
        if not taxon_path.exists(): continue
            
        species_dirs = [d for d in taxon_path.iterdir() if d.is_dir() and not d.name.startswith('.')]
        print(f"\n📂 CLADE: {taxon}")
        
        for s_idx, s_dir in enumerate(species_dirs, 1):
            species_name = s_dir.name
            standard_name = species_name.replace(" ", "_")
            g4_file = g4_path / f"{standard_name}_G4.bed"
            coord_file = s_dir / feat_file
            
            if not g4_file.exists() or not coord_file.exists():
                print(f"  → [{s_idx}/{len(species_dirs)}] {species_name}: SKIP (Missing File)")
                continue

            try:
                # Progress logic: count total lines in the species file
                total_lines = int(subprocess.check_output(f"wc -l < '{coord_file}'", shell=True))
                if total_lines == 0: continue

                bin_totals = {i: [0.0, 0] for i in range(1, 101)}
                processed_regions = 0
                
                # Batch Processing Loop for memory efficiency
                for chunk in pd.read_csv(coord_file, sep='\t', header=None, chunksize=BATCH_SIZE, comment='#'):
                    chunk.columns = ['chrom', 'start', 'end', 'id', 'score', 'strand'] + list(chunk.columns[6:])
                    original_size = len(chunk)
                    
                    # 1. Length Filter (100bp minimum)
                    chunk = chunk[(chunk['end'] - chunk['start']) >= 100].copy()
                    if chunk.empty: 
                        processed_regions += original_size
                        continue
                    
                    with tempfile.TemporaryDirectory() as tmpdir:
                        c_bed, c_bins, c_counts = os.path.join(tmpdir, "ch.bed"), os.path.join(tmpdir, "bi.bed"), os.path.join(tmpdir, "co.bed")
                        chunk.to_csv(c_bed, sep='\t', index=False, header=False)
                        
                        # 2. bedtools scaling and intersection
                        subprocess.run(f"bedtools makewindows -b '{c_bed}' -n 100 -i srcwinnum > '{c_bins}'", shell=True, capture_output=True)
                        subprocess.run(f"bedtools intersect -c -a '{c_bins}' -b '{g4_file}' > '{c_counts}'", shell=True, capture_output=True)
                        
                        if not os.path.exists(c_counts) or os.path.getsize(c_counts) == 0:
                            processed_regions += original_size
                            continue

                        # 3. Load Results
                        counts = pd.read_csv(c_counts, sep='\t', header=None, names=['chr', 's', 'e', 'name_bin', 'cnt'], dtype={0: str})
                        
                        # Robust splitting for ID and bin number
                        counts[['feat_id', 'bin_num']] = counts['name_bin'].str.rsplit('_', n=1, expand=True)
                        counts['bin_num'] = counts['bin_num'].astype(int)
                        
                        strand_map = chunk[['id', 'strand']].drop_duplicates()
                        counts = counts.merge(strand_map, left_on='feat_id', right_on='id')
                        
                        # 4. Flip bins for negative strand (Biological 5' to 3' alignment)
                        # Bin 1 will always represent the Start Codon
                        counts.loc[counts['strand'] == '-', 'bin_num'] = 101 - counts['bin_num']
                        
                        # Update running totals
                        batch_sum = counts.groupby('bin_num')['cnt'].sum()
                        batch_n = counts.groupby('bin_num')['feat_id'].nunique()

                        for b_num in range(1, 101):
                            if b_num in batch_sum:
                                bin_totals[b_num][0] += batch_sum[b_num]
                                bin_totals[b_num][1] += batch_n[b_num]
                        
                        processed_regions += original_size

                    # Real-time percentage update
                    percent = (processed_regions / total_lines) * 100
                    sys.stdout.write(f"\r  → {species_name}: {percent:.1f}% ({processed_regions}/{total_lines})")
                    sys.stdout.flush()

                if processed_regions > 0:
                    # --- FIX: ALIGN NAME WITH METADATA ---
                    lookup_name = species_name.replace("_", " ")
                    species_gc_val = gc_map.get(lookup_name, 50)
                    species_gc_ratio = species_gc_val / 100.0
                    
                    final_data = []
                    for b in range(1, 101):
                        avg_cnt = bin_totals[b][0] / bin_totals[b][1] if bin_totals[b][1] > 0 else 0
                        norm_dens = avg_cnt / species_gc_ratio
                        final_data.append([b, avg_cnt, norm_dens, species_name, taxon])
                    
                    with open(out_path, 'a') as f:
                        pd.DataFrame(final_data).to_csv(f, sep='\t', header=False, index=False)
                        f.flush()
                        os.fsync(f.fileno())
                        
                    sys.stdout.write(f" ✅ DONE (Norm: {species_gc_ratio:.4f})\n")
                else:
                    sys.stdout.write(" ⚠️ SKIPPED\n")

            except Exception as e:
                sys.stdout.write(f" ❌ ERROR: {e}\n")

print(f"\n{'='*80}\n🎉 PIPELINE FINISHED: CDS results synced to disk.\n{'='*80}")
