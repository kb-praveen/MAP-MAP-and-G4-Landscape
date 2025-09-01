#!/bin/bash

# Input and output paths
IN_BEDGRAPH="Path_to_begraph/hg38_windows_10bp_g4_signal.bedGraph"
CHROMSIZES="Path_to_chr_file/hg38.all_chrom.sizes1"
OUT_BASE="Path_to_scaled_bigwigs/G4_bigwigs/"

declare -a profiles=("No_mutations" "Transversion" "Transition" "C_to_T_transition")
declare -a samples=("Lc" "Ac" "Gg" "Hs" "Ev" "Input")

# Define scaling factors as an associative array
declare -A scale
scale["No_mutations,Lc"]=0.829777465117616
scale["No_mutations,Ac"]=1
scale["No_mutations,Gg"]=0.588400382898829
scale["No_mutations,Hs"]=0.640173147620756
scale["No_mutations,Ev"]=0.612092441360728
scale["No_mutations,Input"]=0.339737093880258

scale["Transversion,Lc"]=0.938569457787444
scale["Transversion,Ac"]=1
scale["Transversion,Gg"]=0.587556265026758
scale["Transversion,Hs"]=0.542329495455455
scale["Transversion,Ev"]=0.650111789914111
scale["Transversion,Input"]=0.557243091932816

scale["Transition,Lc"]=0.784941708928731
scale["Transition,Ac"]=1
scale["Transition,Gg"]=0.561628672979184
scale["Transition,Hs"]=0.467064007106275
scale["Transition,Ev"]=0.571687485322135
scale["Transition,Input"]=0.541182797508809

scale["C_to_T_transition,Lc"]=0.760204783734271
scale["C_to_T_transition,Ac"]=1
scale["C_to_T_transition,Gg"]=0.614313353622038
scale["C_to_T_transition,Hs"]=0.613136141775366
scale["C_to_T_transition,Ev"]=0.661736052751268
scale["C_to_T_transition,Input"]=0.626236898156554

for profile in "${profiles[@]}"; do
    outdir="${OUT_BASE}/${profile}"
    mkdir -p "$outdir"
    for sample in "${samples[@]}"; do
        factor="${scale[${profile},${sample}]}"
        outbg="${outdir}/hg38_G4_scaled_${sample}.bedGraph"
        outbw="${outdir}/hg38_G4_scaled_${sample}.bw"
        awk -v s="$factor" 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\t%.6f\n", $1, $2, $3, $4*s}' "$IN_BEDGRAPH" > "$outbg"
        bedGraphToBigWig "$outbg" "$CHROMSIZES" "$outbw"
        echo "Saved: $outbg and $outbw (scaled by $factor)"
    done
done

echo "All scaled bedGraph and bigWig files generated and organized."
