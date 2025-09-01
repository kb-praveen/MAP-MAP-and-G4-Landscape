# 🧬 Vertebrate Genomes G4-Potential Analysis 🌀  
### & GC-Retention Analysis Using MAP-MAP

---

## 🔬 Pipeline Overview:  
This pipeline analyzes the **potential G4 regions** in vertebrate genomes and their promoters across taxa—including five major vertebrate classes:  
🐟 Fish | 🐸 Amphibia | 🦎 Reptilia | 🐦 Aves | 🦁 Mammalia  

It is coupled with GC-retention analysis, employing the MAP-MAP framework to uncover methylation associated point mutation in potential G4-forming regions in Human genome and influence of various forms of CGGBP1.  

---

## ⚙️ Stepwise Pipeline Description

**# =============================================================**  
**# 🚩 STEP01 - G4 PREDICTION 🌀**  
🔍 Predict G4-forming regions in vertebrate genomes and their promoters across taxa using `pqsfinder`. Also generates summary of associated features.  
*Custom R script*  

**# =============================================================**  

- **🗂️ STEP02:** Extract bed coordinates of potential G4 regions.  
  *Custom R script*  

- **🔄 STEP03:** Single base bed file conversion of G4 coordinates.  
  *Custom Bash script*  

- **⚙️ STEP04:** Conversion of bed files to extract mutation profiles from MeDIP data.  
  *Custom Bash script*  

- **🔎 STEP05:** Grep mutation profile from MeDIP data.  
  *Custom Bash script*  

- **📊 STEP06:** Count mutation types in G4-containing region for each sample.  
  *Custom R script*  

- **🧬 STEP07:** Extract C to T mutations from MeDIP data for each sample.  
  *Custom R script*  

- **⚡ STEP08:** Extract Transition and Transversion mutations from raw files.  
  *Custom Bash script*  

- **🎯 STEP09:** Extract Transition and Transversion mutations ±1 kb upstream and downstream from sequence files.  
  *Custom Bash script*  

- **📈 STEP10:** Generate hg38 G4 bigwig files from potential G4-forming regions predicted by `pqsfinder`.  
  *Custom Bash script*  

- **🗃️ STEP11:** Convert G4 coordinate beds for C to T Transition, Transversion, and mutation-free regions.  
  *Custom Bash script*  

- **🌐 STEP12:** Generate scaled bigwig genome-wide for G4 and GC skew signals.  
  *Custom Bash script*  

- **🧮 STEP13:** Extract all transition mutations from raw files and generate corresponding bed files.  
  *Custom R script*  

- **🎨 STEP14:** Perform mutation segregation and generate bed files for promoter G4-forming regions.  
  *Custom Bash script*  

---

## 📌 Notes

- All scripts follow a modular approach for independent execution and easy integration.  
- The pipeline blends R and Bash scripts tailored for genomic data processing.  
- Visualization and summary steps complement raw data extraction to provide in-depth G4 landscape insights.

---

*Explore the pipeline scripts folder for detailed documentation and usage examples for each step!* 🚀
