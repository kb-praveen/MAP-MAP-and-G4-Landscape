
# 🧬 Methylation-Associated Point-Mutation Assessment Pipeline 🧬
# ( MAP-MAP )
**Authors: Praveen Kumar and Umashankar Singh**
**HoMeCellLab - IIT GANDHINAGAR**

*Kumar, P., Morbia, I., Satish, A. L., Datta, S., & Singh, U. (2025). CGGBP1 from higher amniotes restricts cytosine methylation and drives a GC-bias in transcription factor-binding sites at repressed promoters. Transcription, 1–36. [https://doi.org/10.1080/21541264.2025.2533598]*
---
## 📊 Stepwise visualization of the pipeline flowchart

**A 🚀 STEP1:** FASTQ Header Modification: (Mate1/Mate2 labelling) - Python script  
**B 🧹 STEP2:** FASTQ Quality Control and Deduplication using fastp - Bash script  
**C 🎯 STEP3:** Paired end alignment using Bowtie2 - Bash script  
**D 🛠 STEP4:** Alignment Post-Processing: SAM to BAM Conversion, Sorting, Indexing, and BAM to BED Conversion - Bash script  
**E 🔀 STEP5:** Strand segregation (Plus and Minus) after Paired-end alignment both for Mate1 and Mate2 from raw FASTQ file - Bash script  
**F 🔄 STEP6:** Reverse complementing the reads falling on minus (-) strand after paired end alignment from raw FASTQ file (Mate1/Mate2) - Bash script  
**G 🧩 STEP7:** Merge all FASTQ mapping to PLUS (+) strand - Bash script  
**H 🎯 STEP8:** Single end alignment using Bowtie2 - Bash script  
**I 🛠 STEP9:** Single end alignment Post-Processing: SAM to BAM Conversion, Sorting, Indexing, and BAM to BED Conversion - Bash script  
**J 🔎 STEP10:** SAM file comparison (Paired end vs single end) to make sure alignment occurred correctly - Bash script  
**K 📊 STEP11:** Compare reads between two text files - Python script  
**L ✍️ STEP12:** Modify header of commonly occurring reads between two files - Python script  
**M ✅ STEP13:** Filter correctly aligned common reads from the FASTQ file - Python script  
**N 🚫 STEP14:** Filtering faulty reads after single end alignment - Bash script  
**O 🧬 STEP15:** Adding CIGAR values to bed files - Python script  
**P ➕ STEP16:** Appending sequence from fastq files using read IDs - Python script  
**Q 🔧 STEP17:** Using CIGAR values considering INDELs convert sequence appended files to single base coordinate - Python script  
**R ♻️ STEP18:** Script for handling Duplicate entries occurred due to chunk process - Python script  
**S 📉 STEP19:** Filter coverage from 4 to 1k from these files - Bash script  
**T 📝 STEP20:** Modify coverage files for further processing - Bash script  
**U 📈 STEP21:** Script for converting Base calls into percentage - Bash script  
**V 🔡 STEP22:** Script for Converting Base call Sequences to unique strings - Python script  
**W ⚙️ STEP23:** Classifying mutation type in individual sample files - Bash script  
**X 📊 STEP24:** Mutation rate calculation in individual sample files - Bash script  

A -> B -> C -> D -> E -> F -> G -> H -> I -> J -> K -> L -> M -> N -> O -> P -> Q -> R -> S -> T -> U -> V -> W -> X

---

## 📌 Notes

- 🐍 Python scripts are highlighted with 🧬 and 🐍 emojis.
- 🐚 Bash scripts use tools like 🧹 (cleanup), 🎯 (alignment), 🔄 (processing), etc.
- The color shading improves readability and distinction among steps.
- This flowchart is rendered automatically on GitHub when Mermaid support is enabled.
- Emojis next to step titles make the flow visually engaging and easier to identify step types.

---

*End of Methylation-Associated Point-Mutation Analysis Pipeline 🚀*

---


