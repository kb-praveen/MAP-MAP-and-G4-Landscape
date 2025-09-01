***
# ========***********************************************************========
# @@@@@-- Methylation Associated - Point Mutation Assessment Pipeline --@@@@@
# ========***********************************************************========
***

---

$$$$$$$$$$$$$$$$$$$--------------- MAP-MAP ---------------$$$$$$$$$$$$$$$$$$$

---

# ================
Kumar, P., Morbia, I., Satish, A. L., Datta, S., & Singh, U. (2025). CGGBP1  
from higher amniotes restricts cytosine methylation and drives a GC-bias in  
transcription factor-binding sites at repressed promoters. Transcription,  
1–36. [https://doi.org/10.1080/21541264.2025.2533598]
# ================

---

!!!!!!!!!!!!!!!!!!!!!!  HoMeCellLab - IIT GANDHINAGAR  !!!!!!!!!!!!!!!!!!!!!!

---

***
# ================
# -----------------START OF METHLATION ASSOCIATED POINT MUTATION ASSESSMENT PIPELINE-----------------
# ================
***

# ================
## STEP1: FASTQ Header Modification (Mate1/Mate2 labelling): Python script
# ================

---

# ================
## STEP2: FASTQ Quality Control and Deduplication using fastp: Bash script
# ================

---

# ================
## STEP3: Paired end alignment using Bowtie2: Bash script
# ================

---

# ================
## STEP4: Alignment Post-Processing: SAM to BAM Conversion, Sorting, Indexing, and BAM to BED Conversion: Bash script
# ================

---

# ================
## STEP5: Strand segregation (Plus and Minus) after Paired-end alignment both for Mate1 and Mate2 from raw FASTQ file: Bash script
# ================

---

# ================
## STEP6: Reverse complementing the reads falling on the minus (-) strand after paired end alignment from raw FASTQ file (Mate1/Mate2): Bash script
# ================

---

# ================
## STEP7: Merge all FASTQ mapping to PLUS (+) strand: Bash script
# ================

---

# ================
## STEP8: Single end alignment using Bowtie2: Bash script
# ================

---

# ================
## STEP9: Single end alignment Post-Processing: SAM to BAM Conversion, Sorting, Indexing, and BAM to BED Conversion: Bash script
# ================

---

# ================
## STEP10: SAM file comparison (Paired end vs single end) to ensure correct alignment: Bash script
# ================

---

# ================
## STEP11: Compare reads between two text files: Python script
# ================

---

# ================
## STEP12: Modify header of commonly occurring reads between two files: Python script
# ================

---

# ================
## STEP13: Filter correctly aligned common reads from the FASTQ file: Python script
# ================

---

# ================
## STEP14: Filtering faulty reads after single end alignment: Bash script
# ================

---

# ================
## STEP15: Adding CIGAR values to BED files: Python script
# ================

---

# ================
## STEP16: Appending sequence from FASTQ files using read IDs: Python script
# ================

---

# ================
## STEP17: Using CIGAR values considering INDELs, convert sequence appended files to single base coordinate: Python script
# ================

---

# ================
## STEP18: Script for handling duplicate entries occurred due to chunk processing: Python script
# ================

---

# ================
## STEP19: Filter coverage from 4 to 1k from these files: Bash script
# ================

---

# ================
## STEP20: Modify coverage files for further processing: Bash script
# ================

---

# ================
## STEP21: Script for converting base calls into percentage: Bash script
# ================

---

# ================
## STEP22: Script for converting base call sequences to unique strings: Python script
# ================

---

# ================
## STEP23: Classifying mutation type in individual sample files: Bash script
# ================

---

# ================
## STEP24: Mutation rate calculation in individual sample files: Bash script
# ================

***

# ==================***************************************************************==================
# ------------------END OF METHLATION ASSOCIATED POINT MUTATION ASSESSMENT PIPELINE------------------
# ==================***************************************************************==================
***
