# 🧬 MotifBreakR Motif Disruption Analysis Toolkit

A lightweight and reproducible toolkit to perform motif disruption analysis using the MotifBreakR R package.
It reads a BED-style variant file, extracts alleles automatically, evaluates TF binding disruption using JASPAR2024 motifs, and exports a comprehensive results table — ready for downstream interpretation and visualization

---

## 📂 Scripts Overview

| Script | Purpose |
|:--------|:---------|
| `Input BED file` | Contains variants with genomic positions and embedded REF/ALT alleles. |
| `snps.from.file()` | Extracts SNPs, REF/ALT, and rsID using dbSNP GRCh38. |
| `motifbreakR()` | Predicts TF motif disruption across JASPAR2024 PWM library.. |
| `calculatePvalue()` | Computes p-values for disruption strength per motif. |

---

## ⚙️ Requirements
- `motifbreakR` 
- `BSgenome.Hsapiens.UCSC.hg38`  
- `MotifDb`  
- BiocParallel

---

## 🚀 Workflow

### **Step 1 — Load Required Libraries**
```R
./varID_query_gnomAD.sh --vl variants.txt --data-type <exome|genome|joint>


🔹 Example
./gnomad_query.sh --vl my_variants.txt --data-type exome


⚙️ Options
Option	Description
--vl	Variant list file (one variant per line)
--data-type	Dataset type: exome, genome, or joint
--help	Show help message



📄 Input Example
1-55516888-G-A
2-1234567-T-C
10-8965432-C-T

📁 Output

All queried variants are saved in:

varStore/<variant_id>.json


Example:

varStore/1-55516888-G-A.json
