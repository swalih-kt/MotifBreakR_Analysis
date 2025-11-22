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



### **Steps**

library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(MotifDb)
library(BiocParallel)
library(BSgenome)

### **Step 1 — Load Required Libraries**
**Step 3 — Read Variants from BED File**

snps <- snps.from.file(
  file = "Motif_Breaker - SKT.bed",
  format = "bed",
  dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
  search.genome = BSgenome.Hsapiens.UCSC.hg38,
  check.unnamed.for.rsid = TRUE
)

**Step 4 — Load Motif PWMs (JASPAR 2024)**
human.jaspar2024 <- query(MotifDb, c("jaspar2024", "Hsapiens"))

Step 5 — Run MotifBreakR
results <- motifbreakR(
  snpList = snps,
  filterp = TRUE,
  pwmList = human.jaspar2024,
  threshold = 1e-4,
  method = "log",
  bkg = c(A=.25, C=.25, G=.25, T=.25),
  BPPARAM = BiocParallel::SerialParam()
)

**Step 6 — Export Results to TSV**
df_results <- as.data.frame(results, row.names = NULL)
df_results[] <- lapply(df_results, function(col)
  if (is.list(col)) sapply(col, paste, collapse = ";") else col)
write.table(df_results, "motifbreakr_results.tsv", sep="\t",
            quote = FALSE, row.names = FALSE)

🔍 Optional: Examine / Visualize a Single Variant
rsV <- results[names(results) %in% "rs1479475149"]
rsV <- calculatePvalue(rsV)
plotMB(results = results, rsid = "rs1479475149", effect = "strong")
