📘 MotifBreakR Analysis — Documentation & Pipeline

This repository provides a complete workflow for performing motif disruption analysis using the MotifBreakR package in R.

📥 1. Example Input Format

The input file contains four columns:

Column	Description
CHR	Chromosome
START	1-based start genomic coordinate
END	End position of the variant
ID	Variant identifier in the format chr:position:REF:ALT
🔹 Example Input (BED format)
chr13	51936209	51936210	chr13:51936210:G:A
chr13	51938272	51938273	chr13:51938273:C:A
chr13	51943391	51943392	chr13:51943392:A:G
chr13	51943461	51943462	chr13:51943462:G:T
chr13	51944785	51944786	chr13:51944786:A:T
chr13	51945935	51945936	chr13:51945936:C:T
chr13	51959329	51959330	chr13:51959330:G:C
chr13	51960074	51960075	chr13:51960075:G:C
chr13	52011454	52011455	chr13:52011455:T:CGGCG
...

📂 2. Working Directory Setup
setwd("/home/srishti-sharma/Desktop/Srishti/motif/final_analysis/")

🔧 3. Load Required Libraries
library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(MotifDb)
library(BiocParallel)
library(BSgenome)

📤 4. Read Input BED File & Extract SNPs
bed_file <- "Motif_Breaker - SKT.bed"

snps <- snps.from.file(
  file = bed_file,
  format = "bed",
  dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
  search.genome = BSgenome.Hsapiens.UCSC.hg38,
  check.unnamed.for.rsid = TRUE
)


MotifBreakR automatically extracts:
✔ Chromosomal positions
✔ REF and ALT alleles
✔ rsID (if matched in dbSNP)
✔ Mappings in genome database

🧬 5. Load JASPAR2024 Human Motifs
human.jaspar2024 <- query(MotifDb, c("jaspar2024", "Hsapiens"))

🚀 6. Run MotifBreakR
results <- motifbreakR(
  snpList = snps,
  filterp = TRUE,
  pwmList = human.jaspar2024,
  threshold = 1e-4,
  method = "log",
  bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
  BPPARAM = BiocParallel::SerialParam()
)

📊 7. Convert Results to DataFrame
df_results <- as.data.frame(results, row.names = NULL)

df_results[] <- lapply(df_results, function(col) {
  if (is.list(col)) sapply(col, function(x) paste(x, collapse = ";")) else col
})

💾 8. Save Results
write.table(
  df_results,
  "motifbreakr_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

🔍 Optional: Examine a Single Variant
rs1479475149 <- results[names(results) %in% "rs1479475149"]
rs1479475149 <- calculatePvalue(rs1479475149)

📈 Plot Motif Disruption Effect
plotMB(results = results, rsid = "rs1479475149", effect = "strong")
