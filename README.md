🧬 MotifBreakR Analysis Pipeline

This repository provides a complete and reproducible workflow for performing motif disruption analysis using the MotifBreakR package in R.
It explains the input format, required dependencies, analysis pipeline, visualization, and exporting results.

📥 Input Format

The analysis requires a BED-like input file with four columns:

Column	Description
CHR	Chromosome
START	1-based start genomic coordinate of variant
END	End coordinate of variant
ID	Variant ID in the format chr:position:REF:ALT
🔹 Example
chr13	51936209	51936210	chr13:51936210:G:A
chr13	51938272	51938273	chr13:51938273:C:A
chr13	51943391	51943392	chr13:51943392:A:G
chr13	51943461	51943462	chr13:51943462:G:T
...


⚠ REF and ALT alleles must be embedded in the fourth column (ID) for MotifBreakR to interpret them correctly.

📂 Working Directory
setwd("/home/srishti-sharma/Desktop/Srishti/motif/final_analysis/")

🔧 Load Required Libraries
library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(MotifDb)
library(BiocParallel)
library(BSgenome)

🧮 Extract Variants from BED
bed_file <- "Motif_Breaker - SKT.bed"

snps <- snps.from.file(
  file = bed_file,
  format = "bed",
  dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
  search.genome = BSgenome.Hsapiens.UCSC.hg38,
  check.unnamed.for.rsid = TRUE
)


MotifBreakR automatically retrieves:

Chromosomal coordinates

REF and ALT alleles

dbSNP rsID (if available)

🧬 Load Motifs (JASPAR 2024 – Human)
human.jaspar2024 <- query(MotifDb, c("jaspar2024", "Hsapiens"))

🚀 Run MotifBreakR
results <- motifbreakR(
  snpList = snps,
  filterp = TRUE,
  pwmList = human.jaspar2024,
  threshold = 1e-4,
  method = "log",
  bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
  BPPARAM = BiocParallel::SerialParam()
)

📊 Convert Results to DataFrame
df_results <- as.data.frame(results, row.names = NULL)

df_results[] <- lapply(df_results, function(col) {
  if (is.list(col)) sapply(col, function(x) paste(x, collapse = ";")) else col
})

💾 Export Output
write.table(
  df_results,
  "motifbreakr_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

🔍 Examine / Visualize a Single Variant (Optional)
rsVar <- results[names(results) %in% "rs1479475149"]
rsVar <- calculatePvalue(rsVar)
plotMB(results = results, rsid = "rs1479475149", effect = "strong")
