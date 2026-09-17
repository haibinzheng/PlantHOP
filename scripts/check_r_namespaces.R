#!/usr/bin/env Rscript

packages <- c(
  "Seurat", "SeuratObject", "Signac", "BiocGenerics", "S4Vectors",
  "IRanges", "GenomeInfoDbData", "GenomeInfoDb", "XVector",
  "Biostrings", "GenomicRanges", "Rsamtools"
)

for (package in packages) {
  available <- requireNamespace(package, quietly = TRUE)
  version <- if (available) as.character(packageVersion(package)) else ""
  cat(package, available, version, sep = "\t")
  cat("\n")
}
