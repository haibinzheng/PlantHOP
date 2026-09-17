#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: inspect_gse232863_rna_assay.R INPUT_RDS OUTPUT_JSON")
}

input_rds <- normalizePath(args[[1]], mustWork = TRUE)
output_json <- args[[2]]
obj <- readRDS(input_rds)
rna <- obj@assays[["RNA"]]

matrix_record <- function(name, matrix) {
  values <- matrix@x
  list(
    name = name,
    class = class(matrix),
    rows = nrow(matrix),
    columns = ncol(matrix),
    nonzero_entries = length(values),
    finite_values = sum(is.finite(values)),
    minimum_stored_value = if (length(values)) min(values, na.rm = TRUE) else NA_real_,
    maximum_stored_value = if (length(values)) max(values, na.rm = TRUE) else NA_real_,
    nonnegative = if (length(values)) all(values >= 0, na.rm = TRUE) else TRUE,
    integer_valued = if (length(values)) all(abs(values - round(values)) < 1e-8, na.rm = TRUE) else TRUE
  )
}

record <- list(
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  source_rds = input_rds,
  source_rds_modified = FALSE,
  assay = "RNA",
  assay_class = class(rna),
  assay_slots = slotNames(rna),
  counts = matrix_record("counts", rna@counts),
  data = matrix_record("data", rna@data),
  scale_data_dimensions = dim(rna@scale.data),
  feature_count = nrow(rna),
  cell_count = ncol(rna),
  candidate_programs_loaded = FALSE,
  program_scores_computed = FALSE
)
write_json(record, output_json, pretty = TRUE, auto_unbox = TRUE)

rm(obj, rna)
invisible(gc())
