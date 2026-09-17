#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: audit_gse232863_structure.R INPUT_RDS OUTPUT_DIR")
}

input_rds <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

started_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
obj <- readRDS(input_rds)
metadata <- obj[[]]
feature_ids <- rownames(obj)

safe_unique_count <- function(x) {
  length(unique(x[!is.na(x)]))
}

safe_examples <- function(x, n = 8L) {
  values <- unique(as.character(x[!is.na(x)]))
  paste(head(values, n), collapse = " | ")
}

column_summary <- data.table(
  column = colnames(metadata),
  storage_class = vapply(metadata, function(x) paste(class(x), collapse = ";"), character(1)),
  non_missing = vapply(metadata, function(x) sum(!is.na(x)), integer(1)),
  missing = vapply(metadata, function(x) sum(is.na(x)), integer(1)),
  unique_non_missing = vapply(metadata, safe_unique_count, integer(1)),
  examples = vapply(metadata, safe_examples, character(1))
)
fwrite(column_summary, file.path(output_dir, "metadata_columns.tsv"), sep = "\t")

candidate_pattern <- "organ|tissue|cell.?type|celltype|annotation|cluster|ident|sample|library|batch|stage|genotype|treatment|replicate"
candidate_columns <- column_summary[
  grepl(candidate_pattern, column, ignore.case = TRUE) & unique_non_missing <= 500L,
  column
]

count_tables <- lapply(candidate_columns, function(column_name) {
  values <- as.character(metadata[[column_name]])
  values[is.na(values)] <- "<NA>"
  result <- as.data.table(table(values, useNA = "ifany"))
  setnames(result, c("value", "cell_count"))
  result[, metadata_column := column_name]
  setcolorder(result, c("metadata_column", "value", "cell_count"))
  result[order(-cell_count, value)]
})
if (length(count_tables) > 0L) {
  fwrite(rbindlist(count_tables), file.path(output_dir, "categorical_counts.tsv"), sep = "\t")
} else {
  fwrite(data.table(metadata_column = character(), value = character(), cell_count = integer()),
         file.path(output_dir, "categorical_counts.tsv"), sep = "\t")
}

if (all(c("sample", "tissue", "cluster_names") %in% colnames(metadata))) {
  sample_cluster_counts <- as.data.table(metadata)[, .(cell_count = .N), by = .(sample, tissue, cluster_names)]
  setorder(sample_cluster_counts, tissue, sample, -cell_count, cluster_names)
  fwrite(sample_cluster_counts, file.path(output_dir, "sample_cluster_counts.tsv"), sep = "\t")
}

feature_patterns <- data.table(
  pattern_name = c(
    "MSU7_LOC_Os",
    "RAPDB_Os_gene",
    "Ensembl_rice_gene",
    "Arabidopsis_AGI",
    "contains_version_suffix"
  ),
  regex = c(
    "^LOC_Os[0-9]{2}g[0-9]+",
    "^Os[0-9]{2}g[0-9]+",
    "^Osativa",
    "^AT[1-5CM]G[0-9]+",
    "\\.[0-9]+$"
  )
)
feature_patterns[, matching_features := vapply(regex, function(pattern) {
  sum(grepl(pattern, feature_ids, perl = TRUE))
}, integer(1))]
fwrite(feature_patterns, file.path(output_dir, "feature_id_patterns.tsv"), sep = "\t")
fwrite(data.table(feature_id = feature_ids), file.path(output_dir, "feature_ids.tsv"), sep = "\t")

assay_names <- names(obj@assays)
assay_summary <- rbindlist(lapply(assay_names, function(assay_name) {
  assay <- obj@assays[[assay_name]]
  data.table(
    assay = assay_name,
    assay_class = paste(class(assay), collapse = ";"),
    note = "Dimensions and layers intentionally not accessed during pre-freeze structure audit"
  )
}))
fwrite(assay_summary, file.path(output_dir, "assays.tsv"), sep = "\t")

summary_record <- list(
  started_utc = started_utc,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  input_path = input_rds,
  input_bytes = unname(file.info(input_rds)$size),
  object_class = class(obj),
  object_size_bytes_in_memory = as.numeric(object.size(obj)),
  cells = nrow(metadata),
  features = length(feature_ids),
  active_assay = obj@active.assay,
  assays = assay_names,
  reductions = names(obj@reductions),
  graphs = names(obj@graphs),
  metadata_columns = colnames(metadata),
  candidate_categorical_columns = candidate_columns,
  feature_id_examples = head(feature_ids, 20L),
  r_version = R.version.string,
  seurat_version = as.character(packageVersion("Seurat")),
  seuratobject_version = as.character(packageVersion("SeuratObject")),
  expression_values_inspected = FALSE,
  candidate_programs_loaded = FALSE,
  source_rds_modified = FALSE
)
write_json(summary_record, file.path(output_dir, "structure_summary.json"), pretty = TRUE, auto_unbox = TRUE)

rm(obj, metadata)
invisible(gc())
