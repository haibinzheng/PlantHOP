#!/usr/bin/env Rscript

# Expression-blind feasibility audit for the author-annotated Eutrema object.
# This script never reads assay values or computes a cell-program score.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: audit_gse268881_esa_feasibility.R INPUT_RDS OUTPUT_DIR")
}

input_rds <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- args[[2L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required <- c("SeuratObject", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) {
  stop("missing R packages: ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages(library(SeuratObject))

started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
obj <- readRDS(input_rds)

metadata <- obj[[]]
feature_assay <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
feature_ids <- rownames(obj[[feature_assay]])
if (is.null(feature_ids)) feature_ids <- character()

label_pattern <- paste(
  c("cell.?type", "annotation", "(^|_)ann($|_)", "cluster", "ident"),
  collapse = "|"
)
sample_pattern <- paste(
  c("orig.ident", "sample", "replicate", "treatment", "condition", "species", "batch"),
  collapse = "|"
)

column_inventory <- data.frame(
  column = names(metadata),
  class = vapply(metadata, function(x) paste(class(x), collapse = ";"), character(1L)),
  non_missing = vapply(metadata, function(x) sum(!is.na(x)), integer(1L)),
  unique_non_missing = vapply(metadata, function(x) length(unique(x[!is.na(x)])), integer(1L)),
  label_candidate = grepl(label_pattern, names(metadata), ignore.case = TRUE, perl = TRUE),
  sample_candidate = grepl(sample_pattern, names(metadata), ignore.case = TRUE, perl = TRUE),
  stringsAsFactors = FALSE
)
write.table(
  column_inventory,
  file.path(output_dir, "metadata_column_inventory.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

candidate_columns <- column_inventory$column[
  column_inventory$label_candidate | column_inventory$sample_candidate
]
count_tables <- lapply(candidate_columns, function(column_name) {
  values <- as.character(metadata[[column_name]])
  values[is.na(values)] <- "<NA>"
  counts <- sort(table(values), decreasing = TRUE)
  data.frame(
    column = column_name,
    value = names(counts),
    n_cells = as.integer(counts),
    stringsAsFactors = FALSE
  )
})
if (length(count_tables)) {
  write.table(
    do.call(rbind, count_tables),
    file.path(output_dir, "candidate_metadata_value_counts.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )
}

required_cross_columns <- c(
  "intspace_sample", "intspace_treatment", "intspace_replicate", "intspace_celltype"
)
if (all(required_cross_columns %in% names(metadata))) {
  cross_counts <- as.data.frame(table(
    sample = metadata$intspace_sample,
    treatment = metadata$intspace_treatment,
    replicate = metadata$intspace_replicate,
    celltype = metadata$intspace_celltype,
    useNA = "ifany"
  ), stringsAsFactors = FALSE)
  cross_counts <- cross_counts[cross_counts$Freq > 0L, , drop = FALSE]
  names(cross_counts)[names(cross_counts) == "Freq"] <- "n_cells"
  write.table(
    cross_counts,
    file.path(output_dir, "sample_treatment_replicate_celltype_counts.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "<NA>"
  )
}

writeLines(feature_ids, file.path(output_dir, "feature_ids.txt"), useBytes = TRUE)
write.table(
  utils::head(data.frame(feature_id = feature_ids, stringsAsFactors = FALSE), 200L),
  file.path(output_dir, "feature_id_examples.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

assay_inventory <- lapply(Assays(obj), function(assay_name) {
  assay_obj <- obj[[assay_name]]
  data.frame(
    assay = assay_name,
    class = paste(class(assay_obj), collapse = ";"),
    n_features = nrow(assay_obj),
    n_cells = ncol(assay_obj),
    feature_metadata_columns = paste(names(assay_obj[[]]), collapse = ";"),
    stringsAsFactors = FALSE
  )
})
assay_inventory <- do.call(rbind, assay_inventory)
write.table(
  assay_inventory,
  file.path(output_dir, "assay_inventory.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

compact_metadata <- metadata[, unique(candidate_columns), drop = FALSE]
saveRDS(compact_metadata, file.path(output_dir, "author_metadata_candidates.rds"), compress = "xz")

audit <- list(
  status = "expression_blind_metadata_audit_complete",
  started = started,
  completed = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  input_rds = input_rds,
  input_bytes = unname(file.info(input_rds)$size),
  object_class = class(obj),
  cells = ncol(obj),
  features = nrow(obj),
  assays = Assays(obj),
  default_assay = DefaultAssay(obj),
  metadata_columns = ncol(metadata),
  candidate_metadata_columns = candidate_columns,
  feature_assay = feature_assay,
  feature_id_first = if (length(feature_ids)) feature_ids[[1L]] else NA_character_,
  feature_id_last = if (length(feature_ids)) feature_ids[[length(feature_ids)]] else NA_character_,
  expression_values_read = FALSE,
  program_scores_computed = FALSE,
  auroc_computed = FALSE,
  model_training_started = FALSE
)
jsonlite::write_json(
  audit,
  file.path(output_dir, "audit.json"),
  auto_unbox = TRUE, pretty = TRUE, na = "null"
)

writeLines("audit_complete", file.path(output_dir, "STATUS"))
