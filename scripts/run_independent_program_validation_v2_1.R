#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

project <- "."
meta_dir <- file.path(project, "metadata")
report_dir <- file.path(project, "reports")
source_rds <- "data/external/GSE232863/GSE232863_scRNA_omics.Rds"
orthogroups_path <- "data/runs/PhyloOpenCell/orthofinder_full7_v1/results/Results_Sep10/Orthogroups/Orthogroups.tsv"
run_dir <- "data/runs/PhyloOpenCell/independent_program_validation_v2_1"
partial_dir <- paste0(run_dir, ".partial")
derived_dir <- "data/derived/PhyloOpenCell/external_validation_v2"
top_k <- 128L
random_sets <- 500L
seed <- 20260914L

if (dir.exists(run_dir) || dir.exists(partial_dir)) {
  stop("Refusing to overwrite existing v2.1 run or partial directory")
}
dir.create(partial_dir, recursive = TRUE)
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

freeze <- fread(file.path(meta_dir, "external_validation_v2_freeze_v1.tsv"))
candidates <- fread(file.path(meta_dir, "conserved_program_candidates_v2.tsv"))
feature_map <- fread(file.path(meta_dir, "external_validation_v2_rice_feature_orthogroup_map_v1.tsv"))
addendum_path <- file.path(report_dir, "external_validation_v2_execution_addendum_v1.json")
if (!file.exists(addendum_path)) stop("Execution addendum must exist before scoring")

off_target_keys <- c(
  rice_leaf_photosynthetic_v1 = "leaf:epidermal_system",
  rice_leaf_epidermal_v1 = "leaf:photosynthetic_ground_tissue",
  rice_root_epidermal_v1 = "root:stele_lineage",
  rice_root_ground_v1 = "root:epidermal_system",
  rice_root_stele_v1 = "root:epidermal_system",
  rice_root_cap_v1 = "root:ground_tissue"
)

split_semicolon <- function(value) {
  parts <- strsplit(value, ";", fixed = TRUE)[[1]]
  parts[nzchar(parts)]
}

roc_auc <- function(y, score) {
  n_pos <- sum(y == 1L)
  n_neg <- sum(y == 0L)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

average_precision <- function(y, score) {
  if (sum(y == 1L) == 0L) return(NA_real_)
  ordering <- order(-score, seq_along(score))
  sorted_y <- y[ordering]
  precision <- cumsum(sorted_y == 1L) / seq_along(sorted_y)
  mean(precision[sorted_y == 1L])
}

standardized_effect <- function(y, score) {
  positive <- score[y == 1L]
  negative <- score[y == 0L]
  pooled <- sqrt((var(positive) + var(negative)) / 2)
  if (!is.finite(pooled) || pooled == 0) return(0)
  (mean(positive) - mean(negative)) / pooled
}

metric_record <- function(y, score) {
  list(
    roc_auc = roc_auc(y, score),
    average_precision = average_precision(y, score),
    standardized_mean_difference = standardized_effect(y, score),
    mean_positive = mean(score[y == 1L]),
    mean_negative = mean(score[y == 0L])
  )
}

gene_count <- function(value) {
  if (is.na(value) || !nzchar(trimws(value))) return(0L)
  length(strsplit(value, ",", fixed = TRUE)[[1]])
}

started_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
obj <- readRDS(source_rds)
rna <- obj@assays[["RNA"]]
counts <- rna@counts
cell_meta <- as.data.table(obj@meta.data)
cell_meta[, source_cell_index := .I]

map_lookup <- setNames(feature_map$orthogroup_id, feature_map$source_feature_id)
feature_ogs <- unname(map_lookup[rownames(counts)])
feature_ogs[is.na(feature_ogs) | feature_ogs == ""] <- NA_character_
og_levels <- sort(unique(feature_ogs[!is.na(feature_ogs)]))
og_to_index <- setNames(seq_along(og_levels), og_levels)
feature_og_index <- unname(og_to_index[feature_ogs])

all_samples <- unique(unlist(lapply(freeze$samples, split_semicolon)))
selected_source_indices <- which(cell_meta$sample %in% all_samples)
selected_meta <- cell_meta[selected_source_indices, .(source_cell_index, sample, tissue, cluster_names)]

j_parts <- vector("list", length(selected_source_indices))
x_parts <- vector("list", length(selected_source_indices))
for (row_index in seq_along(selected_source_indices)) {
  source_col <- selected_source_indices[[row_index]]
  start <- counts@p[[source_col]] + 1L
  end <- counts@p[[source_col + 1L]]
  if (start > end) next
  feature_indices <- counts@i[start:end] + 1L
  values <- counts@x[start:end]
  mapped_indices <- feature_og_index[feature_indices]
  keep <- !is.na(mapped_indices) & is.finite(values) & values > 0
  if (!any(keep)) next
  feature_indices <- feature_indices[keep]
  mapped_indices <- mapped_indices[keep]
  values <- values[keep]
  ordering <- order(-values, feature_indices)
  ordering <- head(ordering, top_k)
  mapped_indices <- mapped_indices[ordering]
  weights <- 1 / log2(seq_along(ordering) + 1)
  best <- tapply(weights, mapped_indices, max)
  norm_value <- sqrt(sum(best ^ 2))
  if (!is.finite(norm_value) || norm_value == 0) next
  j_parts[[row_index]] <- as.integer(names(best))
  x_parts[[row_index]] <- as.numeric(best) / norm_value
}

entry_counts <- lengths(j_parts)
encoded <- sparseMatrix(
  i = rep.int(seq_along(selected_source_indices), entry_counts),
  j = unlist(j_parts, use.names = FALSE),
  x = unlist(x_parts, use.names = FALSE),
  dims = c(length(selected_source_indices), length(og_levels)),
  dimnames = list(NULL, og_levels)
)

rm(j_parts, x_parts)
invisible(gc())

orthogroups <- fread(orthogroups_path)
setkey(orthogroups, Orthogroup)
species_columns <- setdiff(names(orthogroups), "Orthogroup")
species_names <- gsub("_", " ", sub("\\.representative_proteins$", "", species_columns))
species_names <- paste0(toupper(substr(species_names, 1, 1)), substr(species_names, 2, nchar(species_names)))
names(species_columns) <- species_names

candidate_sets <- split(candidates$orthogroup_id, paste(candidates$panel, candidates$coarse_label, sep = ":"))
candidate_sets <- lapply(candidate_sets, unique)

strict_single_copy <- function(target_rows) {
  result <- character()
  for (row_index in seq_len(nrow(target_rows))) {
    orthogroup <- target_rows$orthogroup_id[[row_index]]
    supporting_species <- trimws(split_semicolon(target_rows$species[[row_index]]))
    required_species <- unique(c("Oryza sativa", supporting_species))
    og_row <- orthogroups[J(orthogroup)]
    if (nrow(og_row) != 1L || !all(required_species %in% names(species_columns))) next
    counts_by_species <- vapply(required_species, function(species) {
      gene_count(og_row[[species_columns[[species]]]][[1]])
    }, integer(1))
    if (all(counts_by_species == 1L)) result <- c(result, orthogroup)
  }
  unique(result)
}

score_program <- function(matrix, orthogroups) {
  columns <- match(orthogroups, og_levels)
  columns <- columns[!is.na(columns)]
  if (!length(columns)) return(rep(0, nrow(matrix)))
  Matrix::rowSums(matrix[, columns, drop = FALSE])
}

set.seed(seed)
pooled_rows <- list()
sample_rows <- list()
summary_rows <- list()

for (contrast_index in seq_len(nrow(freeze))) {
  frozen <- freeze[contrast_index]
  contrast_id <- frozen$contrast_id[[1]]
  panel_value <- frozen$panel[[1]]
  coarse_value <- frozen$coarse_label[[1]]
  tissues <- split_semicolon(frozen$tissues[[1]])
  samples <- split_semicolon(frozen$samples[[1]])
  positives <- split_semicolon(frozen$positive_labels[[1]])
  negatives <- split_semicolon(frozen$negative_labels[[1]])
  target_key <- paste(panel_value, coarse_value, sep = ":")
  off_key <- off_target_keys[[contrast_id]]
  target_ogs <- candidate_sets[[target_key]]
  off_target_ogs <- candidate_sets[[off_key]]
  target_rows <- candidates[candidates$panel == panel_value & candidates$coarse_label == coarse_value]

  rice_column <- species_columns[["Oryza sativa"]]
  rice_counts <- vapply(target_ogs, function(orthogroup) {
    og_row <- orthogroups[J(orthogroup)]
    if (nrow(og_row) != 1L) return(0L)
    gene_count(og_row[[rice_column]][[1]])
  }, integer(1))
  single_ogs <- target_ogs[rice_counts == 1L]
  expanded_ogs <- target_ogs[rice_counts > 1L]
  strict_ogs <- strict_single_copy(target_rows)

  eligible_rows <- which(
    selected_meta$sample %in% samples &
      selected_meta$tissue %in% tissues &
      selected_meta$cluster_names %in% c(positives, negatives)
  )
  contrast_matrix <- encoded[eligible_rows, , drop = FALSE]
  contrast_meta <- selected_meta[eligible_rows]
  y <- as.integer(contrast_meta$cluster_names %in% positives)

  methods <- list(
    complete_orthogroup_family_program = target_ogs,
    target_species_single_gene_orthogroups = single_ogs,
    target_species_expanded_family_orthogroups = expanded_ogs,
    strict_single_copy_across_supporting_species = strict_ogs,
    off_target_program_negative_control = off_target_ogs
  )
  method_scores <- list()
  for (method_name in names(methods)) {
    method_ogs <- methods[[method_name]]
    scores <- score_program(contrast_matrix, method_ogs)
    method_scores[[method_name]] <- scores
    metrics <- metric_record(y, scores)
    pooled_rows[[length(pooled_rows) + 1L]] <- data.table(
      contrast_id = contrast_id,
      panel = panel_value,
      coarse_label = coarse_value,
      role = frozen$role[[1]],
      method = method_name,
      candidate_orthogroups = length(method_ogs),
      mapped_candidate_orthogroups = sum(method_ogs %in% og_levels),
      positive_cells = sum(y == 1L),
      negative_cells = sum(y == 0L),
      roc_auc = metrics$roc_auc,
      average_precision = metrics$average_precision,
      standardized_mean_difference = metrics$standardized_mean_difference,
      mean_positive = metrics$mean_positive,
      mean_negative = metrics$mean_negative
    )
    for (sample_name in samples) {
      sample_mask <- contrast_meta$sample == sample_name
      sample_y <- y[sample_mask]
      sample_scores <- scores[sample_mask]
      sample_metrics <- metric_record(sample_y, sample_scores)
      sample_rows[[length(sample_rows) + 1L]] <- data.table(
        contrast_id = contrast_id,
        sample = sample_name,
        method = method_name,
        positive_cells = sum(sample_y == 1L),
        negative_cells = sum(sample_y == 0L),
        roc_auc = sample_metrics$roc_auc,
        average_precision = sample_metrics$average_precision,
        standardized_mean_difference = sample_metrics$standardized_mean_difference
      )
    }
  }

  excluded <- unique(c(target_ogs, off_target_ogs))
  random_pool <- setdiff(og_levels, excluded)
  random_programs <- replicate(
    random_sets,
    sample(random_pool, min(length(target_ogs), length(random_pool)), replace = FALSE),
    simplify = FALSE
  )
  random_i <- unlist(lapply(seq_along(random_programs), function(index) {
    rep.int(index, length(random_programs[[index]]))
  }), use.names = FALSE)
  random_j <- match(unlist(random_programs, use.names = FALSE), og_levels)
  indicator <- sparseMatrix(
    i = random_j,
    j = random_i,
    x = 1,
    dims = c(length(og_levels), random_sets)
  )
  random_scores <- as.matrix(contrast_matrix %*% indicator)
  random_aucs <- vapply(seq_len(random_sets), function(index) roc_auc(y, random_scores[, index]), numeric(1))
  complete_scores <- method_scores[["complete_orthogroup_family_program"]]
  complete_metrics <- metric_record(y, complete_scores)
  accumulated_sample_rows <- rbindlist(sample_rows)
  current_contrast_id <- contrast_id
  complete_sample_rows <- accumulated_sample_rows[
    accumulated_sample_rows[["contrast_id"]] == current_contrast_id &
      accumulated_sample_rows[["method"]] == "complete_orthogroup_family_program"
  ]
  required_positive_samples <- if (panel_value == "leaf") 3L else 2L
  positive_sample_directions <- sum(complete_sample_rows$roc_auc > 0.5, na.rm = TRUE)
  empirical_p <- (1 + sum(random_aucs >= complete_metrics$roc_auc, na.rm = TRUE)) / (random_sets + 1)
  contrast_pass <- (
    frozen$role[[1]] == "primary" &&
      complete_metrics$roc_auc > 0.60 &&
      complete_metrics$standardized_mean_difference > 0.25 &&
      empirical_p <= 0.05 &&
      positive_sample_directions >= required_positive_samples
  )
  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    contrast_id = contrast_id,
    panel = panel_value,
    coarse_label = coarse_value,
    role = frozen$role[[1]],
    positive_cells = sum(y == 1L),
    negative_cells = sum(y == 0L),
    target_candidate_orthogroups = length(target_ogs),
    mapped_target_orthogroups = sum(target_ogs %in% og_levels),
    complete_program_roc_auc = complete_metrics$roc_auc,
    complete_program_average_precision = complete_metrics$average_precision,
    complete_program_standardized_mean_difference = complete_metrics$standardized_mean_difference,
    random_program_auc_median = median(random_aucs, na.rm = TRUE),
    random_program_auc_q95 = as.numeric(quantile(random_aucs, 0.95, na.rm = TRUE)),
    empirical_p_random_auc_ge_complete = empirical_p,
    positive_sample_directions = positive_sample_directions,
    required_positive_sample_directions = required_positive_samples,
    frozen_contrast_pass = contrast_pass
  )
}

pooled_table <- rbindlist(pooled_rows)
sample_table <- rbindlist(sample_rows)
summary_table <- rbindlist(summary_rows)
complete_pooled <- pooled_table[
  method == "complete_orthogroup_family_program",
  .(
    contrast_id,
    complete_family_roc_auc = roc_auc,
    complete_family_mapped_orthogroups = mapped_candidate_orthogroups
  )
]
single_pooled <- pooled_table[
  method == "target_species_single_gene_orthogroups",
  .(
    contrast_id,
    single_gene_roc_auc = roc_auc,
    single_gene_mapped_orthogroups = mapped_candidate_orthogroups
  )
]
complete_sample_medians <- sample_table[
  method == "complete_orthogroup_family_program",
  .(complete_family_sample_auc_median = median(roc_auc, na.rm = TRUE)),
  by = contrast_id
]
single_sample_medians <- sample_table[
  method == "target_species_single_gene_orthogroups",
  .(single_gene_sample_auc_median = median(roc_auc, na.rm = TRUE)),
  by = contrast_id
]
method_advantage_table <- Reduce(
  function(left, right) merge(left, right, by = "contrast_id", all = TRUE),
  list(complete_pooled, single_pooled, complete_sample_medians, single_sample_medians)
)
method_advantage_table[, pooled_auc_delta_complete_minus_single := complete_family_roc_auc - single_gene_roc_auc]
method_advantage_table[, frozen_method_advantage_supported :=
  single_gene_mapped_orthogroups > 0L &
    pooled_auc_delta_complete_minus_single >= 0.02 &
    complete_family_sample_auc_median > single_gene_sample_auc_median
]
fwrite(pooled_table, file.path(partial_dir, "method_comparison_pooled.tsv"), sep = "\t")
fwrite(sample_table, file.path(partial_dir, "method_comparison_by_sample.tsv"), sep = "\t")
fwrite(summary_table, file.path(partial_dir, "validation_summary.tsv"), sep = "\t")
fwrite(method_advantage_table, file.path(partial_dir, "method_advantage.tsv"), sep = "\t")

primary_summary <- summary_table[role == "primary"]
passing_systems <- primary_summary[frozen_contrast_pass == TRUE]
panels_with_pass <- unique(passing_systems$panel)
majority_reversal <- primary_summary[
  positive_sample_directions < ceiling(required_positive_sample_directions / 2)
]
overall_status <- if (
  nrow(passing_systems) >= 2L &&
    all(c("leaf", "root") %in% panels_with_pass) &&
    nrow(majority_reversal) == 0L
) {
  "positive_external_replication"
} else if (nrow(passing_systems) >= 1L) {
  "mixed_external_support"
} else {
  "external_nonreplication"
}

derived_path <- file.path(derived_dir, "GSE232863_top128_orthogroup_v1_1.rds")
saveRDS(
  list(
    matrix = encoded,
    cell_metadata = selected_meta,
    orthogroups = og_levels,
    top_k = top_k,
    seed = seed,
    source_rds_sha256 = "7a308241e07575770f227057b1cd716b4cb9c255200c4ff4189c567bc2fbc958"
  ),
  derived_path,
  compress = "gzip"
)

audit <- list(
  started_utc = started_utc,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  status = "independent_program_validation_v2_1_complete",
  overall_interpretation = overall_status,
  source_rds = source_rds,
  source_rds_modified = FALSE,
  source_assay = "RNA",
  source_layer = "counts",
  source_expression_class = "nonnegative_integer",
  selected_cells_encoded = nrow(encoded),
  mapped_orthogroups = ncol(encoded),
  top_k = top_k,
  random_sets = random_sets,
  seed = seed,
  primary_contrasts = nrow(primary_summary),
  passing_primary_contrasts = nrow(passing_systems),
  passing_contrast_ids = passing_systems$contrast_id,
  majority_reversal_contrast_ids = majority_reversal$contrast_id,
  derived_path = derived_path,
  implementation_correction = "v2.1 corrects the per-contrast sample-direction filter; the frozen protocol, seed, candidate programs, contrasts, and thresholds are unchanged. v2 is retained for audit but is not authoritative.",
  interpretation_boundary = "Cell-level metrics are descriptive; frozen sample-direction consistency is required. Results validate transfer of pre-discovered orthogroup programs, not causal function or paralog replacement.",
  safety = list(
    source_rds_modified = FALSE,
    freeze_modified_after_scoring = FALSE,
    model_training_started = FALSE,
    v1_outputs_modified = FALSE,
    other_task_written = FALSE
  )
)
write_json(audit, file.path(partial_dir, "audit.json"), pretty = TRUE, auto_unbox = TRUE)

rm(obj, rna, counts)
invisible(gc())
if (!file.rename(partial_dir, run_dir)) stop("Failed to atomically finalize v2.1 run directory")

file.copy(file.path(run_dir, "validation_summary.tsv"), file.path(meta_dir, "independent_program_validation_summary_v2_1.tsv"), overwrite = FALSE)
file.copy(file.path(run_dir, "method_comparison_pooled.tsv"), file.path(meta_dir, "independent_program_method_comparison_pooled_v2_1.tsv"), overwrite = FALSE)
file.copy(file.path(run_dir, "method_comparison_by_sample.tsv"), file.path(meta_dir, "independent_program_method_comparison_by_sample_v2_1.tsv"), overwrite = FALSE)
file.copy(file.path(run_dir, "method_advantage.tsv"), file.path(meta_dir, "independent_program_method_advantage_v2_1.tsv"), overwrite = FALSE)
file.copy(file.path(run_dir, "audit.json"), file.path(report_dir, "independent_program_validation_v2_1_audit.json"), overwrite = FALSE)
cat(toJSON(audit, pretty = TRUE, auto_unbox = TRUE), "\n")
