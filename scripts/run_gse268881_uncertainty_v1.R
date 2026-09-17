#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

project <- "/workspace/projects/PhyloOpenCell"
output_dir <- "/data/runs/PhyloOpenCell/gse268881_root_xylem_uncertainty_v1"
partial_dir <- paste0(output_dir, ".partial")
plan_path <- file.path(project, "reports/gse268881_root_xylem_validation_v1/gse268881_uncertainty_analysis_plan_v1.md")
mapping_path <- file.path(project, "metadata/gse268881_three_species_mapping_freeze_v1.tsv")
bridge_path <- file.path(project, "metadata/gse268881_project_orthogroup_bridge_v1.tsv")
coarse_path <- file.path(project, "metadata/coarse_marker_orthogroups_v2.tsv")
top_k <- 128L
bootstrap_n <- 2000L
bootstrap_seed <- 20260916L

sha256 <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE)
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

roc_auc <- function(y, score) {
  n_pos <- sum(y == 1L); n_neg <- sum(y == 0L)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

average_precision <- function(y, score) {
  ordering <- order(-score, seq_along(score))
  sorted_y <- y[ordering]
  precision <- cumsum(sorted_y == 1L) / seq_along(sorted_y)
  mean(precision[sorted_y == 1L])
}

smd <- function(y, score) {
  positive <- score[y == 1L]; negative <- score[y == 0L]
  pooled <- sqrt((var(positive) + var(negative)) / 2)
  if (!is.finite(pooled) || pooled == 0) return(0)
  (mean(positive) - mean(negative)) / pooled
}

metrics <- function(y, score) c(auc = roc_auc(y, score), ap = average_precision(y, score), smd = smd(y, score))

encode_counts <- function(counts, feature_og, representative_feature_og, selected_indices) {
  og_levels <- sort(unique(feature_og[!is.na(feature_og) & feature_og != ""]))
  og_to_index <- setNames(seq_along(og_levels), og_levels)
  feature_og_index <- unname(og_to_index[feature_og])
  output_i <- vector("list", length(selected_indices))
  output_j <- vector("list", length(selected_indices))
  output_x <- vector("list", length(selected_indices))
  single_score <- numeric(length(selected_indices))
  for (row_index in seq_along(selected_indices)) {
    source_col <- selected_indices[[row_index]]
    start <- counts@p[[source_col]] + 1L; end <- counts@p[[source_col + 1L]]
    if (start > end) next
    feature_indices <- counts@i[start:end] + 1L
    values <- counts@x[start:end]
    keep <- is.finite(values) & values > 0
    feature_indices <- feature_indices[keep]; values <- values[keep]
    if (!length(values)) next
    ordering <- head(order(-values, feature_indices), top_k)
    ranked_features <- feature_indices[ordering]
    weights <- 1 / log2(seq_along(ordering) + 1)
    mapped_indices <- feature_og_index[ranked_features]
    keep_mapped <- !is.na(mapped_indices)
    if (!any(keep_mapped)) next
    best <- tapply(weights[keep_mapped], mapped_indices[keep_mapped], max)
    norm_value <- sqrt(sum(best ^ 2))
    if (!is.finite(norm_value) || norm_value == 0) next
    output_i[[row_index]] <- rep.int(row_index, length(best))
    output_j[[row_index]] <- as.integer(names(best))
    output_x[[row_index]] <- as.numeric(best) / norm_value
    representative <- representative_feature_og[ranked_features]
    keep_representative <- !is.na(representative) & representative != ""
    if (any(keep_representative)) single_score[[row_index]] <- sum(weights[keep_representative]) / norm_value
  }
  encoded <- sparseMatrix(
    i = unlist(output_i, use.names = FALSE), j = unlist(output_j, use.names = FALSE),
    x = unlist(output_x, use.names = FALSE), dims = c(length(selected_indices), length(og_levels))
  )
  list(matrix = encoded, orthogroups = og_levels, single_score = single_score)
}

stratified_indices <- function(y, replicate) {
  strata <- interaction(replicate, y, drop = TRUE, lex.order = TRUE)
  unlist(lapply(split(seq_along(y), strata), function(index) sample(index, length(index), replace = TRUE)), use.names = FALSE)
}

summarize_scope <- function(species, scope, index, y, replicate, complete, single) {
  yy <- y[index]; rr <- replicate[index]; cc <- complete[index]; ss <- single[index]
  complete_est <- metrics(yy, cc); single_est <- metrics(yy, ss)
  estimate <- c(complete_est, single_est, complete_est - single_est)
  names(estimate) <- c(paste0("complete_", names(complete_est)), paste0("single_", names(single_est)), paste0("delta_", names(complete_est)))
  boot <- matrix(NA_real_, nrow = bootstrap_n, ncol = length(estimate), dimnames = list(NULL, names(estimate)))
  for (iteration in seq_len(bootstrap_n)) {
    sampled <- stratified_indices(yy, rr)
    cm <- metrics(yy[sampled], cc[sampled]); sm <- metrics(yy[sampled], ss[sampled])
    boot[iteration, ] <- c(cm, sm, cm - sm)
  }
  intervals <- apply(boot, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  data.table(
    species = species, scope = scope, estimand = names(estimate), estimate = as.numeric(estimate),
    lower_95 = as.numeric(intervals[1, ]), upper_95 = as.numeric(intervals[2, ]),
    bootstrap_n = bootstrap_n, interval_type = "conditional_stratified_cell_bootstrap_percentile"
  )
}

if (dir.exists(output_dir) || dir.exists(partial_dir)) stop("refusing to overwrite uncertainty output")
dir.create(partial_dir, recursive = TRUE)
mapping <- fread(mapping_path)
bridge <- fread(bridge_path)
coarse <- fread(coarse_path)
target_candidates <- coarse[
  panel == "vascular" & coarse_label == "xylem_lineage" & species == "Arabidopsis thaliana" & within_species_rank <= 18,
  orthogroup_id
]
if (length(target_candidates) != 18L) stop("target definition changed")

datasets <- list(
  esa = list(rds = "/data/datasets/PhyloOpenCell/external/GSE268881/extracted/210705_Esa_DouRe_Combined_wAnn.RDS", gene_column = "esa_gene"),
  sir = list(rds = "/data/datasets/PhyloOpenCell/external/GSE268881/extracted/Sir/210822_Sir_DouRe_Combined_wAnn.RDS", gene_column = "sir_gene")
)
all_intervals <- list(); replicate_points <- list(); set.seed(bootstrap_seed)

for (species in names(datasets)) {
  spec <- datasets[[species]]
  obj <- readRDS(spec$rds); counts <- obj[["RNA"]]@counts
  metadata <- as.data.table(obj@meta.data)
  selected <- which(metadata$intspace_treatment == "Ctrl" & !is.na(metadata$intspace_celltype) & metadata$intspace_celltype != "")
  y <- as.integer(metadata$intspace_celltype[selected] == "Xylem")
  replicate <- as.character(metadata$intspace_replicate[selected])
  gene_map <- setNames(bridge$project_orthogroup, bridge[[spec$gene_column]])
  feature_og <- unname(gene_map[rownames(counts)])
  representatives <- mapping[as.logical(single_gene_representative) & get(spec$gene_column) != ""]
  representative_lookup <- setNames(representatives$project_orthogroup, representatives[[spec$gene_column]])
  representative_feature_og <- unname(representative_lookup[rownames(counts)])
  encoded <- encode_counts(counts, feature_og, representative_feature_og, selected)
  frozen_present <- mapping$project_orthogroup[mapping[[spec$gene_column]] != "" & mapping[[spec$gene_column]] %in% rownames(counts)]
  target <- sort(intersect(target_candidates, unique(frozen_present)))
  if (length(target) != 16L) stop(paste(species, "target mapping changed"))
  target_cols <- match(target, encoded$orthogroups)
  if (anyNA(target_cols)) stop(paste(species, "target absent from encoded matrix"))
  complete <- as.numeric(rowSums(encoded$matrix[, target_cols, drop = FALSE]))
  single <- encoded$single_score
  if (!any(single > 0)) stop(paste(species, "empty single comparator"))
  scopes <- c(list(pooled = seq_along(y)), lapply(sort(unique(replicate)), function(value) which(replicate == value)))
  names(scopes)[-1L] <- sort(unique(replicate))
  for (scope in names(scopes)) {
    all_intervals[[length(all_intervals) + 1L]] <- summarize_scope(species, scope, scopes[[scope]], y, replicate, complete, single)
  }
  for (rep_value in sort(unique(replicate))) {
    index <- which(replicate == rep_value)
    cm <- metrics(y[index], complete[index]); sm <- metrics(y[index], single[index])
    replicate_points[[length(replicate_points) + 1L]] <- data.table(
      species = species, replicate = rep_value, positive_cells = sum(y[index] == 1L), negative_cells = sum(y[index] == 0L),
      complete_auc = cm[["auc"]], complete_ap = cm[["ap"]], complete_smd = cm[["smd"]],
      single_auc = sm[["auc"]], single_ap = sm[["ap"]], single_smd = sm[["smd"]],
      delta_auc = cm[["auc"]] - sm[["auc"]], delta_ap = cm[["ap"]] - sm[["ap"]], delta_smd = cm[["smd"]] - sm[["smd"]]
    )
  }
  rm(obj, counts, encoded); invisible(gc())
}

intervals <- rbindlist(all_intervals); points <- rbindlist(replicate_points)
fwrite(intervals, file.path(partial_dir, "conditional_bootstrap_intervals.tsv"), sep = "\t")
fwrite(points, file.path(partial_dir, "biological_replicate_point_estimates.tsv"), sep = "\t")
audit <- list(
  status = "gse268881_uncertainty_v1_complete", completed = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  plan_sha256 = sha256(plan_path), script_sha256 = sha256(commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))] |> sub("^--file=", "", x = _)),
  mapping_sha256 = sha256(mapping_path), bridge_sha256 = sha256(bridge_path),
  bootstrap_n = bootstrap_n, bootstrap_seed = bootstrap_seed,
  interval_scope = "conditional_on_observed_biological_replicates",
  biological_replicate_population_ci_claimed = FALSE,
  source_objects_modified = FALSE, model_training_started = FALSE, deep_learning_used = FALSE
)
write_json(audit, file.path(partial_dir, "audit.json"), pretty = TRUE, auto_unbox = TRUE)
if (!file.rename(partial_dir, output_dir)) stop("failed to finalize uncertainty output")
cat(toJSON(audit, pretty = TRUE, auto_unbox = TRUE), "\n")
