#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) == 1L) args[[1L]] else ""
if (!mode %in% c("smoke", "formal", "formal_v1_1")) stop("usage: run_gse268881_root_xylem_validation_v1.R smoke|formal|formal_v1_1")

project <- "."
protocol <- file.path(project, "reports/gse268881_root_xylem_validation_v1/gse268881_root_xylem_validation_protocol_freeze_v1.md")
protocol_expected <- "7d93765f465e43f067679e6864cc5495b114867908e724e7a033f52c8611e2b1"
mapping_freeze <- file.path(project, "metadata/gse268881_three_species_mapping_freeze_v1.tsv")
mapping_expected <- "f1410cae32cd687a540647e120d811294697da0cccc6837646f2a68fe77270d4"
bridge_path <- file.path(project, "metadata/gse268881_project_orthogroup_bridge_v1.tsv")
candidates_path <- file.path(project, "metadata/conserved_program_candidates_v2.tsv")
coarse_markers_path <- file.path(project, "metadata/coarse_marker_orthogroups_v2.tsv")
run_dir <- if (mode == "formal_v1_1") {
  "data/runs/PhyloOpenCell/gse268881_root_xylem_validation_v1_1"
} else {
  "data/runs/PhyloOpenCell/gse268881_root_xylem_validation_v1"
}
partial_dir <- paste0(run_dir, ".partial")
smoke_dir <- "data/runs/PhyloOpenCell/gse268881_root_xylem_validation_smoke_v1"
top_k <- 128L
random_programs <- 500L
random_seed <- 20260915L
overlap_ogs <- c("OG0000478", "OG0001133", "OG0002797", "OG0003252", "OG0003815")

sha256 <- function(path) {
  unname(tools::md5sum(path))
}

sha256_external <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE)
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

roc_auc <- function(y, score) {
  y <- as.integer(y)
  n_pos <- sum(y == 1L)
  n_neg <- sum(y == 0L)
  if (!n_pos || !n_neg) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

average_precision <- function(y, score) {
  y <- as.integer(y)
  n_pos <- sum(y == 1L)
  if (!n_pos) return(NA_real_)
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
    start <- counts@p[[source_col]] + 1L
    end <- counts@p[[source_col + 1L]]
    if (start > end) next
    feature_indices <- counts@i[start:end] + 1L
    values <- counts@x[start:end]
    keep_positive <- is.finite(values) & values > 0
    feature_indices <- feature_indices[keep_positive]
    values <- values[keep_positive]
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
    if (any(keep_representative)) {
      single_score[[row_index]] <- sum(weights[keep_representative]) / norm_value
    }
  }

  encoded <- sparseMatrix(
    i = unlist(output_i, use.names = FALSE),
    j = unlist(output_j, use.names = FALSE),
    x = unlist(output_x, use.names = FALSE),
    dims = c(length(selected_indices), length(og_levels))
  )
  list(matrix = encoded, orthogroups = og_levels, single_score = single_score)
}

run_smoke <- function() {
  if (dir.exists(smoke_dir)) stop("refusing to overwrite smoke output")
  dir.create(smoke_dir, recursive = TRUE)
  y <- c(0L, 0L, 1L, 1L)
  score <- c(0.1, 0.2, 0.8, 0.9)
  stopifnot(abs(roc_auc(y, score) - 1) < 1e-12)
  stopifnot(abs(average_precision(y, score) - 1) < 1e-12)
  stopifnot(abs(roc_auc(y, rep(1, 4)) - 0.5) < 1e-12)
  toy <- sparseMatrix(
    i = c(1, 2, 3, 1, 2, 3, 4),
    j = c(1, 1, 1, 2, 2, 3, 4),
    x = c(9, 7, 1, 2, 8, 6, 10),
    dims = c(4, 4)
  )
  encoded <- encode_counts(
    toy,
    c("OG1", "OG1", "OG2", NA_character_),
    c("OG1", NA_character_, "OG2", NA_character_),
    1:4
  )
  norms <- sqrt(rowSums(encoded$matrix ^ 2))
  stopifnot(all(abs(norms[norms > 0] - 1) < 1e-12))
  set.seed(random_seed)
  draw1 <- sample(letters, 5)
  set.seed(random_seed)
  draw2 <- sample(letters, 5)
  stopifnot(identical(draw1, draw2))
  audit <- list(
    status = "gse268881_root_xylem_validation_smoke_passed",
    completed = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    tests = c("auc_direction", "average_precision", "auc_ties", "top_rank_encoding", "l2_normalization", "seed_reproducibility"),
    real_expression_read = FALSE,
    real_labels_read = FALSE,
    program_scores_computed = FALSE
  )
  write_json(audit, file.path(smoke_dir, "audit.json"), pretty = TRUE, auto_unbox = TRUE)
  cat(toJSON(audit, pretty = TRUE, auto_unbox = TRUE), "\n")
}

run_formal <- function() {
  if (dir.exists(run_dir) || dir.exists(partial_dir)) stop("refusing to overwrite formal output")
  if (!file.exists(file.path(smoke_dir, "audit.json"))) stop("smoke audit is required")
  if (sha256_external(protocol) != protocol_expected) stop("protocol hash changed")
  if (sha256_external(mapping_freeze) != mapping_expected) stop("mapping freeze hash changed")
  dir.create(partial_dir, recursive = TRUE)

  bridge <- fread(bridge_path)
  mapping <- fread(mapping_freeze)
  candidates <- fread(candidates_path)
  coarse_markers <- fread(coarse_markers_path)
  target_candidates <- coarse_markers[
    panel == "vascular" & coarse_label == "xylem_lineage" &
      species == "Arabidopsis thaliana" & within_species_rank <= 18,
    orthogroup_id
  ]
  off_target_candidates <- unique(candidates[panel == "vascular" & coarse_label == "phloem_lineage", orthogroup_id])
  if (length(target_candidates) != 18L) stop("frozen xylem candidate count changed")

  datasets <- list(
    esa = list(
      rds = "data/PhyloOpenCell/external/GSE268881/extracted/210705_Esa_DouRe_Combined_wAnn.RDS",
      gene_column = "esa_gene", expected_replicates = c("R1", "R2")
    ),
    sir = list(
      rds = "data/PhyloOpenCell/external/GSE268881/extracted/Sir/210822_Sir_DouRe_Combined_wAnn.RDS",
      gene_column = "sir_gene", expected_replicates = c("R1", "R2")
    )
  )

  pooled_rows <- list()
  replicate_rows <- list()
  summary_rows <- list()
  random_rows <- list()
  manifests <- list()
  set.seed(random_seed)

  for (species in names(datasets)) {
    spec <- datasets[[species]]
    obj <- readRDS(spec$rds)
    metadata <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
    counts <- obj[["RNA"]]@counts
    if (!identical(rownames(obj@meta.data), colnames(counts))) stop(paste(species, "metadata/count order mismatch"))
    if (any(!is.finite(counts@x)) || any(counts@x < 0)) stop(paste(species, "RNA counts are not finite nonnegative"))
    selected <- which(metadata$intspace_treatment == "Ctrl" & !is.na(metadata$intspace_celltype) & metadata$intspace_celltype != "")
    selected_metadata <- metadata[selected, .(
      cell_id,
      replicate = as.character(intspace_replicate),
      source_label = as.character(intspace_celltype)
    )]
    if (!setequal(unique(selected_metadata$replicate), spec$expected_replicates)) stop(paste(species, "unexpected Ctrl replicates"))

    gene_map <- setNames(bridge$project_orthogroup, bridge[[spec$gene_column]])
    feature_og <- unname(gene_map[rownames(counts)])
    representative_flag <- as.logical(mapping$single_gene_representative)
    representative_rows <- mapping[representative_flag & get(spec$gene_column) != ""]
    representative_lookup <- setNames(representative_rows$project_orthogroup, representative_rows[[spec$gene_column]])
    representative_feature_og <- unname(representative_lookup[rownames(counts)])
    encoded <- encode_counts(counts, feature_og, representative_feature_og, selected)
    matrix <- encoded$matrix
    og_to_col <- setNames(seq_along(encoded$orthogroups), encoded$orthogroups)
    available_ogs <- encoded$orthogroups

    frozen_species_genes <- mapping[[spec$gene_column]]
    frozen_present <- mapping$project_orthogroup[frozen_species_genes != "" & frozen_species_genes %in% rownames(counts)]
    target <- sort(intersect(target_candidates, unique(frozen_present)))
    anti_target <- sort(setdiff(target, overlap_ogs))
    off_target <- sort(intersect(off_target_candidates, available_ogs))
    if (length(target) != 16L || length(anti_target) != 11L) stop(paste(species, "frozen target mapping changed"))
    single_present_ogs <- unique(representative_rows[
      get(spec$gene_column) %in% rownames(counts), project_orthogroup
    ])
    if (!length(single_present_ogs) || !any(encoded$single_score > 0)) {
      stop(paste(species, "single-gene comparator is unexpectedly empty"))
    }

    labels <- selected_metadata$source_label
    replicates <- selected_metadata$replicate
    y <- as.integer(labels == "Xylem")
    target_cols <- unname(og_to_col[target])
    anti_cols <- unname(og_to_col[anti_target])
    off_cols <- unname(og_to_col[off_target])
    complete_score <- as.numeric(rowSums(matrix[, target_cols, drop = FALSE]))
    anti_score <- as.numeric(rowSums(matrix[, anti_cols, drop = FALSE]))
    off_score <- if (length(off_cols)) as.numeric(rowSums(matrix[, off_cols, drop = FALSE])) else rep(0, length(y))
    binary <- matrix[, target_cols, drop = FALSE]
    if (length(binary@x)) binary@x[] <- 1
    equal_score <- as.numeric(rowSums(binary)) / length(target)
    method_scores <- list(
      complete_orthogroup_family_program = complete_score,
      deterministic_single_gene_comparator = encoded$single_score,
      anti_circularity_program = anti_score,
      phloem_off_target_program = off_score,
      omg_style_equal_orthogroup_presence = equal_score
    )
    method_sizes <- c(length(target), length(single_present_ogs), length(anti_target), length(off_target), length(target))
    names(method_sizes) <- names(method_scores)

    species_replicate_aucs <- list()
    for (method in names(method_scores)) {
      score <- method_scores[[method]]
      record <- metric_record(y, score)
      pooled_rows[[length(pooled_rows) + 1L]] <- data.table(
        species = species, method = method, orthogroups = method_sizes[[method]],
        positive_cells = sum(y == 1L), negative_cells = sum(y == 0L),
        roc_auc = record$roc_auc, average_precision = record$average_precision,
        standardized_mean_difference = record$standardized_mean_difference,
        mean_positive = record$mean_positive, mean_negative = record$mean_negative
      )
      replicate_aucs <- numeric()
      for (replicate in spec$expected_replicates) {
        mask <- replicates == replicate
        rep_record <- metric_record(y[mask], score[mask])
        replicate_rows[[length(replicate_rows) + 1L]] <- data.table(
          species = species, replicate = replicate, method = method,
          positive_cells = sum(y[mask] == 1L), negative_cells = sum(y[mask] == 0L),
          roc_auc = rep_record$roc_auc, average_precision = rep_record$average_precision,
          standardized_mean_difference = rep_record$standardized_mean_difference,
          mean_positive = rep_record$mean_positive, mean_negative = rep_record$mean_negative
        )
        replicate_aucs <- c(replicate_aucs, rep_record$roc_auc)
      }
      species_replicate_aucs[[method]] <- replicate_aucs
    }

    random_pool <- sort(setdiff(available_ogs, union(target, off_target)))
    if (length(random_pool) < length(target)) stop("random pool too small")
    random_aucs <- numeric(random_programs)
    for (draw in seq_len(random_programs)) {
      chosen <- sample(random_pool, length(target), replace = FALSE)
      score <- as.numeric(rowSums(matrix[, unname(og_to_col[chosen]), drop = FALSE]))
      random_aucs[[draw]] <- roc_auc(y, score)
      random_rows[[length(random_rows) + 1L]] <- data.table(
        species = species, draw = draw, orthogroups = paste(chosen, collapse = ";"), roc_auc = random_aucs[[draw]]
      )
    }
    primary <- metric_record(y, complete_score)
    empirical_p <- (1 + sum(random_aucs >= primary$roc_auc)) / (random_programs + 1)
    primary_rep_aucs <- species_replicate_aucs$complete_orthogroup_family_program
    species_pass <- primary$roc_auc > 0.60 &&
      primary$standardized_mean_difference > 0.25 &&
      empirical_p <= 0.05 && all(primary_rep_aucs > 0.50)
    single_pooled <- metric_record(y, encoded$single_score)
    single_rep_aucs <- species_replicate_aucs$deterministic_single_gene_comparator
    family_advantage <- primary$roc_auc - single_pooled$roc_auc >= 0.02 &&
      median(primary_rep_aucs) > median(single_rep_aucs)
    summary_rows[[length(summary_rows) + 1L]] <- data.table(
      species = species, mapped_target_orthogroups = length(target),
      positive_cells = sum(y == 1L), negative_cells = sum(y == 0L),
      complete_program_roc_auc = primary$roc_auc,
      complete_program_average_precision = primary$average_precision,
      complete_program_standardized_mean_difference = primary$standardized_mean_difference,
      empirical_p_random_auc_ge_complete = empirical_p,
      replicate_1_auc = primary_rep_aucs[[1L]], replicate_2_auc = primary_rep_aucs[[2L]],
      species_primary_pass = species_pass,
      single_gene_roc_auc = single_pooled$roc_auc,
      complete_minus_single_auc = primary$roc_auc - single_pooled$roc_auc,
      family_advantage_supported = family_advantage
    )
    manifests[[species]] <- list(
      source_rds = spec$rds, source_rds_bytes = unname(file.info(spec$rds)$size),
      ctrl_cells = length(selected), xylem_cells = sum(y == 1L),
      mapped_orthogroups = length(available_ogs), target_orthogroups = length(target),
      anti_circularity_orthogroups = length(anti_target), off_target_orthogroups = length(off_target),
      top_k = top_k, source_modified = FALSE
    )
    rm(obj, counts, matrix, encoded)
    invisible(gc())
  }

  pooled <- rbindlist(pooled_rows)
  by_replicate <- rbindlist(replicate_rows)
  summary <- rbindlist(summary_rows)
  random <- rbindlist(random_rows)
  fwrite(pooled, file.path(partial_dir, "method_comparison_pooled.tsv"), sep = "\t")
  fwrite(by_replicate, file.path(partial_dir, "method_comparison_by_replicate.tsv"), sep = "\t")
  fwrite(summary, file.path(partial_dir, "validation_summary.tsv"), sep = "\t")
  fwrite(random, file.path(partial_dir, "random_program_draws.tsv"), sep = "\t")
  passing <- summary[species_primary_pass == TRUE, species]
  family_passing <- summary[family_advantage_supported == TRUE, species]
  overall <- if (length(passing) == 2L) "positive_multi_species_xylem_replication" else if (length(passing)) "mixed_multi_species_support" else "multi_species_nonreplication"
  audit <- list(
    status = if (mode == "formal_v1_1") "gse268881_root_xylem_formal_validation_v1_1_complete" else "gse268881_root_xylem_formal_validation_complete",
    completed = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    overall_interpretation = overall,
    eligible_species = names(datasets),
    passing_species = passing,
    family_advantage_species = family_passing,
    protocol_sha256 = sha256_external(protocol),
    mapping_freeze_sha256 = sha256_external(mapping_freeze),
    bridge_sha256 = sha256_external(bridge_path),
    supersedes_for_method_comparison = if (mode == "formal_v1_1") "data/runs/PhyloOpenCell/gse268881_root_xylem_validation_v1" else NULL,
    top_k = top_k, random_programs = random_programs, random_seed = random_seed,
    datasets = manifests,
    safety = list(
      source_objects_modified = FALSE,
      frozen_protocol_modified = FALSE,
      all_eligible_results_retained = TRUE,
      model_training_started = FALSE,
      deep_learning_used = FALSE
    )
  )
  write_json(audit, file.path(partial_dir, "audit.json"), pretty = TRUE, auto_unbox = TRUE)
  if (!file.rename(partial_dir, run_dir)) stop("failed to atomically finalize run")
  cat(toJSON(audit, pretty = TRUE, auto_unbox = TRUE), "\n")
}

if (mode == "smoke") run_smoke() else run_formal()
