#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

project <- "."
run_dir <- "data/runs/PhyloOpenCell/label_resolution_circularity_robustness_v1"
partial_dir <- paste0(run_dir, ".partial")
report_dir <- file.path(project, "reports/label_resolution_circularity_robustness_v1")
freeze_dir <- report_dir

paths <- list(
  protocol = file.path(project, "reports/gse268881_root_xylem_validation_v1/gse268881_root_xylem_validation_protocol_freeze_v1.md"),
  mapping = file.path(project, "metadata/gse268881_three_species_mapping_freeze_v1.tsv"),
  bridge = file.path(project, "metadata/gse268881_project_orthogroup_bridge_v1.tsv"),
  markers = file.path(project, "metadata/coarse_marker_orthogroups_v2.tsv"),
  authority_script = file.path(project, "scripts/run_gse268881_root_xylem_validation_v1.R"),
  authority_summary = "data/runs/PhyloOpenCell/gse268881_root_xylem_validation_v1_1/validation_summary.tsv",
  authority_replicates = "data/runs/PhyloOpenCell/gse268881_root_xylem_validation_v1_1/method_comparison_by_replicate.tsv",
  marker_overlap = file.path(project, "reports/gse268881_esa_feasibility_v1/gse268881_author_marker_overlap_v1.tsv")
)
expected_hashes <- c(
  protocol = "7d93765f465e43f067679e6864cc5495b114867908e724e7a033f52c8611e2b1",
  mapping = "f1410cae32cd687a540647e120d811294697da0cccc6837646f2a68fe77270d4",
  bridge = "69c4da08c8d6415ce479391600ec4fc46e34a2388ae7b4feeadad682e18ea05d",
  markers = "4793e68ffe55e93ca7f18a57058afeb9616bf1d1f6457bba42116995ab60df13",
  authority_script = "1b141a34344e882d8575c4a00a4f62f4784b381781c55719c4dc1ea92eb7bcdb",
  authority_summary = "ddbf41e4f41ca573baab5f53a53116731b614001fe040794f77e9d97df499447",
  authority_replicates = "9bb6a7d91aafe6ec9af08ac4722e4111424b653c44443467524b0757b2c30a13",
  marker_overlap = "63622c11403b07e0107c2a5797b34ef0dd6c96b98605c4d34c6b3dd2be732cf5"
)
datasets <- list(
  esa = list(rds = "data/PhyloOpenCell/external/GSE268881/extracted/210705_Esa_DouRe_Combined_wAnn.RDS", gene_column = "esa_gene"),
  sir = list(rds = "data/PhyloOpenCell/external/GSE268881/extracted/Sir/210822_Sir_DouRe_Combined_wAnn.RDS", gene_column = "sir_gene")
)
k_values <- c(64L, 128L, 256L)
base_seed <- 20260915L
permutation_ids <- 0:10
partial_excluded_ogs <- c("OG0000478", "OG0001133", "OG0002797", "OG0003252", "OG0003815")
all_labels <- c("Atrichoblast", "Columella", "Cortex", "Endodermis", "LRC", "Pericycle", "Phloem", "Procambium", "Trichoblast", "Xylem")
endpoint_definitions <- list(
  xylem_vs_all_other_root = list(positive = "Xylem", negative = setdiff(all_labels, "Xylem")),
  xylem_vs_other_stele = list(positive = "Xylem", negative = c("Pericycle", "Phloem", "Procambium")),
  stele_vascular_vs_non_stele = list(positive = c("Pericycle", "Phloem", "Procambium", "Xylem"),
                                      negative = c("Atrichoblast", "Columella", "Cortex", "Endodermis", "LRC", "Trichoblast"))
)

sha256 <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE)
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

roc_auc <- function(y, score) {
  n_pos <- sum(y == 1L); n_neg <- sum(y == 0L)
  if (!n_pos || !n_neg) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

average_precision <- function(y, score) {
  n_pos <- sum(y == 1L)
  if (!n_pos) return(NA_real_)
  ordering <- order(-score, seq_along(score))
  sorted_y <- y[ordering]
  precision <- cumsum(sorted_y == 1L) / seq_along(sorted_y)
  mean(precision[sorted_y == 1L])
}

standardized_effect <- function(y, score) {
  positive <- score[y == 1L]; negative <- score[y == 0L]
  pooled <- sqrt((var(positive) + var(negative)) / 2)
  if (!is.finite(pooled) || pooled == 0) return(0)
  (mean(positive) - mean(negative)) / pooled
}

metric_record <- function(y, score) {
  list(auroc = roc_auc(y, score), average_precision = average_precision(y, score),
       smd = standardized_effect(y, score), mean_positive = mean(score[y == 1L]), mean_negative = mean(score[y == 0L]))
}

score_order_condition <- function(counts, feature_og, selected_indices, tie_rank, full_target, partial_target) {
  og_levels <- sort(unique(feature_og[!is.na(feature_og) & feature_og != ""]))
  og_to_index <- setNames(seq_along(og_levels), og_levels)
  feature_og_index <- unname(og_to_index[feature_og])
  result_full <- matrix(0, nrow = length(selected_indices), ncol = length(k_values), dimnames = list(NULL, k_values))
  result_partial <- result_full
  for (row_index in seq_along(selected_indices)) {
    source_col <- selected_indices[[row_index]]
    start <- counts@p[[source_col]] + 1L
    end <- counts@p[[source_col + 1L]]
    if (start > end) next
    feature_indices <- counts@i[start:end] + 1L
    values <- counts@x[start:end]
    keep <- is.finite(values) & values > 0
    feature_indices <- feature_indices[keep]; values <- values[keep]
    if (!length(values)) next
    ordered <- order(-values, tie_rank[feature_indices])
    for (k_index in seq_along(k_values)) {
      take <- head(ordered, k_values[[k_index]])
      ranked_features <- feature_indices[take]
      weights <- 1 / log2(seq_along(take) + 1)
      mapped_indices <- feature_og_index[ranked_features]
      keep_mapped <- !is.na(mapped_indices)
      if (!any(keep_mapped)) next
      best <- tapply(weights[keep_mapped], mapped_indices[keep_mapped], max)
      norm_value <- sqrt(sum(best ^ 2))
      if (!is.finite(norm_value) || norm_value == 0) next
      og_names <- og_levels[as.integer(names(best))]
      result_full[row_index, k_index] <- sum(best[og_names %in% full_target]) / norm_value
      result_partial[row_index, k_index] <- sum(best[og_names %in% partial_target]) / norm_value
    }
  }
  list(full = result_full, partial = result_partial)
}

write_failure <- function(message, stage, input_hashes = NULL) {
  dir.create(partial_dir, recursive = TRUE, showWarnings = FALSE)
  write_json(list(status = "FAILED_GATE", completed = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                  stage = stage, message = message, input_hashes = input_hashes,
                  safety = list(source_modified = FALSE, old_results_overwritten = FALSE, model_training_started = FALSE)),
             file.path(partial_dir, "audit.json"), pretty = TRUE, auto_unbox = TRUE)
  stop(message, call. = FALSE)
}

if (dir.exists(run_dir) || dir.exists(partial_dir)) stop("refusing to overwrite run or partial directory")
if (!dir.exists(report_dir)) stop("freeze/report directory must be uploaded before execution")
dir.create(partial_dir, recursive = TRUE)

actual_hashes <- vapply(names(expected_hashes), function(name) sha256(paths[[name]]), character(1))
if (!identical(unname(actual_hashes), unname(expected_hashes))) write_failure("frozen input hash mismatch", "input_hash_gate", as.list(actual_hashes))

mapping <- fread(paths$mapping)
bridge <- fread(paths$bridge)
markers <- fread(paths$markers)
authority <- fread(paths$authority_summary)
authority_replicates <- fread(paths$authority_replicates)
program <- unique(markers[panel == "vascular" & coarse_label == "xylem_lineage" & species == "Arabidopsis thaliana" & within_species_rank <= 18,
                          .(frozen_rank = within_species_rank, orthogroup_id)])
setorder(program, frozen_rank, orthogroup_id)
if (nrow(program) != 18L || uniqueN(program$orthogroup_id) != 18L) write_failure("frozen 18-OG program changed", "program_gate", as.list(actual_hashes))
if (!all(partial_excluded_ogs %in% program$orthogroup_id)) write_failure("partial circularity OG set not contained in frozen program", "circularity_gate", as.list(actual_hashes))

metric_rows <- list(); layer_rows <- list(); dataset_audits <- list(); authority_checks <- list()
for (species_name in names(datasets)) {
  spec <- datasets[[species_name]]
  before <- file.info(spec$rds)[, c("size", "mtime")]
  obj <- readRDS(spec$rds)
  metadata <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
  if (!"RNA" %in% names(obj@assays) || !"counts" %in% Layers(obj[["RNA"]])) write_failure(paste(species_name, "RNA counts missing"), "expression_layer_gate", as.list(actual_hashes))
  counts <- obj[["RNA"]]@counts
  if (!identical(rownames(obj@meta.data), colnames(counts))) write_failure(paste(species_name, "metadata/count order mismatch"), "alignment_gate", as.list(actual_hashes))
  if (any(!is.finite(counts@x)) || any(counts@x < 0)) write_failure(paste(species_name, "RNA counts not finite nonnegative"), "expression_value_gate", as.list(actual_hashes))
  selected <- which(metadata$intspace_treatment == "Ctrl" & !is.na(metadata$intspace_celltype) & metadata$intspace_celltype != "")
  selected_metadata <- metadata[selected, .(cell_id, replicate = as.character(intspace_replicate), source_label = as.character(intspace_celltype))]
  if (!setequal(unique(selected_metadata$replicate), c("R1", "R2"))) write_failure(paste(species_name, "replicate set changed"), "label_gate", as.list(actual_hashes))
  if (!setequal(unique(selected_metadata$source_label), all_labels)) write_failure(paste(species_name, "label set changed"), "label_gate", as.list(actual_hashes))

  gene_map <- setNames(bridge$project_orthogroup, bridge[[spec$gene_column]])
  feature_og <- unname(gene_map[rownames(counts)])
  frozen_species_genes <- mapping[[spec$gene_column]]
  frozen_present <- mapping$project_orthogroup[frozen_species_genes != "" & frozen_species_genes %in% rownames(counts)]
  full_target <- sort(intersect(program$orthogroup_id, unique(frozen_present)))
  partial_target <- setdiff(full_target, partial_excluded_ogs)
  if (length(full_target) != 16L || length(partial_target) != 11L) write_failure(paste(species_name, "mapped full/partial target count changed"), "mapping_gate", as.list(actual_hashes))

  layer_rows[[length(layer_rows) + 1L]] <- data.table(
    species = species_name, source_rds = spec$rds, default_assay = DefaultAssay(obj),
    rna_layers = paste(Layers(obj[["RNA"]]), collapse = ";"),
    integrated_present = "integrated" %in% names(obj@assays),
    integrated_layers = if ("integrated" %in% names(obj@assays)) paste(Layers(obj[["integrated"]]), collapse = ";") else "",
    scoring_assay = "RNA", scoring_layer = "counts", scaled_matrix_used = FALSE
  )

  for (permutation_id in permutation_ids) {
    if (permutation_id == 0L) {
      tie_rank <- seq_len(nrow(counts)); order_seed <- NA_integer_; order_label <- "immutable_source_order"
    } else {
      order_seed <- base_seed + permutation_id
      set.seed(order_seed)
      permutation <- sample(seq_len(nrow(counts)), nrow(counts), replace = FALSE)
      tie_rank <- integer(nrow(counts)); tie_rank[permutation] <- seq_along(permutation)
      order_label <- "deterministic_global_feature_permutation"
    }
    scores <- score_order_condition(counts, feature_og, selected, tie_rank, full_target, partial_target)
    for (program_name in c("complete_18og_family", "partial_published_marker_overlap_exclusion")) {
      score_matrix <- if (program_name == "complete_18og_family") scores$full else scores$partial
      for (k_index in seq_along(k_values)) {
        k <- k_values[[k_index]]; score <- score_matrix[, k_index]
        for (endpoint_id in names(endpoint_definitions)) {
          endpoint <- endpoint_definitions[[endpoint_id]]
          endpoint_mask <- selected_metadata$source_label %in% c(endpoint$positive, endpoint$negative)
          endpoint_y <- as.integer(selected_metadata$source_label[endpoint_mask] %in% endpoint$positive)
          endpoint_score <- score[endpoint_mask]
          endpoint_replicate <- selected_metadata$replicate[endpoint_mask]
          for (scope in c("pooled", "R1", "R2")) {
            scope_mask <- if (scope == "pooled") rep(TRUE, length(endpoint_y)) else endpoint_replicate == scope
            record <- metric_record(endpoint_y[scope_mask], endpoint_score[scope_mask])
            metric_rows[[length(metric_rows) + 1L]] <- data.table(
              species = species_name, endpoint_id = endpoint_id, program = program_name,
              k = k, permutation_id = permutation_id, order_seed = ifelse(is.na(order_seed), "", as.character(order_seed)),
              order_condition = order_label, scope = scope,
              positive_cells = sum(endpoint_y[scope_mask] == 1L), negative_cells = sum(endpoint_y[scope_mask] == 0L),
              auroc = record$auroc, average_precision = record$average_precision, smd = record$smd,
              mean_positive = record$mean_positive, mean_negative = record$mean_negative,
              label_evidence = "author_assigned_label_transfer_derived",
              scoring_layer = "RNA_counts"
            )
          }
        }
      }
    }
  }

  species_metrics <- rbindlist(metric_rows)[species == species_name]
  baseline <- species_metrics[endpoint_id == "xylem_vs_all_other_root" & program == "complete_18og_family" & k == 128L & permutation_id == 0L]
  pooled <- baseline[scope == "pooled"]
  auth <- authority[species == species_name]
  if (nrow(pooled) != 1L || nrow(auth) != 1L) write_failure(paste(species_name, "authority row missing"), "authority_crosscheck", as.list(actual_hashes))
  pooled_ok <- abs(pooled$auroc - auth$complete_program_roc_auc) <= 1e-12 &&
    abs(pooled$average_precision - auth$complete_program_average_precision) <= 1e-12 &&
    abs(pooled$smd - auth$complete_program_standardized_mean_difference) <= 1e-12
  replicate_ok <- all(vapply(c("R1", "R2"), function(rep_name) {
    observed <- baseline[scope == rep_name, auroc]
    expected <- authority_replicates[species == species_name & replicate == rep_name & method == "complete_orthogroup_family_program", roc_auc]
    length(observed) == 1L && length(expected) == 1L && abs(observed - expected) <= 1e-12
  }, logical(1)))
  if (!pooled_ok || !replicate_ok) write_failure(paste(species_name, "K128 source-order v1.1 reproduction failed"), "authority_crosscheck", as.list(actual_hashes))
  authority_checks[[species_name]] <- list(pooled_exact_1e12 = pooled_ok, replicate_auroc_exact_1e12 = replicate_ok)
  rm(obj, counts); invisible(gc())
  after <- file.info(spec$rds)[, c("size", "mtime")]
  dataset_audits[[species_name]] <- list(source_rds = spec$rds, bytes = unname(before$size), source_modified = !identical(before, after),
                                         ctrl_cells = nrow(selected_metadata), mapped_full_target_ogs = length(full_target), mapped_partial_target_ogs = length(partial_target))
}

metrics <- rbindlist(metric_rows)
layers <- rbindlist(layer_rows)
fwrite(metrics, file.path(partial_dir, "nested_label_parameter_metrics.tsv"), sep = "\t")
fwrite(layers, file.path(partial_dir, "expression_layer_audit.tsv"), sep = "\t")

circularity <- data.table(
  sensitivity_id = "partial_published_marker_overlap_exclusion",
  evidence_source = paths$marker_overlap,
  evidence_source_sha256 = expected_hashes[["marker_overlap"]],
  excluded_orthogroups = paste(partial_excluded_ogs, collapse = ";"),
  excluded_frozen_orthogroups = length(partial_excluded_ogs), frozen_program_orthogroups = 18L,
  frozen_program_coverage_fraction = length(partial_excluded_ogs) / 18,
  full_label_transfer_feature_set_available = FALSE,
  residual_risk = "complete label-transfer gene set unavailable; integration features are not treated as label-assignment markers",
  allowed_claim = "robustness to removal of documented published-marker overlap only",
  prohibited_claim = "complete exclusion of label circularity"
)
fwrite(circularity, file.path(partial_dir, "circularity_coverage.tsv"), sep = "\t")

audit <- list(
  status = "COMPLETE", completed = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), version = "label_resolution_circularity_robustness_v1",
  metric_rows = nrow(metrics), expected_metric_rows = 2L * 3L * 2L * 3L * 11L * 3L,
  input_hashes = as.list(actual_hashes), authority_checks = authority_checks, datasets = dataset_audits,
  frozen = list(k_values = k_values, permutation_ids = permutation_ids, base_seed = base_seed,
                endpoints = names(endpoint_definitions), circularity_excluded_ogs = partial_excluded_ogs),
  safety = list(source_objects_modified = any(vapply(dataset_audits, function(x) x$source_modified, logical(1))),
                frozen_protocol_modified = FALSE, labels_modified = FALSE, thresholds_modified = FALSE,
                program_members_modified = FALSE, old_results_overwritten = FALSE, scaled_matrix_used = FALSE,
                model_training_started = FALSE, deep_learning_used = FALSE, new_data_downloaded = FALSE)
)
write_json(audit, file.path(partial_dir, "audit.json"), pretty = TRUE, auto_unbox = TRUE)
if (nrow(metrics) != audit$expected_metric_rows) write_failure("unexpected metric row count", "output_gate", as.list(actual_hashes))
if (!file.rename(partial_dir, run_dir)) stop("atomic finalization failed")
cat(toJSON(audit, pretty = TRUE, auto_unbox = TRUE), "\n")
