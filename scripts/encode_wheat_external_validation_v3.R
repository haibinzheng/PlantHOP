#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1L) args[[1]] else "smoke"
if (!mode %in% c("smoke", "formal")) stop("mode must be smoke or formal")

project <- "."
source_rds <- "data/external/GSE270342_Triticum_aestivum/GSE270342_seuratObj_for_publication.rds.gz"
mapping_path <- file.path(project, "metadata", "external_validation_v3_wheat_feature_orthogroup_map_v1.tsv")
addendum_path <- file.path(project, "reports", "external_validation_v3_execution_addendum_v1.json")
output_dir <- file.path("data/derived/PhyloOpenCell/external_validation_v3", paste0("wheat_top128_", mode))
partial_dir <- paste0(output_dir, ".partial")
top_k <- 128L

if (dir.exists(output_dir) || dir.exists(partial_dir)) stop("refusing to overwrite wheat encoding output")
if (!file.exists(addendum_path)) stop("execution addendum is required")
dir.create(partial_dir, recursive = TRUE)

outer_connection <- gzfile(source_rds, open = "rb")
connection <- gzcon(outer_connection)
object <- readRDS(connection)
close(connection)
counts <- object[["RNA"]]@counts
metadata <- as.data.table(object@meta.data, keep.rownames = "cell_id")

mapping <- fread(mapping_path)
accepted <- mapping[grepl("^unanimous_", mapping_status) & nzchar(project_orthogroup_id)]
map_lookup <- setNames(accepted$project_orthogroup_id, accepted$wheat_gene_id)
feature_ogs <- unname(map_lookup[rownames(counts)])
feature_ogs[is.na(feature_ogs) | feature_ogs == ""] <- NA_character_
og_levels <- sort(unique(feature_ogs[!is.na(feature_ogs)]))
og_to_index <- setNames(seq_along(og_levels), og_levels)
feature_og_index <- unname(og_to_index[feature_ogs])

selected_indices <- if (mode == "smoke") {
  unique(unlist(lapply(split(seq_len(nrow(metadata)), metadata$orig.ident), head, 80L)))
} else {
  seq_len(ncol(counts))
}

j_parts <- vector("list", length(selected_indices))
x_parts <- vector("list", length(selected_indices))
for (row_index in seq_along(selected_indices)) {
  source_col <- selected_indices[[row_index]]
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
  ordering <- head(order(-values, feature_indices), top_k)
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
  i = rep.int(seq_along(selected_indices), entry_counts),
  j = unlist(j_parts, use.names = FALSE),
  x = unlist(x_parts, use.names = FALSE),
  dims = c(length(selected_indices), length(og_levels))
)

writeMM(encoded, file.path(partial_dir, "encoded.mtx"))
writeLines(og_levels, file.path(partial_dir, "orthogroups.txt"), useBytes = TRUE)
selected_metadata <- metadata[selected_indices, .(
  row_index = seq_along(selected_indices),
  cell_id,
  replicate = orig.ident,
  source_label = as.character(annotation)
)]
fwrite(selected_metadata, file.path(partial_dir, "cell_metadata.tsv"), sep = "\t")
manifest <- list(
  status = paste0("wheat_external_validation_v3_", mode, "_encoded"),
  mode = mode,
  source_rds = source_rds,
  source_rds_sha256 = "7c0054b23b782e538ead563fb1cbe1688db89f185753cf5ff14ccd2cb72fd58f",
  mapping_path = mapping_path,
  mapping_sha256 = "20b2fad7145d838a91dea634fe611b76a40892e89c3f3f3957eb943aee4acaf8",
  cells = nrow(encoded),
  mapped_orthogroups = ncol(encoded),
  nonzero = length(encoded@x),
  top_k = top_k,
  expression_values_used_for_mapping_or_selection = FALSE,
  source_modified = FALSE
)
write_json(manifest, file.path(partial_dir, "manifest.json"), pretty = TRUE, auto_unbox = TRUE)

rm(object, counts)
invisible(gc())
if (!file.rename(partial_dir, output_dir)) stop("failed to finalize wheat encoding")
cat(toJSON(manifest, pretty = TRUE, auto_unbox = TRUE), "\n")
