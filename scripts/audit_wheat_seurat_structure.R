suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[[1]]
output_path <- args[[2]]

outer_connection <- gzfile(input_path, open = "rb")
connection <- gzcon(outer_connection)
on.exit(close(connection), add = TRUE)
object <- readRDS(connection)
metadata <- object[[]]

matrix_summary <- function(matrix) {
  values <- if (inherits(matrix, "sparseMatrix")) matrix@x else as.vector(matrix)
  sampled <- if (length(values) > 1000000) values[seq_len(1000000)] else values
  list(
    dimensions = as.integer(dim(matrix)),
    class = class(matrix),
    nonzero = if (inherits(matrix, "sparseMatrix")) as.integer(length(matrix@x)) else as.integer(sum(matrix != 0)),
    min_nonzero = if (length(values)) as.numeric(min(values, na.rm = TRUE)) else NULL,
    max = if (length(values)) as.numeric(max(values, na.rm = TRUE)) else NULL,
    nonnegative = if (length(values)) isTRUE(min(values, na.rm = TRUE) >= 0) else TRUE,
    integer_like_sample = if (length(sampled)) isTRUE(all(abs(sampled - round(sampled)) < 1e-8)) else TRUE
  )
}

table_to_named_lists <- function(tab) {
  result <- lapply(seq_len(nrow(tab)), function(index) {
    values <- as.list(as.integer(tab[index, ]))
    names(values) <- colnames(tab)
    values
  })
  names(result) <- rownames(tab)
  result
}

assay_summary <- list()
for (assay_name in Assays(object)) {
  assay <- object[[assay_name]]
  layer_names <- Layers(assay)
  layers <- list()
  for (layer_name in layer_names) {
    layers[[layer_name]] <- matrix_summary(LayerData(assay, layer = layer_name))
  }
  assay_summary[[assay_name]] <- list(
    class = class(assay),
    dimensions = as.integer(dim(assay)),
    layers = layers
  )
}

metadata_summary <- list()
for (column in colnames(metadata)) {
  values <- metadata[[column]]
  unique_count <- length(unique(values))
  entry <- list(class = class(values), unique = as.integer(unique_count), missing = as.integer(sum(is.na(values))))
  if (unique_count <= 60) {
    counts <- sort(table(as.character(values), useNA = "ifany"), decreasing = TRUE)
    entry$counts <- as.list(as.integer(counts))
    names(entry$counts) <- names(counts)
  }
  metadata_summary[[column]] <- entry
}

report <- list(
  input_path = input_path,
  file_size_bytes = as.numeric(file.info(input_path)$size),
  object_class = class(object),
  object_size_bytes_in_memory = as.numeric(object.size(object)),
  cells = as.integer(ncol(object)),
  features_default_assay = as.integer(nrow(object)),
  default_assay = DefaultAssay(object),
  assays = assay_summary,
  reductions = Reductions(object),
  metadata_columns = colnames(metadata),
  metadata = metadata_summary,
  cell_name_examples = head(colnames(object), 10),
  gene_name_examples_RNA = head(rownames(object[["RNA"]]), 20),
  duplicate_gene_names_RNA = as.integer(sum(duplicated(rownames(object[["RNA"]])))),
  annotation_by_replicate = unclass(as.data.frame.matrix(table(metadata$annotation, metadata$orig.ident))),
  annotation_by_replicate_named = table_to_named_lists(table(metadata$annotation, metadata$orig.ident)),
  annotation_by_cluster = unclass(as.data.frame.matrix(table(metadata$seurat_clusters, metadata$annotation))),
  RNA_feature_metadata_columns = colnames(object[["RNA"]]@meta.features),
  RNA_feature_metadata_examples = head(object[["RNA"]]@meta.features, 20)
)

write_json(report, output_path, pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)
cat(output_path, "\n")
