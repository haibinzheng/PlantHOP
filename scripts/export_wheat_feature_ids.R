suppressPackageStartupMessages(library(SeuratObject))

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[[1]]
output_path <- args[[2]]

outer_connection <- gzfile(input_path, open = "rb")
connection <- gzcon(outer_connection)
on.exit(close(connection), add = TRUE)
object <- readRDS(connection)

feature_ids <- rownames(object[["RNA"]])
stopifnot(length(feature_ids) == nrow(object[["RNA"]]))
stopifnot(!anyDuplicated(feature_ids))

writeLines(c("gene_id", feature_ids), output_path, useBytes = TRUE)
cat(output_path, "\n")
