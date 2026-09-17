cat("R_VERSION\t", R.version.string, "\n", sep = "")
for (package in c("SeuratObject", "Seurat", "Matrix", "jsonlite")) {
  if (requireNamespace(package, quietly = TRUE)) {
    cat(package, "\t", as.character(packageVersion(package)), "\n", sep = "")
  } else {
    cat(package, "\tMISSING\n", sep = "")
  }
}
