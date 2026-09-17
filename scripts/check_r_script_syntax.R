#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Expected one R script path")
parse(file = args[[1]])
cat("PARSE_OK\n")
