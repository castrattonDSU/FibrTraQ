#!/usr/bin/env Rscript

usage <- function() {
  cat("Usage: fibrtraq-ismrm2015-manifests.R <scoring-dir> <manifest-output-dir> [--run-output-root <path>] [--bundles <comma-list>] [--force]\n")
  quit(status = 1)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || "--help" %in% args || "-h" %in% args) usage()

scoring_dir <- args[1]
output_dir <- args[2]
opts <- args[-c(1, 2)]

value_after <- function(flag) {
  idx <- which(opts == flag)
  if (!length(idx)) return(NULL)
  if (idx[1] == length(opts)) stop(flag, " requires a value.")
  opts[idx[1] + 1L]
}

run_output_root <- value_after("--run-output-root")
if (is.null(run_output_root)) {
  matched <- grep("^--run-output-root=", opts, value = TRUE)
  if (length(matched)) run_output_root <- sub("^--run-output-root=", "", matched[1])
}

bundle_arg <- value_after("--bundles")
if (is.null(bundle_arg)) {
  matched <- grep("^--bundles=", opts, value = TRUE)
  if (length(matched)) bundle_arg <- sub("^--bundles=", "", matched[1])
}
bundles <- if (!is.null(bundle_arg) && nzchar(bundle_arg)) {
  trimws(strsplit(bundle_arg, ",", fixed = TRUE)[[1]])
} else {
  NULL
}

suppressPackageStartupMessages(library(FibrTraQ))

out <- ftq_make_ismrm2015_manifests(
  scoring_dir = scoring_dir,
  output_dir = output_dir,
  bundles = bundles,
  run_output_root = run_output_root,
  force = "--force" %in% opts
)

cat("Generated ISMRM 2015 FibrTraQ manifests\n")
cat("  manifests: ", nrow(out), "\n", sep = "")
cat("  output:    ", normalizePath(output_dir, mustWork = FALSE), "\n", sep = "")
print(out[, c("bundle", "manifest_path")], row.names = FALSE)
