#!/usr/bin/env Rscript

usage <- function() {
  cat("Usage: fibrtraq-validation-suite.R <manifest-file-or-dir> [--output-csv <path>] [--preflight-only] [--overwrite] [--no-maps] [--no-report] [--stop-on-error] [--verbose]\n")
  quit(status = 1)
}

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || "--help" %in% args || "-h" %in% args) usage()

manifest_arg <- args[1]
opts <- args[-1]

value_after <- function(flag) {
  idx <- which(opts == flag)
  if (!length(idx)) return(NULL)
  if (idx[1] == length(opts)) stop(flag, " requires a value.")
  opts[idx[1] + 1L]
}

output_csv <- value_after("--output-csv")
if (is.null(output_csv)) {
  matched <- grep("^--output-csv=", opts, value = TRUE)
  if (length(matched)) output_csv <- sub("^--output-csv=", "", matched[1])
}

suppressPackageStartupMessages(library(FibrTraQ))

summary <- ftq_run_validation_suite(
  manifests = manifest_arg,
  output_csv = output_csv,
  preflight_only = "--preflight-only" %in% opts,
  overwrite = "--overwrite" %in% opts,
  write_maps = if ("--no-maps" %in% opts) FALSE else NULL,
  write_report = if ("--no-report" %in% opts) FALSE else NULL,
  continue_on_error = !("--stop-on-error" %in% opts),
  verbose = "--verbose" %in% opts
)

cat("FibrTraQ validation suite complete\n")
cat("  manifests: ", nrow(summary), "\n", sep = "")
if (!is.null(output_csv)) cat("  summary:   ", output_csv, "\n", sep = "")
print(summary[, intersect(c("prefix", "run_status", "preflight_passed",
                            "validation_passed", "focal_voxel_streamline_associations",
                            "leave_no_target_associations", "error"),
                          names(summary))],
      row.names = FALSE)
