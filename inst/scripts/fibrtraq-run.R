#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || args[1] %in% c("-h", "--help")) {
  cat("Usage: fibrtraq-run.R <manifest.yml|manifest.json> [--overwrite] [--no-maps] [--no-report] [--preflight-only] [--verbose]\n")
  quit(status = if (length(args) < 1L) 1L else 0L)
}

manifest <- args[1]
overwrite <- "--overwrite" %in% args
write_maps <- !("--no-maps" %in% args)
write_report <- !("--no-report" %in% args)
preflight_only <- "--preflight-only" %in% args
verbose <- "--verbose" %in% args

suppressPackageStartupMessages(library(FibrTraQ))

if (preflight_only) {
  preflight <- ftq_preflight(
    manifest = manifest,
    overwrite = overwrite,
    write_maps = write_maps
  )
  cat("FibrTraQ preflight ", if (isTRUE(preflight$passed)) "PASS" else "FAIL", "\n", sep = "")
  print(preflight$checks, row.names = FALSE)
  quit(status = if (isTRUE(preflight$passed)) 0L else 2L)
}

out <- ftq_run(
  manifest = manifest,
  overwrite = overwrite,
  write_maps = write_maps,
  write_report = write_report,
  verbose = verbose
)

cat("FibrTraQ complete\n")
cat("  result:     ", out$outputs$result_csv, "\n", sep = "")
cat("  rsc:        ", out$outputs$rsc_csv, "\n", sep = "")
cat("  validation: ", out$outputs$validation_json, "\n", sep = "")
if (!is.null(out$report$html)) {
  cat("  report:     ", out$report$html, "\n", sep = "")
}
cat("  validation status: ", if (isTRUE(out$validation$passed)) "PASS" else "FAIL", "\n", sep = "")

quit(status = if (isTRUE(out$validation$passed)) 0L else 2L)
