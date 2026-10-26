#!/usr/bin/env Rscript

usage <- function() {
  cat(
    "Usage: fibrtraq-exploratory-report.R <result-csv> <template-image> <output-dir> [--background <image>] [--prefix <prefix>] [--targets <comma-list>] [--assignment-mode <mode>] [--search-radius-mm <radius>] [--slices <comma-list>]\n",
    "\n",
    "Builds a target report and five figures in <output-dir>/figures/ and\n",
    "a primary-target lookup CSV in <output-dir>.\n",
    "Keep the matching run-metadata JSON beside the CSV to recover its method.\n",
    "Otherwise provide the original --assignment-mode and, for local_search,\n",
    "--search-radius-mm. Slice numbers are one-based axial array indices.\n",
    sep = ""
  )
  quit(status = 1L)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || args[1] %in% c("-h", "--help")) usage()

result_csv <- args[[1]]
template_image <- args[[2]]
output_dir <- args[[3]]
background_image <- NULL
prefix <- "fibrtraq-exploratory"
targets <- NULL
assignment_mode <- NULL
search_radius_mm <- NULL
report_slices <- NULL

i <- 4L
while (i <= length(args)) {
  arg <- args[[i]]
  if (arg == "--background") {
    i <- i + 1L
    if (i > length(args)) usage()
    background_image <- args[[i]]
  } else if (arg == "--prefix") {
    i <- i + 1L
    if (i > length(args)) usage()
    prefix <- args[[i]]
  } else if (arg == "--targets") {
    i <- i + 1L
    if (i > length(args)) usage()
    targets <- trimws(strsplit(args[[i]], ",", fixed = TRUE)[[1]])
    targets <- targets[nzchar(targets)]
  } else if (arg == "--assignment-mode") {
    i <- i + 1L
    if (i > length(args)) usage()
    assignment_mode <- args[[i]]
  } else if (arg == "--search-radius-mm") {
    i <- i + 1L
    if (i > length(args)) usage()
    search_radius_mm <- suppressWarnings(as.numeric(args[[i]]))
    if (length(search_radius_mm) != 1L || !is.finite(search_radius_mm) || search_radius_mm <= 0) {
      stop("--search-radius-mm must be a finite, positive number.")
    }
  } else if (arg == "--slices") {
    i <- i + 1L
    if (i > length(args)) usage()
    report_slices <- suppressWarnings(as.numeric(strsplit(args[[i]], ",", fixed = TRUE)[[1L]]))
    if (!length(report_slices) || any(!is.finite(report_slices)) ||
        any(report_slices < 1 | report_slices != floor(report_slices))) {
      stop("--slices must contain positive, whole-number axial slice indices.")
    }
  } else {
    stop("Unknown argument: ", arg)
  }
  i <- i + 1L
}

if (!requireNamespace("FibrTraQ", quietly = TRUE)) {
  stop("FibrTraQ must be installed or loaded before running this script.")
}

report <- FibrTraQ::ftq_exploratory_report(
  result_csv = result_csv,
  output_dir = output_dir,
  template_image = template_image,
  prefix = prefix,
  target_labels = targets,
  background_image = background_image,
  assignment_mode = assignment_mode,
  search_radius_mm = search_radius_mm,
  report_slices = report_slices,
  overwrite = TRUE
)

cat("FibrTraQ target report complete\n")
cat("  html:            ", report$html, "\n", sep = "")
cat("  target lookup:   ", report$summary$outputs$lut_csv, "\n", sep = "")
cat("  figures:         ", length(report$figures), "\n", sep = "")
