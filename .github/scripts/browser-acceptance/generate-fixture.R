#!/usr/bin/env Rscript
# Use an installed candidate and only the package-owned synthetic fixture.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
if (dir.exists(root)) stop("Choose a new browser acceptance output directory.")
dir.create(root, recursive = TRUE)
root <- normalizePath(root, winslash = "/", mustWork = TRUE)
original <- file.path(root, "original report # %20 caf\u00e9")
relocated <- file.path(root, "relocated report # %20 caf\u00e9")

library(FibrTraQ)
example <- ftq_run_example(output_dir = original, backend = "cpp", verbose = FALSE)
stopifnot(isTRUE(example$validation$passed),
          nrow(example$validation$tables) == 7L,
          all(example$validation$tables$passed),
          nrow(example$validation$claims) == 18L,
          all(example$validation$claims$passed),
          all(example$validation$reports$passed))
original <- normalizePath(original, winslash = "/", mustWork = TRUE)
files <- list.files(original, recursive = TRUE, all.files = TRUE,
                    no.. = TRUE, full.names = FALSE)
files <- files[!dir.exists(file.path(original, files))]
dir.create(relocated)
for (relative in files) {
  destination <- file.path(relocated, relative)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(file.path(original, relative), destination, overwrite = FALSE)) {
    stop("Could not copy report asset: ", relative)
  }
}
stopifnot(identical(unname(tools::md5sum(file.path(original, files))),
                    unname(tools::md5sum(file.path(relocated, files)))))
relocated <- normalizePath(relocated, winslash = "/", mustWork = TRUE)
relative_path <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  stopifnot(startsWith(path, paste0(original, "/")))
  substring(path, nchar(original) + 2L)
}
reports <- list(
  index = relative_path(example$report_index),
  results = relative_path(example$primary_report),
  target = relative_path(example$exploratory_report),
  viewer = relative_path(example$interactive_3d_qc)
)
metadata <- example$analysis$index$metadata_json[[1L]]
manifest <- list(
  fixture = "package-owned synthetic known-answer example; not anatomical MRI",
  package_version = as.character(packageVersion("FibrTraQ")),
  source_commit = Sys.getenv("GITHUB_SHA", unset = "local-uncommitted-review"),
  r_version = R.version.string,
  os = unname(Sys.info()[["sysname"]]),
  original_root = original,
  relocated_root = relocated,
  reports = reports,
  metadata = relative_path(metadata),
  html_files = as.list(files[grepl("[.]html$", files, ignore.case = TRUE)]),
  file_count = length(files),
  copy_checksums_identical = TRUE,
  scientific_tables_passed = 7L,
  declared_values_passed = 18L,
  original_outputs_preserved = TRUE
)
jsonlite::write_json(manifest, file.path(root, "browser-fixture.json"),
                     auto_unbox = TRUE, pretty = TRUE)
cat("Browser fixture generated with ", length(files), " files.\n", sep = "")
