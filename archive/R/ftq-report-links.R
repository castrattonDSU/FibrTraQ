# Rebuild only the summary from its original CSV and existing report files.
# Searching within the recorded subject/focal folder supports a relocated bundle
# without rewriting scientific provenance or guessing between duplicate files.
ftq_rebuild_report_index <- function(output_dir) {
  output_dir <- normalizePath(path.expand(output_dir), winslash = "/", mustWork = TRUE)
  index_csv <- file.path(output_dir, "fibrtraq-run-index.csv")
  if (!file.exists(index_csv)) stop("Choose the folder containing fibrtraq-run-index.csv.")
  index <- utils::read.csv(index_csv, stringsAsFactors = FALSE, check.names = FALSE,
                           colClasses = c(subject = "character", focal_region = "character"))
  reports <- c("primary_report", "exploratory_report", "interactive_3d_qc")
  required <- c("subject", "focal_region", "run_status", "validation_passed",
                "biological_qc_status", reports)
  if (!all(required %in% names(index))) stop("The run index is missing required columns.")
  # Unified reports already link the copied bundle relatively. Open the
  # existing report instead of creating a second HTML file. Its marker
  # distinguishes the cohort entry point from historical summaries.
  current <- file.path(output_dir, "fibrtraq-report-index.html")
  if (nrow(index) == 1L) {
    run_dir <- file.path(output_dir, ftq_sanitize_label(index$subject[[1L]]),
                         ftq_sanitize_label(index$focal_region[[1L]]))
    current <- c(current, list.files(run_dir, pattern = "_report[.]html$", full.names = TRUE))
  }
  existing <- current[file.exists(current)]
  unified <- vapply(existing, function(file) {
    header <- readLines(file, n = 4L, warn = FALSE, encoding = "UTF-8")
    any(grepl('<meta name="fibrtraq-report-layout" content="unified-v1">', header, fixed = TRUE))
  }, logical(1))
  if (sum(unified) == 1L) return(unname(existing[unified]))
  if (sum(unified) > 1L) stop("More than one unified report was found in this output directory.")
  path <- file.path(output_dir, "fibrtraq-report-index-portable.html")
  ftq_check_overwrite(path, FALSE)
  for (i in seq_len(nrow(index))) {
    run_dir <- file.path(output_dir, ftq_sanitize_label(index$subject[i]),
                         ftq_sanitize_label(index$focal_region[i]))
    candidates <- list.files(run_dir, pattern = "[.]html$", recursive = TRUE,
                             full.names = TRUE, ignore.case = TRUE)
    for (column in reports) {
      old <- index[[column]][i]
      if (is.na(old) || !nzchar(old)) next
      name <- basename(ftq_report_normalize_path(old))
      matches <- candidates[basename(candidates) == name]
      if (length(matches) != 1L) {
        stop("Expected one existing report named ", name, " under ", run_dir,
             "; found ", length(matches), ". Keep the complete results folder together.")
      }
      index[[column]][i] <- matches[[1L]]
    }
  }
  ftq_write_regional_study_index(index, path, output_dir, "regional-study")
  path
}
