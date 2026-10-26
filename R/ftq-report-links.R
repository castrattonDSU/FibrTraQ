# Filesystem paths are raw text, not pre-encoded URLs. Encode reserved filename
# characters (including literal percent escapes), while retaining URL separators.
ftq_url_path <- function(path) {
  gsub("%2F", "/", utils::URLencode(enc2utf8(path), reserved = TRUE,
                                  repeated = TRUE), fixed = TRUE)
}

# png() treats filenames as sprintf templates; reports write one named image.
ftq_png_filename <- function(path) gsub("%", "%%", path, fixed = TRUE)

ftq_report_normalize_path <- function(path) {
  path <- path.expand(path)
  windows_path <- grepl("^[A-Za-z]:[/\\\\]|^\\\\\\\\|^//", path)
  if (windows_path) path <- chartr("\\", "/", path)
  # Explicit Windows paths must also be interpretable when rebuilding a report
  # copied to another OS. Do not turn a POSIX filename's literal '\' into '/'.
  if (!windows_path || identical(.Platform$OS.type, "windows")) {
    if (!grepl("^/|^[A-Za-z]:/", path)) path <- file.path(getwd(), path)
    path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  }
  # On Windows normalizePath can retain the leading UNC backslashes even with
  # winslash = "/". Normalize its result as well as the original input.
  if (windows_path || identical(.Platform$OS.type, "windows")) {
    path <- chartr("\\", "/", path)
  }
  path
}

ftq_report_path_parts <- function(path) {
  path <- ftq_report_normalize_path(path)
  if (grepl("^[A-Za-z]:/", path)) {
    anchor <- substr(path, 1L, 2L)
    rest <- substring(path, 4L)
    windows <- TRUE
  } else if (startsWith(path, "//")) {
    pieces <- strsplit(substring(path, 3L), "/", fixed = TRUE)[[1L]]
    anchor <- paste0("//", paste(utils::head(pieces, 2L), collapse = "/"))
    rest <- paste(utils::tail(pieces, max(0L, length(pieces) - 2L)), collapse = "/")
    windows <- TRUE
  } else {
    anchor <- "/"
    rest <- sub("^/", "", path)
    windows <- FALSE
  }
  parts <- character()
  for (part in strsplit(rest, "/", fixed = TRUE)[[1L]]) {
    if (!nzchar(part) || part == ".") next
    if (part == "..") {
      if (length(parts)) parts <- parts[-length(parts)]
    } else {
      parts <- c(parts, part)
    }
  }
  list(anchor = anchor, parts = parts, windows = windows, path = path)
}

ftq_report_relative_path <- function(path, root) {
  target <- ftq_report_path_parts(path)
  base <- ftq_report_path_parts(root)
  fold <- function(x) if (target$windows && base$windows) tolower(x) else x
  if (!identical(fold(target$anchor), fold(base$anchor))) return(target$path)
  common <- 0L
  while (common < min(length(target$parts), length(base$parts)) &&
         identical(fold(target$parts[common + 1L]), fold(base$parts[common + 1L]))) {
    common <- common + 1L
  }
  relative <- c(rep("..", length(base$parts) - common),
                target$parts[seq_along(target$parts) > common])
  if (!length(relative)) "." else paste(relative, collapse = "/")
}

ftq_file_url <- function(path) {
  path <- ftq_report_normalize_path(path)
  if (grepl("^[A-Za-z]:/", path)) {
    paste0("file:///", substr(path, 1L, 2L), ftq_url_path(substring(path, 3L)))
  } else if (startsWith(path, "//")) {
    paste0("file:", ftq_url_path(path))
  } else {
    paste0("file://", ftq_url_path(path))
  }
}

ftq_report_href <- function(path, root) {
  relative <- ftq_report_relative_path(path, root)
  if (grepl("^/|^[A-Za-z]:/", relative)) ftq_file_url(relative) else ftq_url_path(relative)
}

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
