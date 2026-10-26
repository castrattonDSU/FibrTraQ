#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L && nzchar(args[[1]])) args[[1]] else getwd()
root <- normalizePath(root, mustWork = TRUE)

excluded_dirs <- c(
  ".git",
  ".Rproj.user",
  "FibrTraQ.Rcheck",
  "..Rcheck",
  "manuscript",
  ".local-r-lib",
  ".venv-validation"
)
excluded_ext <- c(".tar.gz", ".o", ".so", ".rds", ".png", ".jpg", ".jpeg", ".gif", ".nii", ".gz", ".trk")

paths <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
paths <- paths[file.exists(paths) & !dir.exists(paths)]
rel <- substring(paths, nchar(root) + 2L)

is_excluded_dir <- vapply(strsplit(rel, .Platform$file.sep, fixed = TRUE), function(parts) {
  any(parts %in% excluded_dirs) || any(grepl("\\.Rcheck$", parts))
}, logical(1))
is_excluded_ext <- vapply(tolower(rel), function(path) {
  any(endsWith(path, excluded_ext))
}, logical(1))
paths <- paths[!is_excluded_dir & !is_excluded_ext]
rel <- rel[!is_excluded_dir & !is_excluded_ext]

patterns <- c(
  mac_user_path = paste0("/", "Users", "/"),
  windows_user_path = paste0("C:[/\\\\]", "Users", "[/\\\\]"),
  hidden_app_attachment = "\\.[^/\\\\]+[/\\\\]attachments"
)

hits <- data.frame(file = character(), pattern = character(), line = integer(), text = character())

for (i in seq_along(paths)) {
  lines <- tryCatch(readLines(paths[[i]], warn = FALSE), error = function(e) character())
  if (!length(lines)) next
  for (nm in names(patterns)) {
    idx <- grep(patterns[[nm]], lines, perl = TRUE, useBytes = TRUE)
    if (length(idx)) {
      hits <- rbind(
        hits,
        data.frame(
          file = rel[[i]],
          pattern = nm,
          line = idx,
          text = trimws(lines[idx]),
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

if (nrow(hits)) {
  write.csv(hits, stdout(), row.names = FALSE)
  stop("Path audit failed: personal absolute paths were found.", call. = FALSE)
}

message("Path audit passed: no personal absolute paths found in package-visible source files.")
