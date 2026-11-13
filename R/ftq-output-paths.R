# Use the same report paths for writing and preflight. In particular, a short
# output root alone is insufficient: the prefix occurs in both directory and
# filename, and target labels also appear in ranked table filenames.
ftq_report_output_paths <- function(output_dir, prefix, target_labels = NULL,
                                    write_endpoint = FALSE,
                                    write_local_search = FALSE) {
  figure <- function(suffix) file.path(output_dir, "report-figures", paste0(prefix, suffix))
  table <- function(suffix) file.path(output_dir, "tables", paste0(prefix, suffix))
  figures <- list(
    fate = figure("_fig_fate-summary.png"),
    rsc = figure("_fig_rsc-summary.png"),
    tracks = figure("_fig_tracks-distribution.png"),
    spatial_tracks = figure("_fig_spatial-tracks.png"),
    spatial_targets = figure("_fig_spatial-targets.png"),
    spatial_leave_no_target = figure("_fig_spatial-leave-no-target.png"),
    local_fraction = figure("_fig_local-fraction.png")
  )
  tables <- list(
    all_voxels = table("_all-voxels.csv"),
    all_voxel_columns = table("_all-voxels-columns.csv"),
    local_fraction = table("_local-fraction.csv"),
    connection_summary = table("_connection-summary.csv"),
    fate_summary = table("_report-fate-summary.csv")
  )
  if (isTRUE(write_endpoint)) {
    figures$endpoint_classification <- figure("_fig_endpoint-classification.png")
    figures$endpoint_comparison <- figure("_fig_endpoint-exact-vs-dilated1.png")
    figures$endpoint_penetration <- figure("_fig_endpoint-penetration.png")
    tables$endpoint_penetration_summary <- table("_endpoint-penetration-summary.csv")
  }
  if (isTRUE(write_local_search)) {
    figures$local_search_distances <- figure("_fig_local-search-end-voxel-distances.png")
    tables$local_search_summary <- table("_local-search-summary.csv")
    tables$local_search_endpoint_status <- table("_local-search-endpoint-status.csv")
  }
  metrics <- c("leave_no_target", target_labels)
  ranked <- stats::setNames(lapply(metrics, function(metric) {
    ftq_ranked_table_path(file.path(output_dir, "tables"), prefix, metric)
  }), metrics)
  list(
    html = file.path(output_dir, paste0(prefix, "_report.html")),
    figures = figures,
    tables = tables,
    ranked = ranked,
    maps = list(
      local_fraction = file.path(output_dir, "maps", paste0(prefix, "_local-fraction.nii.gz")),
      local_fraction_valid = file.path(output_dir, "maps", paste0(prefix, "_local-fraction-valid.nii.gz"))
    ),
    guides = file.path(output_dir, c("OUTPUTS_README.html", "OUTPUTS_README.md"))
  )
}

ftq_ranked_table_path <- function(table_dir, prefix, metric, n = NULL) {
  ranking <- if (is.null(n)) "ranked" else paste0("top", n)
  file.path(table_dir, paste0(prefix, "_", ranking, "_",
                             ftq_sanitize_label(metric), "_voxels.csv"))
}

ftq_absolute_output_path <- function(path) {
  path <- path.expand(path)
  if (!ftq_is_absolute_path(path)) path <- file.path(getwd(), path)
  if (.Platform$OS.type != "windows" && grepl("^([A-Za-z]:[/\\\\]|[/\\\\]{2})", path)) {
    return(gsub("\\", "/", path, fixed = TRUE))
  }
  # normalizePath() alone leaves a nonexistent relative path unresolved. Resolve
  # its existing ancestor without creating any part of the requested output.
  tail <- character()
  ancestor <- path
  while (!file.exists(ancestor) && !identical(dirname(ancestor), ancestor)) {
    tail <- c(basename(ancestor), tail)
    ancestor <- dirname(ancestor)
  }
  ancestor <- normalizePath(ancestor, winslash = "/", mustWork = FALSE)
  if (!length(tail)) return(ancestor)
  # file.path("C:/", "new") can retain the root's slash and add another on
  # Windows. Remove only trailing separators, preserving UNC's leading //.
  paste0(sub("/+$", "", ancestor), "/", paste(tail, collapse = "/"))
}

ftq_run_output_paths <- function(output_dir, prefix, target_labels = NULL,
                                 write_maps = TRUE, write_report = TRUE,
                                 write_exploratory_report = FALSE,
                                 write_qc_viewer = FALSE, write_endpoint = TRUE,
                                 write_local_search = FALSE) {
  output_dir <- ftq_absolute_output_path(output_dir)
  paths <- unlist(ftq_planned_output_paths(
    output_dir, prefix, target_labels, write_maps,
    write_endpoint, write_local_search
  ), use.names = FALSE)
  report_dir <- file.path(output_dir, "reports", prefix)
  if (isTRUE(write_report)) {
    paths <- c(paths, unlist(ftq_report_output_paths(
      report_dir, prefix, target_labels, write_endpoint, write_local_search
    ), use.names = FALSE), file.path(output_dir, "maps",
      paste0(prefix, "_fibrtraq_tractogram-density.nii.gz")))
  }
  if (isTRUE(write_exploratory_report)) {
    suffixes <- c("_fibrtraq_target-report.html", "_fig_exploratory-target-totals.png",
      "_fig_exploratory-primary-target.png", "_fig_exploratory-target-entropy.png",
      "_fig_exploratory-target-diversity.png", "_fig_exploratory-leave-no-target-fraction.png")
    paths <- c(paths, file.path(output_dir, "figures", paste0(prefix, suffixes)),
      file.path(output_dir, paste0(prefix, "_fibrtraq_exploratory-primary-target-lut.csv")))
  }
  if (isTRUE(write_qc_viewer)) {
    paths <- c(paths, file.path(report_dir, "report-figures", c(
      paste0(prefix, "_fig_interactive-3d-qc.html"), "OUTPUTS_README.html", "OUTPUTS_README.md"
    )))
  }
  unique(paths)
}

ftq_regional_study_output_paths <- function(output_dir, subjects, region_labels,
                                            focal_regions, description,
                                            write_maps = TRUE, write_report = TRUE,
                                            write_exploratory_report = TRUE,
                                            write_qc_viewer = TRUE,
                                            write_endpoint = TRUE,
                                            write_local_search = FALSE) {
  output_dir <- ftq_absolute_output_path(output_dir)
  paths <- file.path(output_dir, c("fibrtraq-run-index.csv", "fibrtraq-prepared-inputs.csv",
    "fibrtraq-report-index.html", "OUTPUTS_README.html", "OUTPUTS_README.md"))
  for (subject in subjects) {
    subject_dir <- file.path(output_dir, ftq_sanitize_label(subject))
    # Preparation is conditional on input grids and dimensionality; reserve its
    # possible filenames before opening images or writing any derived inputs.
    paths <- c(paths, file.path(subject_dir, "prepared", paste0(ftq_sanitize_label(subject),
      c("_prepared-reference.nii.gz", "_prepared-reference-volumes.csv",
        "_prepared-background.nii.gz", "_prepared-background-volumes.csv",
        "_regional-union.nii.gz", "_regional-union-inputs.csv", "_regional-union-overlaps.csv"))))
    for (focal in focal_regions) {
      prefix <- ftq_prefix(subject, desc = paste0(ftq_sanitize_label(description), "-",
        ftq_sanitize_label(focal), "-regional-connectivity"))
      paths <- c(paths, ftq_run_output_paths(
        file.path(subject_dir, ftq_sanitize_label(focal)), prefix,
        setdiff(region_labels, focal), write_maps, write_report,
        write_exploratory_report, write_qc_viewer, write_endpoint, write_local_search
      ))
    }
  }
  unique(paths)
}

ftq_manifest_output_paths <- function(manifest) {
  m <- manifest
  paths <- ftq_run_output_paths(
    m$output_dir, m$prefix, m$target_labels, m$write_maps, m$write_report,
    write_endpoint = identical(m$endpoint_analysis, "exact_and_dilated1"),
    write_local_search = identical(m$target_intersection_mode, "local_search")
  )
  cfg <- m$tractogram_context
  has_background <- !is.null(m$background_image) && file.exists(m$background_image)
  needs_context <- identical(m$context_source, "tractogram") ||
    (isTRUE(m$write_report) && identical(m$context_source, "auto") && !has_background)
  if (needs_context && isTRUE(cfg$enabled) && is.null(ftq_context_density_image(cfg))) {
    paths <- c(paths, file.path(
      ftq_absolute_output_path(ftq_value(cfg$output_dir, m$output_dir)), "maps",
      paste0(ftq_value(cfg$prefix, m$prefix), "_fibrtraq_tractogram-density.nii.gz")
    ))
  }
  unique(paths)
}

ftq_windows_path_length <- function(paths) {
  # MAX_PATH is measured in UTF-16 code units, not UTF-8 bytes or Unicode code
  # points. A supplementary character occupies two units on Windows.
  vapply(iconv(enc2utf8(paths), from = "UTF-8", to = "UTF-16LE", toRaw = TRUE),
         function(value) length(value) / 2, numeric(1))
}

ftq_check_windows_output_paths <- function(paths, output_dir,
                                           windows = .Platform$OS.type == "windows") {
  if (!isTRUE(windows) || !length(paths)) return(invisible(paths))
  lengths <- ftq_windows_path_length(paths)
  longest <- which.max(lengths)
  limit <- 259L # MAX_PATH (260), including its terminating NUL.
  if (lengths[[longest]] <= limit) return(invisible(paths))
  output_dir <- ftq_absolute_output_path(output_dir)
  root_budget <- limit - (lengths[[longest]] - ftq_windows_path_length(output_dir))
  below_root <- startsWith(paths[[longest]], paste0(output_dir, "/"))
  recommendation <- if (root_budget > 0 && below_root) {
    paste0("Use a shorter output_dir (at most ", root_budget,
      " UTF-16 code units for these labels and outputs), for example C:/FTQ/run1.")
  } else {
    paste0("Shorten output_dir, any configured tractogram-context output directory, ",
      "or the subject, description, focal, or target labels.")
  }
  stop("Windows output path exceeds the supported limit of 259 UTF-16 code units",
    " (260 including the terminator).\nLongest planned path (", lengths[[longest]],
    "): ", paths[[longest]], "\n", recommendation,
    "\nNo analysis or output files have been started by this call.", call. = FALSE)
}
