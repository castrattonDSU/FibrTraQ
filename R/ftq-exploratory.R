#' Summarize exploratory target connectivity across focal voxels
#'
#' Convert a many-target FibrTraQ per-voxel table into derived exploratory
#' metrics that are easier to interpret biologically: primary target, primary
#' target fraction, target diversity, target entropy, total target associations,
#' and leave-no-target fraction.
#'
#' @param result_csv A FibrTraQ per-voxel CSV path, or a data.frame returned by
#'   [fibrtraq()].
#' @param output_dir Optional directory where summary CSVs and maps are written.
#' @param prefix Output filename prefix. Defaults to `"fibrtraq-exploratory"`.
#' @param target_labels Optional character vector of target columns to include.
#'   If omitted, all non-reserved result columns are used.
#' @param template_image Optional NIfTI/MGH template used to write derived maps.
#' @param write_maps Whether to write derived NIfTI maps. Defaults to TRUE when
#'   `template_image` is supplied.
#' @param overwrite Whether existing exploratory output files may be replaced.
#' @param assignment_mode Optional target-assignment mode represented by
#'   `result_csv`. When `result_csv` is a fresh FibrTraQ data frame, its mode
#'   attribute is used by default. CSV inputs recover the mode from the matching
#'   sibling run-metadata JSON when available. Otherwise the historical
#'   `"full_streamline"` fallback is retained with a warning for CSV inputs and
#'   an explicit unverified-method notice in the report.
#'   Legacy and paper terminology aliases are accepted and normalized.
#' @param search_radius_mm Optional finite, strictly positive scanner-space
#'   search radius in millimeters. A fresh local-search result supplies its
#'   radius attribute or matching run metadata by default. Otherwise local
#'   search retains the historical 2 mm fallback with a warning and report notice.
#'
#' @return A list containing per-voxel exploratory data, target totals,
#'   assignment provenance, and paths to any written outputs.
#' @export
ftq_exploratory_summary <- function(result_csv,
                                    output_dir = NULL,
                                    prefix = "fibrtraq-exploratory",
                                    target_labels = NULL,
                                    template_image = NULL,
                                    write_maps = !is.null(template_image),
                                    overwrite = FALSE,
                                    assignment_mode = NULL,
                                    search_radius_mm = NULL) {
  result <- ftq_read_table(result_csv)
  method <- ftq_exploratory_method(result_csv, result, assignment_mode, search_radius_mm)
  target_intersection_mode <- method$mode
  local_search_radius_mm <- method$radius
  targets <- ftq_result_target_columns(result, target_labels)
  target_mat <- as.matrix(result[targets])
  storage.mode(target_mat) <- "numeric"
  target_mat[!is.finite(target_mat)] <- 0

  tracks <- as.numeric(result$tracks)
  local <- as.numeric(result$local)
  leave_no_target <- as.numeric(result$leave_no_target)
  leave_with_target <- tracks - local - leave_no_target
  total_target_hits <- rowSums(target_mat)
  target_diversity <- rowSums(target_mat > 0)

  primary_index <- max.col(target_mat, ties.method = "first")
  primary_hits <- target_mat[cbind(seq_len(nrow(target_mat)), primary_index)]
  has_target <- total_target_hits > 0 & primary_hits > 0
  primary_label <- rep("none", nrow(target_mat))
  primary_label[has_target] <- targets[primary_index[has_target]]
  primary_index[!has_target] <- 0L
  primary_fraction <- ifelse(total_target_hits > 0, primary_hits / total_target_hits, 0)

  entropy_bits <- ftq_target_entropy(target_mat, total_target_hits)
  entropy_fraction <- ifelse(
    target_diversity > 1,
    entropy_bits / log2(target_diversity),
    0
  )
  entropy_fraction[!is.finite(entropy_fraction)] <- 0
  leave_no_target_fraction <- ifelse(tracks > 0, leave_no_target / tracks, 0)
  leave_with_target_fraction <- ifelse(tracks > 0, leave_with_target / tracks, 0)
  target_hits_per_assigned_track <- ifelse(tracks > 0, total_target_hits / tracks, 0)

  coords <- as.data.frame(ftq_coord_matrix(result$fNode_Coords))
  exploratory <- data.frame(
    fNode_Coords = result$fNode_Coords,
    x = coords$x,
    y = coords$y,
    z = coords$z,
    tracks = tracks,
    local = local,
    leave_no_target = leave_no_target,
    leave_with_target = leave_with_target,
    total_target_hits = total_target_hits,
    target_hits_per_assigned_track = target_hits_per_assigned_track,
    target_diversity = target_diversity,
    primary_target = primary_label,
    primary_target_index = primary_index,
    primary_target_hits = primary_hits,
    primary_target_fraction = primary_fraction,
    target_entropy_bits = entropy_bits,
    target_entropy_fraction = entropy_fraction,
    leave_no_target_fraction = leave_no_target_fraction,
    leave_with_target_fraction = leave_with_target_fraction,
    stringsAsFactors = FALSE
  )

  total_assigned <- sum(tracks)
  target_associations <- colSums(target_mat)
  target_totals <- data.frame(
    target = targets,
    focal_voxel_streamline_target_associations = target_associations,
    percent_of_focal_voxel_streamline_associations = if (total_assigned > 0) 100 * target_associations / total_assigned else 0,
    focal_voxels_with_association = colSums(target_mat > 0),
    stringsAsFactors = FALSE
  )
  target_totals <- target_totals[order(-target_totals$focal_voxel_streamline_target_associations, target_totals$target), , drop = FALSE]
  rownames(target_totals) <- NULL

  provenance <- data.frame(
    target_intersection_mode = target_intersection_mode,
    local_search_radius_mm = if (identical(
      target_intersection_mode,
      "local_search"
    )) {
      local_search_radius_mm
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
  # Keep the on-disk provenance schema stable; display context accompanies the
  # returned table as an attribute rather than becoming repeated CSV columns.
  attr(provenance, "report_context") <- method
  outputs <- list(
    summary_csv = NULL,
    target_totals_csv = NULL,
    lut_csv = NULL,
    provenance_csv = NULL,
    maps = list()
  )
  if (!is.null(output_dir)) {
    output_dir <- normalizePath(path.expand(output_dir), mustWork = FALSE)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    outputs$summary_csv <- file.path(output_dir, paste0(prefix, "_fibrtraq_exploratory-per-voxel.csv"))
    outputs$target_totals_csv <- file.path(output_dir, paste0(prefix, "_fibrtraq_exploratory-target-totals.csv"))
    outputs$lut_csv <- file.path(output_dir, paste0(prefix, "_fibrtraq_exploratory-primary-target-lut.csv"))
    outputs$provenance_csv <- file.path(
      output_dir,
      paste0(prefix, "_fibrtraq_exploratory-provenance.csv")
    )
    ftq_check_overwrite(
      c(
        outputs$summary_csv,
        outputs$target_totals_csv,
        outputs$lut_csv,
        outputs$provenance_csv
      ),
      overwrite
    )
    utils::write.csv(exploratory, outputs$summary_csv, row.names = FALSE)
    utils::write.csv(target_totals, outputs$target_totals_csv, row.names = FALSE)
    utils::write.csv(ftq_exploratory_lut(targets), outputs$lut_csv, row.names = FALSE)
    utils::write.csv(provenance, outputs$provenance_csv, row.names = FALSE)
  }

  if (isTRUE(write_maps)) {
    if (is.null(template_image)) {
      stop("template_image is required when write_maps = TRUE.")
    }
    map_dir <- if (!is.null(output_dir)) file.path(output_dir, "maps") else tempdir()
    outputs$maps <- ftq_write_exploratory_maps(
      exploratory = exploratory,
      template_image = template_image,
      map_dir = map_dir,
      prefix = prefix,
      overwrite = overwrite
    )
  }

  list(
    data = exploratory,
    target_totals = target_totals,
    targets = targets,
    provenance = provenance,
    outputs = outputs
  )
}

#' Build a FibrTraQ target report
#'
#' Write a primary-target lookup CSV, five figures, and a compact HTML report
#' from a many-target FibrTraQ result table. Derived metrics are retained in the
#' returned summary. Use [ftq_exploratory_summary()] for explicit research CSV
#' and NIfTI exports.
#'
#' @param result_csv A FibrTraQ per-voxel CSV path, or data.frame.
#' @param output_dir Focal-result directory. The lookup CSV is written here;
#'   the target report and its images are written in its `figures` subdirectory.
#' @param template_image NIfTI/MGH template defining the spatial figure grid.
#' @param prefix Output filename prefix.
#' @param target_labels Optional target columns to include.
#' @param background_image Optional anatomical image for spatial figures.
#' @param tractogram_context Optional tractogram-derived context output or
#'   density image path used when `background_image` is absent.
#' @param overwrite Whether existing exploratory outputs may be replaced.
#' @param assignment_mode Optional target-assignment mode represented by
#'   `result_csv`. When `result_csv` is a fresh FibrTraQ data frame, its mode
#'   attribute is used by default. CSV inputs recover the mode from the matching
#'   sibling run-metadata JSON when available. Otherwise the historical
#'   `"full_streamline"` fallback is retained with a warning for CSV inputs and
#'   an explicit unverified-method notice in the report.
#'   Legacy and paper terminology aliases are accepted and normalized.
#' @param search_radius_mm Optional finite, strictly positive scanner-space
#'   search radius in millimeters. A fresh local-search result supplies its
#'   radius attribute or matching run metadata by default. Otherwise local
#'   search retains the historical 2 mm fallback with a warning and report notice.
#' @param report_slices Optional one-based axial slice indices for the spatial
#'   previews. NULL retains automatic slice selection; this does not alter maps.
#'
#' @return A list with summary outputs, figure paths, and HTML path.
#' @export
ftq_exploratory_report <- function(result_csv,
                                   output_dir,
                                   template_image,
                                   prefix = "fibrtraq-exploratory",
                                   target_labels = NULL,
                                   background_image = NULL,
                                   overwrite = TRUE,
                                   tractogram_context = NULL,
                                   assignment_mode = NULL,
                                   search_radius_mm = NULL,
                                   report_slices = NULL) {
  output_dir <- normalizePath(path.expand(output_dir), mustWork = FALSE)
  fig_dir <- file.path(output_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  summary <- ftq_exploratory_summary(
    result_csv = result_csv,
    output_dir = NULL,
    prefix = prefix,
    target_labels = target_labels,
    template_image = template_image,
    write_maps = FALSE,
    overwrite = overwrite,
    assignment_mode = assignment_mode,
    search_radius_mm = search_radius_mm
  )

  summary$outputs$lut_csv <- file.path(output_dir,
    paste0(prefix, "_fibrtraq_exploratory-primary-target-lut.csv"))
  html <- file.path(fig_dir, paste0(prefix, "_fibrtraq_target-report.html"))
  figures <- list(
    target_totals = file.path(fig_dir, paste0(prefix, "_fig_exploratory-target-totals.png")),
    primary_target = file.path(fig_dir, paste0(prefix, "_fig_exploratory-primary-target.png")),
    entropy = file.path(fig_dir, paste0(prefix, "_fig_exploratory-target-entropy.png")),
    diversity = file.path(fig_dir, paste0(prefix, "_fig_exploratory-target-diversity.png")),
    no_target = file.path(fig_dir, paste0(prefix, "_fig_exploratory-leave-no-target-fraction.png"))
  )
  ftq_check_overwrite(c(summary$outputs$lut_csv, html, unlist(figures)), overwrite)
  utils::write.csv(ftq_exploratory_lut(summary$targets), summary$outputs$lut_csv,
    row.names = FALSE)

  dims <- dim(RNifti::readNifti(template_image, internal = FALSE))[1:3]
  # Share the research exporter grid construction without writing intermediate
  # maps. Build one array at a time so report memory does not scale with the
  # number of research metrics.
  map_array <- function(metric) ftq_exploratory_array(summary$data, dims, metric)
  background <- NULL
  background_label <- NULL
  if (!is.null(background_image)) {
    background <- ftq_read_background(background_image, dims)
    if (!is.null(background)) background_label <- "anatomical image"
  }
  if (is.null(background)) {
    context_image <- ftq_context_density_image(tractogram_context)
    background <- ftq_read_context_background(context_image, dims)
    if (!is.null(background)) background_label <- "Tractogram-derived density context"
  }

  ftq_plot_exploratory_target_totals(summary$target_totals, figures$target_totals)
  ftq_plot_categorical_montage(
    arr = map_array("primary_target_index"),
    lut = ftq_exploratory_lut(summary$targets),
    path = figures$primary_target,
    title = "Primary Target Per Focal Voxel",
    background = background,
    context = TRUE,
    background_label = background_label,
    slices = report_slices
  )
  ftq_plot_montage(
    map_array("target_entropy_fraction"),
    figures$entropy,
    "Normalized Target Entropy Per Focal Voxel",
    background = background,
    palette = "Viridis",
    legend_label = "Entropy",
    context = TRUE,
    background_label = background_label,
    slices = report_slices
  )
  ftq_plot_montage(
    map_array("target_diversity"),
    figures$diversity,
    "Number of Assigned Targets Per Focal Voxel",
    background = background,
    palette = "Viridis",
    legend_label = "Targets",
    context = TRUE,
    background_label = background_label,
    slices = report_slices
  )
  ftq_plot_montage(
    map_array("leave_no_target_fraction"),
    figures$no_target,
    "Leaving the Configured Local ROI Without Target",
    background = background,
    palette = "Inferno",
    legend_label = "Fraction",
    context = TRUE,
    background_label = background_label,
    slices = report_slices
  )
  figures <- figures[file.exists(unlist(figures, use.names = FALSE))]

  ftq_write_exploratory_html_report(html, summary, figures, background_label)
  list(
    summary = summary,
    figures = figures,
    html = html
  )
}

ftq_exploratory_method <- function(input, result, assignment_mode, search_radius_mm) {
  scalar <- function(x) {
    if (is.null(x) || length(x) != 1L || is.na(x) || !nzchar(as.character(x))) return(NULL)
    if (is.character(x) && toupper(trimws(x)) %in% c("", "NA", "NULL")) return(NULL)
    x
  }
  metadata <- NULL
  metadata_path <- NULL
  if (is.character(input) && length(input) == 1L &&
      grepl("_fibrtraq_per-voxel-hits\\.csv$", input)) {
    candidate <- sub("_fibrtraq_per-voxel-hits\\.csv$", "_fibrtraq_run-metadata.json", input)
    if (file.exists(candidate)) {
      metadata <- tryCatch(jsonlite::read_json(candidate, simplifyVector = FALSE),
        error = function(e) stop("Cannot read companion run metadata: ", conditionMessage(e)))
      recorded_file <- scalar(metadata$files$result_csv)
      if (!is.null(recorded_file) && !identical(
        basename(ftq_report_normalize_path(recorded_file)), basename(input)
      )) stop("Companion run metadata names a different result CSV.")
      metadata_path <- candidate
    }
  }
  from_result <- scalar(attr(result, "fibrtraq_target_intersection_mode", exact = TRUE))
  from_summary <- scalar(metadata$summary$target_intersection_mode)
  from_manifest <- scalar(metadata$manifest$target_intersection_mode)
  recorded_modes <- Filter(Negate(is.null), list(from_result, from_summary, from_manifest))
  recorded_modes <- vapply(recorded_modes, ftq_normalize_target_intersection_mode, character(1))
  if (length(unique(recorded_modes)) > 1L) stop("Conflicting target-assignment modes in result provenance.")
  recorded_mode <- if (length(recorded_modes)) recorded_modes[[1L]] else NULL
  mode <- if (!is.null(assignment_mode)) ftq_normalize_target_intersection_mode(assignment_mode) else recorded_mode
  if (!is.null(assignment_mode) && !is.null(recorded_mode) && !identical(mode, recorded_mode)) {
    stop("assignment_mode conflicts with recorded result provenance (", recorded_mode, ").")
  }
  notes <- character()
  mode_source <- if (!is.null(assignment_mode)) "Explicit argument" else if (!is.null(from_result)) {
    "Result attribute"
  } else if (!is.null(recorded_mode)) "Companion run metadata" else "Historical default; unverified"
  if (is.null(mode)) {
    mode <- "full_streamline"
    notes <- c(notes, paste0("Target-assignment method was not recorded. full_streamline is the historical ",
      "default used for this derived report, not a verified method for these data. ",
      "Supply assignment_mode or retain the matching run-metadata JSON."))
    if (!is.data.frame(input)) warning(notes[[length(notes)]], call. = FALSE)
  }
  radii <- Filter(Negate(is.null), list(
    scalar(attr(result, "fibrtraq_local_search_radius_mm", exact = TRUE)),
    scalar(metadata$summary$local_search_radius_mm),
    scalar(metadata$manifest$local_search_radius_mm),
    scalar(metadata$local_search_audit$search_radius_mm)
  ))
  radius <- search_radius_mm
  radius_source <- if (!is.null(radius)) "Explicit argument" else "Recorded result provenance"
  if (identical(mode, "local_search")) {
    radii <- vapply(radii, ftq_validate_local_search_radius, numeric(1))
    if (length(unique(radii)) > 1L) stop("Conflicting local-search radii in result provenance.")
    if (!is.null(radius) && length(radii) && !isTRUE(all.equal(as.numeric(radius), radii[[1L]]))) {
      stop("search_radius_mm conflicts with recorded result provenance (", radii[[1L]], " mm).")
    }
    if (is.null(radius) && length(radii)) radius <- radii[[1L]]
    if (is.null(radius)) {
      radius_source <- "Historical default; unverified"
      notes <- c(notes, "Local-search radius was not recorded. The historical 2 mm default is unverified; supply search_radius_mm for these data.")
      warning(notes[[length(notes)]], call. = FALSE)
    }
  }
  radius <- ftq_validate_local_search_radius(ftq_value(radius, 2))
  list(mode = mode, radius = radius, mode_source = mode_source,
    radius_source = radius_source, notes = notes, metadata_path = metadata_path,
    subject = scalar(metadata$manifest$subject),
    focal_assignment_mode = ftq_value(
      scalar(attr(result, "fibrtraq_focal_assignment_mode", exact = TRUE)),
      scalar(metadata$summary$focal_assignment_mode), scalar(metadata$manifest$focal_assignment_mode),
      "Not recorded"))
}

ftq_result_target_columns <- function(result, target_labels = NULL) {
  required <- c("fNode_Coords", "tracks", "local", "leave_no_target", "none")
  missing <- setdiff(required, names(result))
  if (length(missing)) {
    stop("result table is missing required column(s): ", paste(missing, collapse = ", "))
  }
  reserved <- c(required, "cNode_Coords")
  targets <- if (is.null(target_labels)) setdiff(names(result), reserved) else as.character(target_labels)
  missing_targets <- setdiff(targets, names(result))
  if (length(missing_targets)) {
    stop("target column(s) are missing from result table: ", paste(missing_targets, collapse = ", "))
  }
  if (!length(targets)) stop("No target columns were found in result table.")
  targets
}

ftq_target_entropy <- function(target_mat, total_target_hits) {
  out <- numeric(nrow(target_mat))
  nz_rows <- which(total_target_hits > 0)
  for (i in nz_rows) {
    p <- target_mat[i, ] / total_target_hits[i]
    p <- p[is.finite(p) & p > 0]
    out[i] <- -sum(p * log2(p))
  }
  out
}

ftq_exploratory_lut <- function(targets) {
  data.frame(
    index = seq_along(targets),
    target = targets,
    color = grDevices::hcl.colors(length(targets), "Dark 3"),
    stringsAsFactors = FALSE
  )
}

ftq_exploratory_array <- function(exploratory, dims, metric) {
  coords <- as.matrix(exploratory[, c("x", "y", "z")])
  storage.mode(coords) <- "integer"
  arr <- array(0, dim = dims)
  arr[cbind(coords[, 1], coords[, 2], coords[, 3])] <- as.numeric(exploratory[[metric]])
  arr
}

ftq_write_exploratory_maps <- function(exploratory,
                                       template_image,
                                       map_dir,
                                       prefix,
                                       overwrite) {
  template <- RNifti::readNifti(template_image, internal = FALSE)
  dims <- dim(template)[1:3]
  dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
  metrics <- c(
    "primary_target_index",
    "primary_target_hits",
    "primary_target_fraction",
    "target_diversity",
    "target_entropy_bits",
    "target_entropy_fraction",
    "total_target_hits",
    "target_hits_per_assigned_track",
    "leave_no_target_fraction",
    "leave_with_target_fraction"
  )
  paths <- list()
  for (metric in metrics) {
    out <- file.path(map_dir, paste0(prefix, "_fibrtraq_exploratory_", ftq_sanitize_label(metric), ".nii.gz"))
    ftq_check_overwrite(out, overwrite)
    arr <- ftq_exploratory_array(exploratory, dims, metric)
    RNifti::writeNifti(arr, out, template = template, datatype = "float")
    paths[[metric]] <- out
  }
  paths
}

ftq_plot_exploratory_target_totals <- function(target_totals, path, n = 25L) {
  top <- utils::head(target_totals, n)
  grDevices::png(ftq_png_filename(path), width = 3000, height = max(1500, 120 * nrow(top)), res = 240)
  old <- graphics::par(mar = c(5, 17, 4, 4), xaxs = "i")
  on.exit({
    graphics::par(old)
    grDevices::dev.off()
  }, add = TRUE)
  top <- top[order(top$focal_voxel_streamline_target_associations), , drop = FALSE]
  cols <- grDevices::hcl.colors(nrow(top), "Dark 3")
  labels <- gsub("_", " ", top$target, fixed = TRUE)
  labels <- vapply(labels, function(x) paste(strwrap(x, width = 34), collapse = "\n"), character(1))
  bp <- graphics::barplot(
    top$focal_voxel_streamline_target_associations,
    names.arg = labels,
    horiz = TRUE,
    las = 1,
    cex.names = 0.58,
    col = cols,
    border = NA,
    axes = FALSE,
    xlab = "Focal-voxel\u2013streamline\u2013target associations",
    main = "Targets Ranked by Voxel\u2013Streamline Associations",
    xlim = c(0, max(1, top$focal_voxel_streamline_target_associations) * 1.15)
  )
  ticks <- pretty(c(0, max(top$focal_voxel_streamline_target_associations)))
  graphics::axis(
    1,
    at = ticks,
    labels = format(ticks, big.mark = ",", scientific = FALSE, trim = TRUE)
  )
  graphics::text(
    top$focal_voxel_streamline_target_associations,
    bp,
    labels = format(top$focal_voxel_streamline_target_associations, big.mark = ","),
    pos = 4,
    cex = 0.65
  )
}

ftq_plot_categorical_montage <- function(arr,
                                         lut,
                                         path,
                                         title,
                                         nslices = 9L,
                                         background = NULL,
                                         context = FALSE,
                                         background_label = NULL,
                                         slices = NULL) {
  a <- as.array(arr)
  slices <- ftq_report_slices(slices, dim(a)[3])
  if (is.null(slices)) {
    zsum <- apply(a > 0, 3, sum)
    slices <- which(zsum > 0)
    if (!length(slices)) slices <- seq_len(dim(a)[3])
    if (length(slices) > nslices) {
      slices <- unique(round(seq(min(slices), max(slices), length.out = nslices)))
    }
  }
  crop <- ftq_crop_ranges(list(a), margin = 3L)
  use_context <- isTRUE(context) && !is.null(background)
  slices_per_row <- min(3L, length(slices))
  nr <- ceiling(length(slices) / slices_per_row)
  panels_per_row <- if (use_context) slices_per_row * 2L else slices_per_row
  png_width <- if (use_context) 1120L * slices_per_row + 720L else 860L * slices_per_row + 720L
  grDevices::png(ftq_png_filename(path), width = png_width, height = max(900L, 620L * nr), res = 240)
  layout_mat <- cbind(
    matrix(seq_len(nr * panels_per_row), nrow = nr, byrow = TRUE),
    rep(nr * panels_per_row + 1L, nr)
  )
  panel_widths <- if (use_context) {
    rep(c(0.78, 1.08), slices_per_row)
  } else {
    rep(1, slices_per_row)
  }
  graphics::layout(layout_mat, widths = c(panel_widths, 0.85), heights = rep(1, nr))
  old <- graphics::par(mar = c(1.5, 1.5, 3, 1), oma = c(1.4, 0, 3, 0))
  on.exit({
    graphics::par(old)
    graphics::layout(1)
    grDevices::dev.off()
  }, add = TRUE)
  colors <- stats::setNames(lut$color, lut$index)
  for (z in slices) {
    if (use_context) {
      ftq_draw_categorical_slice(
        a[, , z],
        seq_len(dim(a)[1]),
        seq_len(dim(a)[2]),
        colors,
        paste("z =", z, "full slice"),
        background[, , z],
        crop_box = crop
      )
      ftq_draw_categorical_slice(
        a[crop$x, crop$y, z],
        crop$x,
        crop$y,
        colors,
        paste("z =", z, "focal zoom"),
        background[crop$x, crop$y, z]
      )
    } else {
      mat <- a[crop$x, crop$y, z]
      bg_mat <- if (!is.null(background)) background[crop$x, crop$y, z] else NULL
      ftq_draw_categorical_slice(mat, crop$x, crop$y, colors, paste("z =", z), bg_mat)
    }
  }
  filled_panels <- length(slices) * if (use_context) 2L else 1L
  if (filled_panels < nr * panels_per_row) {
    for (i in seq_len(nr * panels_per_row - filled_panels)) graphics::plot.new()
  }
  ftq_draw_categorical_legend(a, lut)
  graphics::mtext(title, outer = TRUE, cex = 1.2, font = 2)
  if (!is.null(background_label) && nzchar(background_label)) {
    graphics::mtext(paste("Background:", background_label), side = 1, outer = TRUE, cex = 0.65)
  }
}

ftq_draw_categorical_slice <- function(map_mat, x, y, colors, main, bg_mat = NULL, crop_box = NULL) {
  if (!is.null(bg_mat)) {
    graphics::image(x, y, ftq_normalize_background(bg_mat), col = grDevices::gray.colors(128, 0, 0.95),
                    zlim = c(0, 1), axes = FALSE, asp = 1, main = main, xlab = "", ylab = "")
  } else {
    graphics::image(x, y, matrix(0, nrow = nrow(map_mat), ncol = ncol(map_mat)),
                    col = "#f8fafc", zlim = c(0, 1), axes = FALSE, asp = 1,
                    main = main, xlab = "", ylab = "")
  }
  overlay <- map_mat
  overlay[!is.finite(overlay) | overlay <= 0] <- NA
  if (any(!is.na(overlay))) {
    graphics::image(x, y, overlay, col = unname(colors), zlim = c(0.5, length(colors) + 0.5),
                    axes = FALSE, asp = 1, add = TRUE, xlab = "", ylab = "")
  }
  if (!is.null(crop_box)) {
    graphics::rect(min(crop_box$x), min(crop_box$y), max(crop_box$x), max(crop_box$y),
                   border = "#1f9bd1", lwd = 2)
  }
  graphics::box()
}

ftq_draw_categorical_legend <- function(arr, lut, n = 14L) {
  used <- sort(unique(as.integer(as.array(arr)[as.array(arr) > 0])))
  used <- used[used %in% lut$index]
  lut <- lut[match(used, lut$index), , drop = FALSE]
  if (nrow(lut) > n) {
    lut <- lut[seq_len(n), , drop = FALSE]
    footer <- paste0("+ ", length(used) - n, " more")
  } else {
    footer <- ""
  }
  old <- graphics::par(mar = c(1.5, 0.5, 3, 0.5))
  on.exit(graphics::par(old), add = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, max(1, nrow(lut) + 1)))
  if (nrow(lut)) {
    y <- rev(seq_len(nrow(lut)))
    graphics::rect(0.05, y - 0.35, 0.18, y + 0.35, col = lut$color, border = NA)
    labels <- paste0(lut$index, ": ", lut$target)
    graphics::text(0.22, y, labels = labels, adj = c(0, 0.5), cex = 0.58)
  }
  if (nzchar(footer)) graphics::text(0.05, 0.25, footer, adj = c(0, 0.5), cex = 0.58)
  graphics::title("Primary target", cex.main = 0.75)
}

ftq_write_exploratory_html_report <- function(path, summary, figures, background_label = NULL) {
  esc <- function(x) {
    x <- as.character(x)
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x <- gsub('"', "&quot;", x, fixed = TRUE)
    x <- gsub("'", "&#39;", x, fixed = TRUE)
    x
  }
  href <- function(p) esc(ftq_report_href(p, dirname(path)))
  link <- function(p, label) {
    if (is.null(p) || !length(p) || !file.exists(p)) return("")
    paste0("<li><a href='", href(p), "' download>", esc(label), "</a></li>")
  }
  table_html <- function(df, n = 20L) {
    df <- utils::head(df, n)
    header <- paste0("<tr>", paste0("<th scope='col'>", esc(names(df)), "</th>", collapse = ""), "</tr>")
    rows <- if (nrow(df)) apply(df, 1, function(r) paste0("<tr>", paste0("<td>", esc(r), "</td>", collapse = ""), "</tr>")) else ""
    paste0("<div class='table-wrap'><table><thead>", header, "</thead><tbody>", paste(rows, collapse = "\n"), "</tbody></table></div>")
  }
  labels <- c(
    target_totals = "Targets ranked by voxel\u2013streamline associations",
    primary_target = "Primary target at each focal voxel",
    entropy = "Normalized target entropy",
    diversity = "Number of assigned targets at each focal voxel",
    no_target = "Fraction leaving the configured local ROI without an assigned target"
  )
  captions <- c(
    target_totals = "Up to 25 targets are shown; the target association table includes every target. Each bar sums target associations across focal voxels. One streamline may contribute at several voxels; these are not unique streamline counts.",
    primary_target = "The target with the largest association count at each voxel. Equal maxima select the first target in the input target-column order. Voxels without target associations have no primary target.",
    entropy = "Entropy of the positive target-association proportions, divided by log2 of the number of positive targets. The stored value is zero for one positive target and for no target associations; distinguish these using target coverage.",
    diversity = "Number of target columns with a positive association count at each focal voxel. This does not count streamlines.",
    no_target = "leave_no_target / tracks at voxels with assigned streamlines. Local refers to the configured local ROI, which can differ from the focal ROI. A stored zero at a zero-assignment voxel is unevaluable."
  )
  fig_tag <- function(p, name) {
    label <- if (name %in% names(labels)) labels[[name]] else gsub("_", " ", name, fixed = TRUE)
    caption <- if (name %in% names(captions)) captions[[name]] else ""
    paste0("<figure><h3>", esc(label), "</h3><a href='", href(p), "'><img src='", href(p),
      "' alt='", esc(label), "'></a><figcaption>", esc(caption), " <a href='", href(p),
      "' download>Download image (PNG)</a></figcaption></figure>")
  }
  value <- function(x) if (!length(x) || !is.finite(x)) "Not available" else format(round(x, 4), trim = TRUE, big.mark = ",")
  median_value <- function(x) if (length(x)) value(stats::median(x)) else "Not available (no eligible voxels)"
  data <- summary$data
  n <- nrow(data)
  target_positive <- data$total_target_hits > 0
  count_percent <- function(k) paste0(format(k, big.mark = ","), " / ", format(n, big.mark = ","),
    if (n > 0) paste0(" (", format(round(100 * k / n, 1), trim = TRUE), "%)") else " (not available)")
  overview <- data.frame(
    Measure = c("All focal voxels", "Voxels with assigned streamlines", "Zero-assignment focal voxels",
      "Voxels with target associations", "Median target diversity \u2014 all focal voxels",
      "Median primary-target fraction \u2014 all focal voxels", "Median primary-target fraction \u2014 target-positive voxels"),
    Value = c(value(n), count_percent(sum(data$tracks > 0)), count_percent(sum(data$tracks == 0)),
      count_percent(sum(target_positive)), median_value(data$target_diversity),
      median_value(data$primary_target_fraction), median_value(data$primary_target_fraction[target_positive])),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  method <- attr(summary$provenance, "report_context", exact = TRUE)
  mode <- as.character(summary$provenance$target_intersection_mode[[1L]])
  radius_note <- if (identical(mode, "local_search")) paste0(
    "<br><b>Local search radius (mm):</b> ", esc(summary$provenance$local_search_radius_mm[[1L]])) else ""
  warning_note <- if (length(method$notes)) paste0("<div class='warning'><b>Method provenance requires attention</b><ul>",
    paste0("<li>", esc(method$notes), "</li>", collapse = ""), "</ul></div>") else ""
  context_note <- if (identical(background_label, "Tractogram-derived density context")) {
    "<div class='meta'><b>Spatial QC note:</b> Figures use Tractogram-derived density context because no anatomical image was used. This context is not anatomical MRI and is not brain, tissue, or lesion volume.</div>"
  } else ""
  totals <- summary$target_totals
  names(totals) <- c("Target", "Voxel\u2013streamline\u2013target associations", "Associations / assigned voxel\u2013streamline pairs (%)", "Focal voxels with target associations")
  # Reuse an existing guide in this bundle rather than duplicating it next to
  # every report. Standalone reports remain useful without an external guide.
  guide <- NULL
  guide_dir <- dirname(normalizePath(path, mustWork = FALSE))
  repeat {
    candidates <- list.files(guide_dir,
      pattern = "^OUTPUTS_README(-[a-f0-9]{12})?[.]html$", full.names = TRUE)
    for (candidate in candidates) {
      markdown <- sub("[.]html$", ".md", basename(candidate))
      if (identical(readLines(candidate, warn = FALSE), ftq_output_guide_html(markdown))) {
        guide <- candidate
        break
      }
    }
    if (!is.null(guide)) break
    parent <- dirname(guide_dir)
    if (identical(parent, guide_dir)) break
    guide_dir <- parent
  }
  guide_nav <- if (!is.null(guide)) paste0("<a href='", href(guide), "'>Output definitions</a>") else ""
  guide_download <- if (!is.null(guide)) paste0("<li><a href='", href(guide),
    "'>Output definitions and interpretation guide (HTML)</a></li>") else ""
  downloads <- paste0(
    link(summary$outputs$summary_csv, "Complete exploratory metrics \u2014 every focal voxel (CSV)"),
    link(summary$outputs$lut_csv, "Primary-target index and color lookup (CSV)"))
  supporting <- paste0(
    link(summary$outputs$target_totals_csv, "Target association totals (CSV)"),
    link(summary$outputs$provenance_csv, "Exploratory assignment provenance (CSV)"),
    link(method$metadata_path, "Original run metadata (JSON)"))
  map_labels <- c(primary_target_index = "Primary-target index", primary_target_hits = "Primary-target association count",
    primary_target_fraction = "Primary-target association fraction", target_diversity = "Number of positive targets",
    target_entropy_bits = "Target entropy (bits)", target_entropy_fraction = "Normalized target entropy",
    total_target_hits = "Total target associations", target_hits_per_assigned_track = "Target associations per assigned streamline",
    leave_no_target_fraction = "Leaving the local ROI without a target \u2014 fraction", leave_with_target_fraction = "Leaving the local ROI with a target \u2014 fraction")
  maps <- paste(vapply(names(summary$outputs$maps), function(name) {
    label <- if (name %in% names(map_labels)) map_labels[[name]] else name
    link(summary$outputs$maps[[name]], paste0(label, " (NIfTI)"))
  }, character(1)), collapse = "")
  method_rows <- data.frame(
    Item = c("Derived metrics", "Target-assignment evidence", "Focal-voxel allocation", "Counting unit", "Local region", "Coordinates", "Primary-target ties"),
    Definition = c(if (!is.null(summary$outputs$summary_csv)) basename(summary$outputs$summary_csv) else "Calculated from the supplied per-voxel assignments and retained in the returned report summary",
      ftq_value(method$mode_source, "Recorded in supplied summary"),
      ftq_value(method$focal_assignment_mode, "Not recorded"),
      "Focal-voxel\u2013streamline\u2013target associations; one streamline can contribute at several focal voxels",
      "Configured local ROI; it need not equal the focal ROI",
      "x, y, z are one-based focal-image array indices, not scanner millimetres",
      paste0("First target among equal maxima in this input order: ", paste(summary$targets, collapse = ", "))),
    stringsAsFactors = FALSE
  )
  if (identical(mode, "local_search")) method_rows <- rbind(method_rows,
    data.frame(Item = "Radius evidence", Definition = ftq_value(method$radius_source, "Recorded in supplied summary")))
  title <- paste0("FibrTraQ target report", if (!is.null(method$subject)) paste0(" \u2014 ", method$subject) else "")
  css <- "<style>body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;max-width:1120px;margin:32px auto;padding:0 20px;line-height:1.5;color:#222}.table-wrap{overflow-x:auto}table{border-collapse:collapse;margin:12px 0 24px;width:100%}th,td{border:1px solid #ddd;padding:8px 12px;text-align:left;vertical-align:top}th{background:#f3f4f6}img{max-width:100%;height:auto;border:1px solid #eee}figure{margin:20px 0 32px}.meta,.warning{padding:12px 16px;margin:16px 0;border-left:4px solid #64748b;background:#f8fafc}.warning{border-color:#ad6800;background:#fff8e6}details{margin:20px 0;padding:12px;border:1px solid #ddd}summary{cursor:pointer;font-weight:600}figcaption{color:#374151}a{color:#075ba8}nav a{margin-right:18px}</style>"
  html <- paste0(
    "<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'><title>",
    esc(title), "</title>", css, "</head><body><h1>", esc(title), "</h1>",
    "<nav><a href='#coverage'>Coverage</a><a href='#figures'>Figures</a><a href='#downloads'>Downloads</a>", guide_nav, "</nav>",
    "<div class='meta'><b>Purpose:</b> inspect how target associations vary across focal voxels.",
    "<br><b>Target assignment mode:</b> ", esc(mode), radius_note, "</div>", warning_note, context_note,
    "<h2 id='coverage'>Coverage and interpretation</h2>", table_html(overview, nrow(overview)),
    "<p>Coverage percentages use all focal voxels as the denominator. Zero assigned streamlines does not always mean no streamline passes through a voxel: first_intersection assigns each streamline to its first focal voxel. Zero target associations can occur even when streamlines are assigned.</p>",
    "<p>Primary-target fractions stored as zero at voxels without target associations are included in the all-voxel median. The target-positive median excludes those voxels and describes the balance among assigned targets.</p>",
    "<h2 id='figures'>Spatial and target-association figures</h2>",
    paste(vapply(names(figures), function(name) fig_tag(figures[[name]], name), character(1)), collapse = "\n"),
    "<details><summary>Target association table</summary><p>All targets, ranked by associations. Percentages divide target associations by all assigned focal-voxel\u2013streamline pairs. Targets can overlap, so these percentages do not partition 100%.</p>",
    table_html(totals, nrow(totals)), "</details>",
    "<h2 id='downloads'>Lookup and definitions</h2><ul>", downloads, supporting,
    guide_download, "</ul>",
    if (nzchar(maps)) paste0("<details><summary>Research maps</summary><p>Use the primary-target lookup to interpret indices and the method definitions to interpret zero values and counting units.</p><ul>", maps, "</ul></details>") else "",
    "<details><summary>Methods and provenance</summary>", table_html(method_rows, nrow(method_rows)),
    "<p>This report summarizes saved assignments; it does not rerun tractography or independently validate anatomical connectivity.</p></details>",
    "</body></html>"
  )
  writeLines(html, path)
}
