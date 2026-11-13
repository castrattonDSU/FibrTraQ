#' Read a FibrTraQ run manifest
#'
#' @param manifest Path to a JSON/YAML manifest, or an R list.
#' @return A normalized manifest list.
#' @export
ftq_read_manifest <- function(manifest) {
  if (is.list(manifest)) {
    x <- manifest
  } else if (is.character(manifest) && length(manifest) == 1L) {
    path <- path.expand(manifest)
    if (!file.exists(path)) stop("Manifest file does not exist: ", manifest)
    ext <- tolower(basename(path))
    x <- if (grepl("\\.json$", ext)) {
      jsonlite::fromJSON(path, simplifyVector = FALSE)
    } else if (grepl("\\.(ya?ml)$", ext)) {
      yaml::read_yaml(path)
    } else {
      stop("Manifest must be a .json, .yml, or .yaml file.")
    }
    x$manifest_path <- normalizePath(path, mustWork = TRUE)
  } else {
    stop("manifest must be a path or list.")
  }

  if (is.null(x$dataset_root)) x$dataset_root <- "."
  x$dataset_root <- normalizePath(path.expand(x$dataset_root), mustWork = FALSE)
  if (is.null(x$subject)) stop("Manifest requires 'subject'.")
  if (is.null(x$description) && !is.null(x$desc)) x$description <- x$desc
  if (is.null(x$description)) stop("Manifest requires 'description' or 'desc'.")
  if (is.null(x$fnode)) stop("Manifest requires 'fnode'.")
  if (is.null(x$tractogram)) stop("Manifest requires 'tractogram'.")
  if (is.null(x$targets) || length(x$targets) == 0L) stop("Manifest requires one or more targets.")

  x$fnode <- ftq_manifest_path(x$fnode, x$dataset_root)
  local_path <- x[["local"]]
  if (!is.null(local_path) && !is.na(local_path) && nzchar(local_path)) {
    x[["local"]] <- ftq_manifest_path(local_path, x$dataset_root)
  } else {
    x[["local"]] <- NULL
  }
  x$tractogram <- ftq_manifest_path(x$tractogram, x$dataset_root)
  if (!is.null(x$ref_image) && !is.na(x$ref_image) && nzchar(x$ref_image)) {
    x$ref_image <- ftq_manifest_path(x$ref_image, x$dataset_root)
  }

  report <- x$report
  if (is.null(report)) report <- list()
  background_image <- ftq_value(report$background_image, x$background_image, NULL)
  if (!is.null(background_image) && !is.na(background_image) && nzchar(background_image)) {
    x$background_image <- ftq_manifest_path(background_image, x$dataset_root)
  } else {
    x$background_image <- NULL
  }
  x$context_source <- tolower(as.character(ftq_value(report$context_source, x$context_source, "auto"))[1])
  context_alias <- c(background = "anatomy", anatomical = "anatomy", map_only = "none", maponly = "none")
  if (x$context_source %in% names(context_alias)) x$context_source <- unname(context_alias[[x$context_source]])
  x$context_source <- match.arg(x$context_source, c("auto", "anatomy", "tractogram", "none"))
  x$tractogram_context <- ftq_parse_tractogram_context(
    ftq_value(report$tractogram_context, x$tractogram_context, NULL),
    x$dataset_root
  )

  parsed_targets <- ftq_parse_targets(x$targets, x$dataset_root)
  x$target_labels <- parsed_targets$labels
  x$target_paths <- parsed_targets$paths
  parsed_endpoint_contexts <- ftq_parse_endpoint_contexts(
    x$endpoint_contexts,
    x$dataset_root
  )
  x$endpoint_context_labels <- parsed_endpoint_contexts$labels
  x$endpoint_context_paths <- parsed_endpoint_contexts$paths

  coord <- x$coordinate
  if (is.null(coord)) coord <- list()
  x$track_space <- ftq_value(coord$track_space, x$track_space, "auto")
  x$discretization <- ftq_value(coord$discretization, x$discretization, "any")
  axis_raw <- ftq_value(coord$axis_flip, x$axis_flip, NULL)
  if (is.null(axis_raw) || (is.character(axis_raw) && identical(tolower(axis_raw), "auto"))) {
    x$axis_flip <- NULL
  } else {
    x$axis_flip <- as.logical(unlist(axis_raw))[1:3]
    if (length(x$axis_flip) != 3L || anyNA(x$axis_flip)) {
      stop("axis_flip must contain three logical values, or be omitted for auto-detection.")
    }
  }

  analysis <- x$analysis
  if (is.null(analysis)) analysis <- list()
  x$focal_assignment_mode <- as.character(
    ftq_value(analysis$focal_assignment_mode, x$focal_assignment_mode, "all_focal_voxels")
  )[1]
  x$focal_assignment_mode <- match.arg(
    x$focal_assignment_mode,
    c("all_focal_voxels", "first_intersection")
  )
  x$streamline_voxelization_mode <- as.character(
    ftq_value(analysis$streamline_voxelization_mode, x$streamline_voxelization_mode, "sample_points")
  )[1]
  x$streamline_voxelization_mode <- match.arg(
    x$streamline_voxelization_mode,
    c("sample_points", "segment_interpolation")
  )
  x$target_intersection_mode <- ftq_normalize_target_intersection_mode(
    ftq_value(
      analysis$target_intersection_mode,
      x$target_intersection_mode,
      "full_streamline"
    )
  )
  x$local_search_radius_mm <- ftq_validate_local_search_radius(
    ftq_value(
      analysis$local_search_radius_mm,
      x$local_search_radius_mm,
      2
    )
  )
  endpoint_analysis_explicit <- if ("endpoint_analysis_explicit" %in% names(x)) {
    isTRUE(x$endpoint_analysis_explicit)
  } else {
    !is.null(analysis$endpoint_analysis) || !is.null(x$endpoint_analysis)
  }
  x$endpoint_analysis <- as.character(
    ftq_value(analysis$endpoint_analysis, x$endpoint_analysis, "exact_and_dilated1")
  )[1]
  x$endpoint_analysis <- match.arg(
    x$endpoint_analysis,
    c("exact_and_dilated1", "none")
  )
  x$endpoint_analysis_explicit <- endpoint_analysis_explicit
  if (identical(x$endpoint_analysis, "none") && length(x$endpoint_context_paths)) {
    stop("endpoint_contexts require analysis.endpoint_analysis = 'exact_and_dilated1'.")
  }
  x$reverse_streamlines <- isTRUE(ftq_value(
    analysis$reverse_streamlines,
    x$reverse_streamlines,
    FALSE
  ))

  execution <- x$execution
  if (is.null(execution)) execution <- list()
  x$track_reader <- tolower(as.character(ftq_value(execution$track_reader, x$track_reader, "auto"))[1])
  x$track_chunk_size <- as.integer(ftq_value(execution$track_chunk_size, x$track_chunk_size, 50000L))[1]

  outputs <- x$outputs
  if (is.null(outputs)) outputs <- list()
  x$derivative_name <- ftq_value(outputs$derivative_name, x$derivative_name, "fibrtraq-v0.3.0")
  x$output_dir <- ftq_value(outputs$output_dir, x$output_dir, NULL)
  if (is.null(x$output_dir)) {
    x$output_dir <- file.path(x$dataset_root, "derivatives", x$derivative_name)
  } else {
    x$output_dir <- ftq_manifest_path(x$output_dir, x$dataset_root)
  }
  x$write_maps <- isTRUE(ftq_value(outputs$write_maps, x$write_maps, TRUE))
  x$write_report <- isTRUE(ftq_value(outputs$write_report, x$write_report, TRUE))
  x$overwrite <- isTRUE(ftq_value(outputs$overwrite, x$overwrite, FALSE))
  x$expected <- ftq_parse_expected(ftq_value(outputs$expected, x$expected, NULL), x$dataset_root)
  x$prefix <- ftq_prefix(x$subject, x$session, x$description)

  x
}

#' Run FibrTraQ from a manifest
#'
#' Endpoint analysis is enabled by default through
#' `analysis.endpoint_analysis: exact_and_dilated1`. The returned `endpoint`
#' member and three dedicated endpoint CSVs are additive; manifests selecting
#' `none` retain the legacy output contract.
#'
#' @param manifest Path to JSON/YAML manifest, or a manifest list.
#' @param output_dir Optional output directory override.
#' @param overwrite Whether existing core output files may be overwritten.
#' @param write_maps Whether to write NIfTI map outputs.
#' @param write_report Whether to build an HTML/PNG report.
#' @param verbose Emit progress messages.
#' @param backend Counting backend passed to [FibrTraQ()].
#' @param progress Emit main tractogram pass progress messages.
#' @param progress_every Number of streamlines between progress messages.
#' @param profile Attach backend stage timing information to the returned object.
#' @param report_slices Optional distinct, positive one-based axial slice indices
#'   for report figures. NULL retains automatic slice selection.
#' @param report_error_action Whether a report-generation error should stop the
#'   call (the default) or return completed scientific outputs with report status.
#' @return A list containing paths, result tables, validation, report outputs,
#'   and an additive `endpoint` member when endpoint analysis is enabled.
#' @export
ftq_run <- function(manifest,
                    output_dir = NULL,
                    overwrite = NULL,
                    write_maps = NULL,
                    write_report = NULL,
                    verbose = FALSE,
                    backend = c("r", "cpp", "auto"),
                    progress = verbose,
                    progress_every = 100000L,
                    profile = FALSE,
                    report_slices = NULL,
                    report_error_action = c("stop", "return")) {
  backend <- match.arg(backend)
  report_error_action <- match.arg(report_error_action)
  m <- ftq_read_manifest(manifest)
  if (!is.null(output_dir)) m$output_dir <- output_dir
  if (!is.null(overwrite)) m$overwrite <- isTRUE(overwrite)
  if (!is.null(write_maps)) m$write_maps <- isTRUE(write_maps)
  if (!is.null(write_report)) m$write_report <- isTRUE(write_report)

  ftq_check_windows_output_paths(ftq_manifest_output_paths(m), m$output_dir)
  preflight <- ftq_preflight(m, output_dir = m$output_dir, overwrite = m$overwrite)
  if (!isTRUE(preflight$passed)) {
    stop("FibrTraQ preflight failed:\n", paste(ftq_preflight_failures(preflight), collapse = "\n"))
  }
  dir.create(m$output_dir, recursive = TRUE, showWarnings = FALSE)

  fnode <- ftq_read_volume_data(m$fnode)
  if (is.null(m$ref_image)) {
    # The array reader deliberately leaves coordinate mapping unchanged. Keep
    # focal-file spacing for the endpoint audit instead of the unit spacing
    # RNifti would synthesize from an array without its image header.
    attr(fnode, "fibrtraq_voxel_size_mm") <- ftq_endpoint_voxel_size(
      fnode, ref_image = m$fnode
    )
  }
  report_slices <- ftq_report_slices(ftq_value(report_slices, m$report$slices), dim(fnode)[3L])
  local <- if (is.null(m[["local"]])) {
    NULL
  } else {
    ftq_read_volume_data(m[["local"]])
  }
  cnodes <- lapply(m$target_paths, ftq_read_volume_data)
  names(cnodes) <- m$target_labels
  endpoint_contexts <- lapply(m$endpoint_context_paths, function(path) {
    mask <- ftq_read_volume_data(path)
    attr(mask, "fibrtraq_source_path") <- path
    attr(mask, "fibrtraq_checksum_md5") <- unname(tools::md5sum(path))
    mask
  })
  names(endpoint_contexts) <- m$endpoint_context_labels
  trk <- ftq_read_tractogram(
    m$tractogram,
    reader = m$track_reader,
    chunk_size = m$track_chunk_size
  )
  m$tractogram_format <- ftq_track_format(trk)
  compute_ref_image <- m$ref_image
  if (is.null(compute_ref_image) &&
      (identical(m$target_intersection_mode, "local_search") ||
       identical(m$tractogram_format, "tck"))) {
    compute_ref_image <- m$fnode
  }
  if (identical(m$tractogram_format, "tck")) {
    coordinate_settings <- ftq_tck_coordinate_settings(
      m$track_space, m$axis_flip, compute_ref_image
    )
    m$track_space <- coordinate_settings$track_space
    m$axis_flip <- coordinate_settings$axis_flip
    m$coordinate_mapping <- coordinate_settings$metadata
  }

  compute_args <- list(
    fNode = fnode,
    localROI = local,
    cNode_list = cnodes,
    cNode_labels = m$target_labels,
    all_tracks = trk,
    discretization = m$discretization,
    focal_assignment_mode = m$focal_assignment_mode,
    streamline_voxelization_mode = m$streamline_voxelization_mode,
    assignment_mode = m$target_intersection_mode,
    search_radius_mm = m$local_search_radius_mm,
    reverse_streamlines = m$reverse_streamlines,
    verbose = verbose,
    track_space = m$track_space,
    ref_image = compute_ref_image,
    backend = backend,
    progress = progress,
    progress_every = progress_every,
    profile = profile,
    endpoint_analysis = m$endpoint_analysis,
    endpoint_context_list = endpoint_contexts,
    endpoint_context_labels = m$endpoint_context_labels
  )
  if (!is.null(m$axis_flip)) compute_args$axis_flip <- m$axis_flip
  result <- do.call(fibrtraq, compute_args)
  rsc <- rsc_summary(result, percent = TRUE, warn_legacy = FALSE)
  context <- ftq_run_tractogram_context(m, result, overwrite = m$overwrite)
  outputs <- ftq_write_outputs(
    result = result,
    rsc = rsc,
    template_image = m$fnode,
    output_dir = m$output_dir,
    prefix = m$prefix,
    target_labels = m$target_labels,
    manifest = m,
    preflight = preflight,
    context = context,
    write_maps = m$write_maps,
    overwrite = m$overwrite
  )
  validation <- ftq_validate_run(
    result_csv = outputs$result_csv,
    rsc_csv = outputs$rsc_csv,
    map_paths = outputs$map_paths,
    expected = m$expected,
    validation_path = outputs$validation_json,
    roi_summary_csv = outputs$roi_summary_csv,
    endpoint_summary = outputs$endpoint_summary_csv,
    endpoint_penetration = outputs$endpoint_penetration_csv,
    endpoint_context = outputs$endpoint_context_csv,
    local_search_assignments = outputs$local_search_assignments_csv
  )
  validation$focal_assignment_mode <- m$focal_assignment_mode
  validation$streamline_voxelization_mode <- m$streamline_voxelization_mode
  validation$target_intersection_mode <- m$target_intersection_mode
  validation$local_search_radius_mm <- m$local_search_radius_mm
  validation$local_search_spatial_units <- ftq_value(
    attr(result, "fibrtraq_local_search_spatial_units"),
    NA_character_
  )
  validation$local_search_spatial_units_scale_to_mm <- ftq_value(
    attr(result, "fibrtraq_local_search_spatial_units_scale_to_mm"),
    NA_real_
  )
  validation$local_search_spatial_units_interpretation <- ftq_value(
    attr(result, "fibrtraq_local_search_spatial_units_interpretation"),
    NA_character_
  )
  validation$endpoint_analysis <- m$endpoint_analysis
  validation$reverse_streamlines <- m$reverse_streamlines
  jsonlite::write_json(validation, outputs$validation_json, auto_unbox = TRUE, pretty = TRUE)
  report <- NULL
  report_status <- "not_requested"
  report_error <- NULL
  if (isTRUE(m$write_report)) {
    report <- tryCatch(ftq_report(
      result_csv = outputs$result_csv,
      rsc_csv = outputs$rsc_csv,
      manifest = m,
      output_dir = file.path(m$output_dir, "reports", m$prefix),
      map_paths = outputs$map_paths,
      validation = validation,
      preflight = preflight,
      tractogram_context = context,
      roi_summary = outputs$roi_summary_csv,
      endpoint_summary = outputs$endpoint_summary_csv,
      endpoint_penetration = outputs$endpoint_penetration_csv,
      endpoint_context = outputs$endpoint_context_csv,
      local_search_assignments = outputs$local_search_assignments_csv,
      overwrite = TRUE,
      report_slices = report_slices
    ), error = function(e) {
      report_error <<- conditionMessage(e)
      NULL
    })
    report_status <- if (is.null(report_error)) "completed" else "failed"
    if (!is.null(report_error) && identical(report_error_action, "stop")) {
      stop("Scientific outputs and validation were saved in ", m$output_dir,
           ", but report generation failed: ", report_error, call. = FALSE)
    }
  }

  list(
    manifest = m,
    preflight = preflight,
    result = result,
    rsc = rsc,
    roi_summary = attr(result, "fibrtraq_roi_summary"),
    endpoint = if (identical(m$endpoint_analysis, "exact_and_dilated1")) {
      ftq_endpoint_tables(result)
    } else {
      NULL
    },
    local_search = if (identical(
      m$target_intersection_mode,
      "local_search"
    )) {
      ftq_local_search_tables(result)
    } else {
      NULL
    },
    outputs = outputs,
    validation = validation,
    report = report,
    report_status = report_status,
    report_error = report_error,
    backend = if (!is.null(attr(result, "fibrtraq_backend"))) attr(result, "fibrtraq_backend") else ftq_resolve_backend(backend),
    profile = attr(result, "fibrtraq_profile"),
    package_identity = ftq_package_identity()
  )
}

#' Preflight-check a FibrTraQ manifest
#'
#' @param manifest Path to JSON/YAML manifest, or a manifest list.
#' @param output_dir Optional output directory override.
#' @param overwrite Whether existing outputs are allowed.
#' @param write_maps Whether map outputs are planned.
#' @return A structured preflight result with blocking errors, warnings, and recommendations.
#' @export
ftq_preflight <- function(manifest,
                          output_dir = NULL,
                          overwrite = NULL,
                          write_maps = NULL) {
  m <- ftq_read_manifest(manifest)
  if (!is.null(output_dir)) m$output_dir <- output_dir
  if (!is.null(overwrite)) m$overwrite <- isTRUE(overwrite)
  if (!is.null(write_maps)) m$write_maps <- isTRUE(write_maps)

  checks <- list()
  add <- function(check, severity, passed, detail = "") {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = check,
      severity = severity,
      passed = isTRUE(passed),
      detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  }
  add_affine_grid_check <- function(check, grid, severity = "error") {
    affine_available <- isTRUE(grid$affine_available)
    add(
      check,
      if (affine_available) severity else "warning",
      affine_available && isTRUE(grid$affine_match),
      grid$detail
    )
  }

  allowed_track_space <- c("auto", "voxelmm", "vox10", "voxel", "world", "trackvis")
  add("track_space_supported", "error", m$track_space %in% allowed_track_space,
      paste("track_space =", m$track_space))
  is_tck <- grepl("\\.tck$", tolower(m$tractogram))
  effective_ref_image <- m$ref_image
  if (is.null(effective_ref_image) &&
      (identical(m$target_intersection_mode, "local_search") || is_tck)) {
    effective_ref_image <- m$fnode
  }
  if (is_tck) {
    settings <- tryCatch(
      ftq_tck_coordinate_settings(m$track_space, m$axis_flip, effective_ref_image),
      error = function(e) e
    )
    add("tck_world_coordinate_mapping", "error", !inherits(settings, "error"),
        if (inherits(settings, "error")) conditionMessage(settings) else
          "Native TCK uses RAS+ world millimeters and the declared reference affine; axis flips are disabled.")
    if (!inherits(settings, "error")) {
      m$track_space <- settings$track_space
      m$axis_flip <- settings$axis_flip
    }
    spatial_paths <- c(
      list(fnode = m$fnode),
      if (!is.null(m[["local"]])) list(local = m[["local"]]),
      stats::setNames(as.list(m$target_paths), paste0("target:", m$target_labels)),
      if (length(m$endpoint_context_paths)) {
        stats::setNames(as.list(m$endpoint_context_paths),
                        paste0("endpoint_context:", m$endpoint_context_labels))
      }
    )
    for (label in names(spatial_paths)) {
      geometry <- tryCatch(
        ftq_tck_validate_reference(spatial_paths[[label]], context = label),
        error = function(e) e
      )
      add(paste0("tck_spatial_geometry:", label), "error",
          !inherits(geometry, "error"),
          if (inherits(geometry, "error")) conditionMessage(geometry) else
            "Declared, unambiguous NIfTI affine with known spatial units.")
    }
  }
  needs_ref_image <- m$track_space %in% c("world", "trackvis")
  add(
    "ref_image_required_for_world_mapping",
    "error",
    !needs_ref_image || !is.null(effective_ref_image),
    paste("track_space =", m$track_space)
  )
  if (identical(m$track_space, "auto")) {
    add(
      "ref_image_recommended_for_auto_mapping",
      "warning",
      !is.null(effective_ref_image),
      "A reference image enables world and TrackVis candidates during automatic coordinate mapping."
    )
  }
  discretization <- as.character(unlist(m$discretization))
  add("discretization_supported", "error",
      length(discretization) > 0L && all(discretization %in% c("any", "round", "floor", "ceil")),
      paste("discretization =", paste(discretization, collapse = ",")))
  add("focal_assignment_mode_supported", "error",
      m$focal_assignment_mode %in% c("all_focal_voxels", "first_intersection"),
      paste("focal_assignment_mode =", m$focal_assignment_mode))
  add("streamline_voxelization_mode_supported", "error",
      m$streamline_voxelization_mode %in% c("sample_points", "segment_interpolation"),
      paste("streamline_voxelization_mode =", m$streamline_voxelization_mode))
  add("target_intersection_mode_supported", "error",
      m$target_intersection_mode %in%
        c("full_streamline", "endpoint_voxels", "local_search"),
      paste("target_intersection_mode =", m$target_intersection_mode))
  add(
    "local_search_radius_mm_positive",
    "error",
    length(m$local_search_radius_mm) == 1L &&
      is.finite(m$local_search_radius_mm) &&
      m$local_search_radius_mm > 0,
    paste("local_search_radius_mm =", m$local_search_radius_mm)
  )
  add("endpoint_analysis_supported", "error",
      m$endpoint_analysis %in% c("exact_and_dilated1", "none"),
      paste("endpoint_analysis =", m$endpoint_analysis))
  add("prefix_path_safe", "error", !grepl("[/\\\\]", m$prefix),
      paste("prefix =", m$prefix))

  reserved <- c("fNode_Coords", "tracks", "local", "leave_no_target", "none", "cNode_Coords")
  add("target_labels_not_reserved", "error", !any(m$target_labels %in% reserved),
      paste(intersect(m$target_labels, reserved), collapse = ","))
  map_names <- paste0("target-", ftq_sanitize_label(m$target_labels))
  add("target_map_names_unique", "error", !anyDuplicated(map_names),
      paste(map_names[duplicated(map_names)], collapse = ","))

  add_path_check <- function(name, path, severity = "error") {
    add(paste0(name, "_exists"), severity, !is.null(path) && file.exists(path),
        if (is.null(path)) "" else path)
  }
  add_path_check("fnode", m$fnode)
  add_path_check("tractogram", m$tractogram)
  add("tractogram_format_supported", "error",
      grepl("\\.(trk|tck)$", tolower(m$tractogram)), m$tractogram)
  if (is_tck && file.exists(m$tractogram)) {
    tck_header <- tryCatch(
      ftq_read_tractogram(m$tractogram, reader = "stream"),
      error = function(e) e
    )
    add("tck_header_readable", "error", !inherits(tck_header, "error"),
        if (inherits(tck_header, "error")) conditionMessage(tck_header) else
          "MRtrix TCK header and payload offset are valid; payload integrity is checked while reading.")
  }
  add("track_reader_supported", "error", m$track_reader %in% c("auto", "memory", "stream"),
      paste("track_reader =", m$track_reader))
  add("track_chunk_size_positive", "error",
      !is.na(m$track_chunk_size) && m$track_chunk_size > 0L,
      paste("track_chunk_size =", m$track_chunk_size))
  if (!is.null(m[["local"]])) add_path_check("local", m[["local"]])
  if (!is.null(m$ref_image)) add_path_check("ref_image", m$ref_image)
  for (i in seq_along(m$target_paths)) {
    add_path_check(paste0("target:", m$target_labels[[i]]), m$target_paths[[i]])
  }
  for (i in seq_along(m$endpoint_context_paths)) {
    add_path_check(
      paste0("endpoint_context:", m$endpoint_context_labels[[i]]),
      m$endpoint_context_paths[[i]]
    )
  }
  if (!is.null(m$background_image)) {
    add_path_check("background_image", m$background_image, severity = "warning")
  } else {
    can_context <- !identical(m$context_source, "none") &&
      isTRUE(m$tractogram_context$enabled) &&
      grepl("\\.(trk|tck)$", tolower(m$tractogram))
    if (can_context) {
      add("spatial_context_fallback_available", "recommendation", TRUE,
          "No anatomical background supplied; reports will use tractogram-derived density context when generated. Anatomy remains recommended for publication figures.")
    } else {
      add("spatial_context_map_only", "recommendation", FALSE,
          "No anatomical background supplied and no supported tractogram-derived context is available; reports will use map-only spatial figures.")
    }
  }
  for (path in ftq_expected_paths(m$expected)) {
    add("expected_output_exists", "error", file.exists(path), path)
  }

  out_paths <- ftq_planned_output_paths(
    m$output_dir,
    m$prefix,
    m$target_labels,
    m$write_maps,
    write_endpoint = identical(m$endpoint_analysis, "exact_and_dilated1"),
    write_local_search = identical(
      m$target_intersection_mode,
      "local_search"
    )
  )
  existing_core <- out_paths$core[file.exists(unlist(out_paths$core, use.names = FALSE))]
  map_files <- unlist(out_paths$maps, use.names = FALSE)
  if (is.null(map_files)) map_files <- character(0)
  existing_maps <- out_paths$maps[file.exists(map_files)]
  add("core_outputs_available", "error", isTRUE(m$overwrite) || length(existing_core) == 0L,
      paste(unlist(existing_core, use.names = FALSE), collapse = "\n"))
  add("map_outputs_available", "error", isTRUE(m$overwrite) || length(existing_maps) == 0L,
      paste(unlist(existing_maps, use.names = FALSE), collapse = "\n"))

  volume_info <- list()
  fnode_mask <- NULL
  local_mask <- NULL
  target_masks <- list()
  fnode_read <- ftq_read_volume_for_preflight(m$fnode)
  if (!is.null(fnode_read$error)) {
    add("fnode_readable", "error", FALSE, fnode_read$error)
  } else {
    volume_info$fnode <- ftq_volume_info(fnode_read$data)
    if (isTRUE(volume_info$fnode$is_3d)) {
      fnode_mask <- fnode_read$data != 0
    }
    add("fnode_readable", "error", TRUE, m$fnode)
    add(
      "fnode_is_3d",
      "error",
      isTRUE(volume_info$fnode$is_3d),
      paste0(
        "dimensions = ",
        ftq_dim_label(volume_info$fnode$full_dims),
        "; dimensionality = ",
        volume_info$fnode$dimensionality
      )
    )
    add("fnode_nonempty", "error", volume_info$fnode$nonzero_voxels > 0L,
        paste("nonzero_voxels =", volume_info$fnode$nonzero_voxels))
    add("fnode_binary_values", "warning", isTRUE(volume_info$fnode$binary_values),
        paste("min =", volume_info$fnode$min, "max =", volume_info$fnode$max))
  }
  explicit_ref_image_is_used <- !is.null(m$ref_image) && (
    m$track_space %in% c("world", "trackvis", "auto") ||
      identical(m$target_intersection_mode, "local_search")
  )
  if (isTRUE(explicit_ref_image_is_used) &&
      file.exists(m$fnode) && file.exists(m$ref_image)) {
    reference_grid <- ftq_preflight_grid_check(
      m$fnode,
      m$ref_image,
      "fnode",
      "ref_image",
      spatial_units_to_mm = TRUE
    )
    reference_check_prefix <- if (identical(
      m$target_intersection_mode,
      "local_search"
    )) {
      "local_search_ref_image"
    } else {
      "mapping_ref_image"
    }
    add(
      paste0(reference_check_prefix, "_dims_match_fnode"),
      "error",
      isTRUE(reference_grid$dims_match),
      reference_grid$detail
    )
    if (isTRUE(reference_grid$dims_match)) {
      add_affine_grid_check(
        paste0(reference_check_prefix, "_affine_matches_fnode"),
        reference_grid
      )
    }
  }
  if (identical(m$target_intersection_mode, "local_search")) {
    geometry <- if (!isTRUE(volume_info$fnode$is_3d)) {
      simpleError("Local-search geometry requires an exactly 3-D fNode mask.")
    } else tryCatch(
      {
        geometry <- ftq_local_search_geometry(
          fNode = fnode_read$data,
          ref_image = effective_ref_image
        )
        ftq_local_search_voxel_radius(
          m$local_search_radius_mm,
          geometry
        )
        geometry
      },
      error = function(e) e
    )
    geometry_available <- !inherits(geometry, "error")
    add(
      "local_search_geometry_available",
      "error",
      geometry_available,
      if (geometry_available) {
        paste0(
          "voxel_size_mm = ",
          paste(signif(geometry$voxel_size_mm, 8), collapse = ",")
        )
      } else {
        conditionMessage(geometry)
      }
    )
    if (geometry_available) {
      add(
        "local_search_spatial_units_known",
        "warning",
        isTRUE(geometry$spatial_units_known),
        paste0(
          "header_spatial_units = ", geometry$spatial_units,
          "; scale_to_mm = ",
          format(geometry$spatial_units_scale_to_mm, scientific = FALSE),
          "; interpretation = ",
          geometry$spatial_units_interpretation
        )
      )
    }
  }

  fnode_dims <- volume_info$fnode$dims
  if (!is.null(m[["local"]])) {
    local_read <- ftq_read_volume_for_preflight(m[["local"]])
    if (!is.null(local_read$error)) {
      add("local_readable", "error", FALSE, local_read$error)
    } else {
      volume_info$local <- ftq_volume_info(local_read$data)
      if (isTRUE(volume_info$local$is_3d)) {
        local_mask <- local_read$data != 0
      }
      add("local_readable", "error", TRUE, m[["local"]])
      add(
        "local_is_3d",
        "error",
        isTRUE(volume_info$local$is_3d),
        paste0(
          "dimensions = ",
          ftq_dim_label(volume_info$local$full_dims),
          "; dimensionality = ",
          volume_info$local$dimensionality
        )
      )
      local_grid <- ftq_preflight_grid_check(
        m$fnode,
        m[["local"]],
        "fnode",
        "local",
        spatial_units_to_mm = is_tck
      )
      add("local_dims_match_fnode", "error", isTRUE(local_grid$dims_match), local_grid$detail)
      if (isTRUE(local_grid$dims_match)) {
        add_affine_grid_check("local_affine_matches_fnode", local_grid)
      }
      add("local_nonempty", "warning", volume_info$local$nonzero_voxels > 0L,
          paste("nonzero_voxels =", volume_info$local$nonzero_voxels))
      add("local_binary_values", "warning", isTRUE(volume_info$local$binary_values),
          paste("min =", volume_info$local$min, "max =", volume_info$local$max))
      if (!is.null(fnode_mask) && !is.null(local_mask) &&
          identical(volume_info$local$dims, fnode_dims)) {
        overlap <- sum(fnode_mask & local_mask, na.rm = TRUE)
        add("focal_local_overlap", "recommendation", TRUE,
            paste("overlap_voxels =", overlap,
                  "; local is reported as a non-exclusive flag."))
      }
    }
  } else {
    local_mask <- fnode_mask
    add("local_defaults_to_fnode", "recommendation", TRUE,
        "No local ROI provided; local classification uses the focal ROI and remains a non-exclusive flag.")
    if (!is.null(fnode_mask)) {
      overlap <- sum(fnode_mask, na.rm = TRUE)
      add("focal_local_overlap", "recommendation", TRUE,
          paste("overlap_voxels =", overlap,
                "; local defaults to focal and is reported as a non-exclusive flag."))
    }
  }

  if (length(m$target_paths)) {
    for (i in seq_along(m$target_paths)) {
      label <- m$target_labels[[i]]
      target_read <- ftq_read_volume_for_preflight(m$target_paths[[i]])
      if (!is.null(target_read$error)) {
        add(paste0("target_readable:", label), "error", FALSE, target_read$error)
      } else {
        info <- ftq_volume_info(target_read$data)
        volume_info$targets[[label]] <- info
        add(paste0("target_readable:", label), "error", TRUE, m$target_paths[[i]])
        add(
          paste0("target_is_3d:", label),
          "error",
          isTRUE(info$is_3d),
          paste0(
            "dimensions = ", ftq_dim_label(info$full_dims),
            "; dimensionality = ", info$dimensionality
          )
        )
        target_grid <- ftq_preflight_grid_check(
          m$fnode, m$target_paths[[i]], "fnode", paste0("target:", label),
          spatial_units_to_mm = is_tck
        )
        add(paste0("target_dims_match_fnode:", label), "error", isTRUE(target_grid$dims_match), target_grid$detail)
        if (isTRUE(target_grid$dims_match)) {
          add_affine_grid_check(
            paste0("target_affine_matches_fnode:", label),
            target_grid
          )
        }
        add(paste0("target_nonempty:", label), "error", info$nonzero_voxels > 0L,
            paste("nonzero_voxels =", info$nonzero_voxels))
        add(paste0("target_binary_values:", label), "warning", isTRUE(info$binary_values),
            paste("min =", info$min, "max =", info$max))
        target_mask <- if (isTRUE(info$is_3d)) {
          target_read$data != 0
        } else {
          NULL
        }
        if (!is.null(target_mask)) target_masks[[label]] <- target_mask
        if (!is.null(fnode_mask) && !is.null(target_mask) &&
            identical(info$dims, fnode_dims)) {
          overlap <- sum(fnode_mask & target_mask, na.rm = TRUE)
          add(paste0("focal_target_overlap:", label), "warning", overlap == 0L,
              paste("overlap_voxels =", overlap))
        }
        if (!is.null(local_mask) && !is.null(target_mask) &&
            identical(dim(target_mask), fnode_dims) &&
            identical(dim(local_mask), fnode_dims)) {
          overlap <- sum(local_mask & target_mask, na.rm = TRUE)
          add(paste0("local_target_overlap:", label), "warning", overlap == 0L,
              paste("overlap_voxels =", overlap,
                    "; local is reported as a non-exclusive flag."))
        }
      }
    }
  }

  if (length(target_masks) > 1L) {
    target_labels <- names(target_masks)
    for (a in seq_len(length(target_masks) - 1L)) {
      for (b in seq.int(a + 1L, length(target_masks))) {
        overlap <- sum(target_masks[[a]] & target_masks[[b]], na.rm = TRUE)
        add(paste0("target_target_overlap:", target_labels[[a]], "__", target_labels[[b]]),
            "warning", overlap == 0L, paste("overlap_voxels =", overlap))
      }
    }
  }

  if (length(m$endpoint_context_paths)) {
    for (i in seq_along(m$endpoint_context_paths)) {
      label <- m$endpoint_context_labels[[i]]
      context_read <- ftq_read_volume_for_preflight(m$endpoint_context_paths[[i]])
      if (!is.null(context_read$error)) {
        add(paste0("endpoint_context_readable:", label), "error", FALSE, context_read$error)
        next
      }
      info <- ftq_volume_info(context_read$data)
      volume_info$endpoint_contexts[[label]] <- info
      add(paste0("endpoint_context_readable:", label), "error", TRUE,
          m$endpoint_context_paths[[i]])
      add(
        paste0("endpoint_context_is_3d:", label),
        "error",
        isTRUE(info$is_3d),
        paste0(
          "dimensions = ", ftq_dim_label(info$full_dims),
          "; dimensionality = ", info$dimensionality
        )
      )
      context_grid <- ftq_preflight_grid_check(
        m$fnode,
        m$endpoint_context_paths[[i]],
        "fnode",
        paste0("endpoint_context:", label),
        spatial_units_to_mm = is_tck
      )
      add(paste0("endpoint_context_dims_match_fnode:", label), "error",
          isTRUE(context_grid$dims_match), context_grid$detail)
      if (isTRUE(context_grid$dims_match)) {
        add_affine_grid_check(
          paste0("endpoint_context_affine_matches_fnode:", label),
          context_grid
        )
      }
      expected_voxels <- if (!length(info$full_dims) ||
                              anyNA(info$full_dims)) {
        NA_real_
      } else {
        prod(info$full_dims)
      }
      add(paste0("endpoint_context_finite:", label), "error",
          is.finite(expected_voxels) && info$finite_voxels == expected_voxels,
          paste("finite_voxels =", info$finite_voxels, "expected =", expected_voxels))
      add(paste0("endpoint_context_binary:", label), "error",
          isTRUE(info$binary_values), paste("min =", info$min, "max =", info$max))
      add(paste0("endpoint_context_nonempty:", label), "warning",
          info$nonzero_voxels > 0L, paste("nonzero_voxels =", info$nonzero_voxels))
    }
  }

  if (!is.null(m$background_image) && file.exists(m$background_image)) {
    background_read <- ftq_read_volume_for_preflight(m$background_image)
    if (!is.null(background_read$error)) {
      add("background_image_readable", "warning", FALSE, background_read$error)
    } else {
      volume_info$background_image <- ftq_volume_info(background_read$data)
      add(
        "background_image_is_3d",
        "warning",
        isTRUE(volume_info$background_image$is_3d),
        paste0(
          "dimensions = ",
          ftq_dim_label(volume_info$background_image$full_dims),
          "; dimensionality = ",
          volume_info$background_image$dimensionality
        )
      )
      background_grid <- ftq_preflight_grid_check(m$fnode, m$background_image, "fnode", "background_image")
      add("background_image_dims_match_fnode", "warning",
          isTRUE(background_grid$dims_match), background_grid$detail)
      if (isTRUE(background_grid$dims_match)) {
        add_affine_grid_check(
          "background_image_affine_matches_fnode",
          background_grid,
          severity = "warning"
        )
      }
    }
  }

  checks_df <- if (length(checks)) do.call(rbind, checks) else data.frame(
    check = character(), severity = character(), passed = logical(), detail = character()
  )
  rownames(checks_df) <- NULL
  list(
    passed = !any(!checks_df$passed & checks_df$severity == "error"),
    checks = checks_df,
    manifest = m,
    volume_info = volume_info,
    checked_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
}

#' Write FibrTraQ tables, maps, metadata, and validation placeholders
#'
#' Results carrying endpoint attributes also write dedicated endpoint summary,
#' penetration, and context CSVs. Legacy results omit those files.
#'
#' @param result Per-voxel FibrTraQ result table.
#' @param rsc RSC summary table.
#' @param template_image Reference image path for NIfTI maps.
#' @param output_dir Output directory.
#' @param prefix File prefix.
#' @param target_labels Target column names.
#' @param manifest Optional manifest list for metadata.
#' @param preflight Optional preflight result from [ftq_preflight()].
#' @param context Optional tractogram-derived context output from
#'   [ftq_tractogram_context()].
#' @param write_maps Whether to write NIfTI maps.
#' @param overwrite Whether existing core output files may be overwritten.
#' @return A list of output paths.
#' @export
ftq_write_outputs <- function(result,
                              rsc,
                              template_image,
                              output_dir,
                              prefix,
                              target_labels = NULL,
                              manifest = NULL,
                              preflight = NULL,
                              context = NULL,
                              write_maps = TRUE,
                              overwrite = FALSE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  endpoint_tables <- list(
    summary = attr(result, "fibrtraq_endpoint_summary"),
    penetration = attr(result, "fibrtraq_endpoint_penetration"),
    context = attr(result, "fibrtraq_endpoint_context"),
    metadata = attr(result, "fibrtraq_endpoint_metadata")
  )
  local_search_assignments <- attr(
    result,
    "fibrtraq_local_search_assignments"
  )
  write_local_search <- !is.null(local_search_assignments)
  write_endpoint <- !is.null(endpoint_tables$summary)
  if (write_endpoint && any(vapply(endpoint_tables[1:3], is.null, logical(1)))) {
    stop("Endpoint-enabled results must contain summary, penetration, and context tables.")
  }
  planned <- ftq_planned_output_paths(
    output_dir,
    prefix,
    target_labels,
    write_maps = FALSE,
    write_endpoint = write_endpoint,
    write_local_search = write_local_search
  )
  result_csv <- planned$core$result_csv
  rsc_csv <- planned$core$rsc_csv
  fate_csv <- planned$core$fate_csv
  roi_summary_csv <- planned$core$roi_summary_csv
  metadata_json <- planned$core$metadata_json
  validation_json <- planned$core$validation_json

  ftq_check_overwrite(unlist(planned$core, use.names = FALSE), overwrite)
  utils::write.csv(result, result_csv, row.names = FALSE)
  utils::write.csv(rsc, rsc_csv, row.names = FALSE)
  fate <- ftq_fate_summary(result)
  utils::write.csv(fate, fate_csv, row.names = FALSE)
  roi_summary <- ftq_roi_summary(result)
  utils::write.csv(roi_summary, roi_summary_csv, row.names = FALSE)
  if (write_endpoint) {
    utils::write.csv(
      endpoint_tables$summary,
      planned$core$endpoint_summary_csv,
      row.names = FALSE
    )
    utils::write.csv(
      endpoint_tables$penetration,
      planned$core$endpoint_penetration_csv,
      row.names = FALSE
    )
    utils::write.csv(
      endpoint_tables$context,
      planned$core$endpoint_context_csv,
      row.names = FALSE
    )
  }
  if (write_local_search) {
    utils::write.csv(
      local_search_assignments,
      planned$core$local_search_assignments_csv,
      row.names = FALSE
    )
  }

  map_paths <- list()
  if (isTRUE(write_maps)) {
    map_dir <- file.path(output_dir, "maps")
    dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
    map_paths <- ftq_write_maps(result, template_image, map_dir, prefix, target_labels, overwrite)
  }

  meta <- list(
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    package = ftq_package_identity(),
    manifest = manifest,
    preflight = ftq_compact_preflight(preflight),
    inputs = if (!is.null(manifest)) ftq_file_inventory(ftq_manifest_input_paths(manifest)) else NULL,
    summary = ftq_result_summary(result, rsc, roi_summary),
    endpoint_analysis = if (write_endpoint) endpoint_tables$metadata else NULL,
    local_search_audit = if (write_local_search) {
      list(
        assignment_rows = nrow(local_search_assignments),
        search_radius_mm = attr(result, "fibrtraq_local_search_radius_mm"),
        spatial_units = attr(result, "fibrtraq_local_search_spatial_units"),
        spatial_units_scale_to_mm = attr(
          result,
          "fibrtraq_local_search_spatial_units_scale_to_mm"
        ),
        spatial_units_interpretation = attr(
          result,
          "fibrtraq_local_search_spatial_units_interpretation"
        )
      )
    } else {
      NULL
    },
    files = list(
      result_csv = result_csv,
      rsc_csv = rsc_csv,
      fate_csv = fate_csv,
      roi_summary_csv = roi_summary_csv,
      endpoint_summary_csv = if (write_endpoint) planned$core$endpoint_summary_csv else NULL,
      endpoint_penetration_csv = if (write_endpoint) planned$core$endpoint_penetration_csv else NULL,
      endpoint_context_csv = if (write_endpoint) planned$core$endpoint_context_csv else NULL,
      local_search_assignments_csv = if (write_local_search) {
        planned$core$local_search_assignments_csv
      } else {
        NULL
      },
      map_paths = map_paths,
      context = ftq_compact_context(context)
    )
  )
  jsonlite::write_json(meta, metadata_json, auto_unbox = TRUE, pretty = TRUE, null = "null")

  out <- list(
    result_csv = result_csv,
    rsc_csv = rsc_csv,
    fate_csv = fate_csv,
    roi_summary_csv = roi_summary_csv,
    metadata_json = metadata_json,
    validation_json = validation_json,
    map_paths = map_paths,
    context = ftq_compact_context(context)
  )
  if (write_endpoint) {
    out$endpoint_summary_csv <- planned$core$endpoint_summary_csv
    out$endpoint_penetration_csv <- planned$core$endpoint_penetration_csv
    out$endpoint_context_csv <- planned$core$endpoint_context_csv
  }
  if (write_local_search) {
    out$local_search_assignments_csv <-
      planned$core$local_search_assignments_csv
  }
  out
}

#' Validate a FibrTraQ run
#'
#' @param result_csv Per-voxel CSV path or data.frame.
#' @param rsc_csv RSC CSV path or data.frame.
#' @param map_paths Optional named list of NIfTI map paths.
#' @param expected Optional list of known-good result/rsc/map paths to compare.
#' @param validation_path Optional JSON output path.
#' @param tolerance Numeric tolerance for summary comparisons.
#' @param roi_summary_csv Optional ROI-summary CSV/data frame used to reconcile
#'   endpoint target-intersection denominators.
#' @param endpoint_summary Optional endpoint-summary CSV/data frame.
#' @param endpoint_penetration Optional endpoint-penetration CSV/data frame.
#' @param endpoint_context Optional endpoint-context CSV/data frame.
#' @param local_search_assignments Optional local-search end-voxel assignment
#'   CSV/data frame. When supplied, it is reconciled against the ROI summary.
#' @return A validation result list.
#' @export
ftq_validate_run <- function(result_csv,
                             rsc_csv,
                             map_paths = NULL,
                             expected = NULL,
                             validation_path = NULL,
                             tolerance = 1e-8,
                             roi_summary_csv = NULL,
                             endpoint_summary = NULL,
                             endpoint_penetration = NULL,
                             endpoint_context = NULL,
                             local_search_assignments = NULL) {
  result <- ftq_read_table(result_csv)
  rsc <- ftq_read_table(rsc_csv)
  checks <- list()
  add_check <- function(name, passed, detail = "") {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name,
      passed = isTRUE(passed),
      detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  }

  required <- c("fNode_Coords", "tracks", "local", "leave_no_target", "none")
  add_check("required_columns", all(required %in% names(result)),
            paste(setdiff(required, names(result)), collapse = ", "))
  reserved <- c("fNode_Coords", "tracks", "local", "leave_no_target", "none", "cNode_Coords")
  target_cols <- setdiff(names(result), reserved)
  add_check("has_target_columns", length(target_cols) > 0L, paste(target_cols, collapse = ", "))

  numeric_cols <- intersect(c("tracks", "local", "leave_no_target", "none", target_cols), names(result))
  nonnegative <- all(vapply(result[numeric_cols], function(x) all(!is.na(x) & x >= 0), logical(1)))
  add_check("nonnegative_counts", nonnegative)
  if (all(c("tracks", "local", "leave_no_target") %in% names(result))) {
    leave_with_target <- result$tracks - result$local - result$leave_no_target
    add_check("leave_with_target_nonnegative", all(leave_with_target >= 0),
              paste("min =", min(leave_with_target)))
  }
  if ("none" %in% names(result)) {
    add_check("none_matches_zero_tracks", identical(as.integer(result$none), as.integer(result$tracks == 0)))
  }

  expected_rsc <- rsc_summary(result, percent = TRUE, warn_legacy = FALSE)
  if (length(target_cols)) {
    expected_rsc <- expected_rsc[match(rsc$node, expected_rsc$node), , drop = FALSE]
    expected_assoc <- ftq_rsc_association_column(expected_rsc)
    observed_assoc <- ftq_rsc_association_column(rsc)
    rsc_ok <- all(abs(expected_assoc - observed_assoc) <= tolerance) &&
      all(abs(expected_rsc$voxels_connected - rsc$voxels_connected) <= tolerance) &&
      all(abs(expected_rsc$rsc_percent - rsc$rsc_percent) <= tolerance) &&
      all(abs(ftq_rsc_total_voxels(expected_rsc) - ftq_rsc_total_voxels(rsc)) <= tolerance)
    add_check("rsc_matches_per_voxel_table", rsc_ok)
  }

  if (!is.null(map_paths) && length(map_paths)) {
    map_check <- ftq_validate_maps(result, map_paths)
    checks <- c(checks, map_check)
  }
  if (!is.null(expected)) {
    expected_check <- ftq_validate_expected(
      result,
      rsc,
      map_paths,
      expected,
      tolerance,
      endpoint_summary = endpoint_summary,
      endpoint_penetration = endpoint_penetration,
      endpoint_context = endpoint_context
    )
    checks <- c(checks, expected_check)
  }

  biological_qc_status <- "unavailable"
  biological_qc_warnings <- character()
  endpoint_inputs <- list(endpoint_summary, endpoint_penetration, endpoint_context)
  endpoint_supplied <- !vapply(endpoint_inputs, is.null, logical(1))
  if (any(endpoint_supplied)) {
    if (!all(endpoint_supplied)) {
      add_check(
        "endpoint_tables_complete",
        FALSE,
        "endpoint_summary, endpoint_penetration, and endpoint_context must be supplied together"
      )
    } else {
      endpoint_validation <- ftq_validate_endpoint_tables(
        endpoint_summary,
        endpoint_penetration,
        endpoint_context,
        roi_summary = roi_summary_csv,
        tolerance = tolerance
      )
      if (nrow(endpoint_validation$checks)) {
        checks <- c(
          checks,
          split(endpoint_validation$checks, seq_len(nrow(endpoint_validation$checks)))
        )
      }
      biological_qc_status <- endpoint_validation$biological_qc_status
      biological_qc_warnings <- endpoint_validation$biological_qc_warnings
    }
  }

  if (!is.null(local_search_assignments)) {
    local_search_checks <- ftq_validate_local_search_assignments(
      local_search_assignments,
      roi_summary = roi_summary_csv,
      target_labels = target_cols,
      tolerance = tolerance
    )
    if (nrow(local_search_checks)) {
      checks <- c(
        checks,
        split(local_search_checks, seq_len(nrow(local_search_checks)))
      )
    }
  }

  checks_df <- do.call(rbind, checks)
  rownames(checks_df) <- NULL
  out <- list(
    passed = all(checks_df$passed),
    checks = checks_df,
    biological_qc_status = biological_qc_status,
    biological_qc_warnings = biological_qc_warnings,
    checked_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  if (!is.null(validation_path)) {
    jsonlite::write_json(out, validation_path, auto_unbox = TRUE, pretty = TRUE)
  }
  out
}

#' Build a FibrTraQ HTML report and figure bundle
#'
#' @param result_csv Per-voxel CSV path or data.frame.
#' @param rsc_csv RSC CSV path or data.frame.
#' @param manifest Optional manifest list.
#' @param output_dir Report output directory.
#' @param map_paths Optional named map path list.
#' @param validation Optional validation result from [ftq_validate_run()].
#' @param preflight Optional preflight result from [ftq_preflight()].
#' @param tractogram_context Optional tractogram-derived context output or
#'   density image path used when no anatomical background is available.
#' @param endpoint_summary Optional endpoint-summary CSV/data frame.
#' @param endpoint_penetration Optional endpoint-penetration CSV/data frame.
#' @param endpoint_context Optional endpoint-context CSV/data frame.
#' @param local_search_assignments Optional local-search end-voxel assignment
#'   CSV/data frame written by [ftq_write_outputs()]. When supplied, the report
#'   summarizes end-to-node distances and partitions primary local-search
#'   connections by whether the full voxelized path also enters the cNode.
#' @param roi_summary Optional ROI-summary CSV/data frame. When supplied, the
#'   report uses ROI-level unique connection counts and the unique-streamline
#'   fate denominator; RSC retains its focal-voxel denominator.
#'   Legacy calls without it remain supported and are
#'   explicitly labeled as association-weighted where necessary.
#' @param overwrite Whether existing report files may be overwritten.
#' @param report_slices Optional distinct, positive one-based axial slice indices.
#'   NULL uses manifest report.slices when supplied, otherwise automatic selection.
#' @return A list of report paths.
#' @export
ftq_report <- function(result_csv,
                       rsc_csv,
                       manifest = NULL,
                       output_dir,
                       map_paths = NULL,
                       validation = NULL,
                       preflight = NULL,
                       overwrite = TRUE,
                       tractogram_context = NULL,
                       endpoint_summary = NULL,
                       endpoint_penetration = NULL,
                       endpoint_context = NULL,
                       local_search_assignments = NULL,
                       roi_summary = NULL,
                       report_slices = NULL) {
  manifest <- ftq_report_resolve_manifest(result_csv, manifest)
  result <- ftq_read_table(result_csv)
  rsc <- ftq_read_table(rsc_csv)
  if (!is.null(manifest)) {
    if (!is.list(manifest)) manifest <- ftq_read_manifest(manifest)
    analysis <- manifest$analysis
    if (is.null(analysis)) analysis <- list()
    report <- manifest$report
    if (is.null(report)) report <- list()
    manifest$focal_assignment_mode <- match.arg(
      as.character(ftq_value(
        analysis$focal_assignment_mode,
        manifest$focal_assignment_mode,
        "all_focal_voxels"
      ))[1L],
      c("all_focal_voxels", "first_intersection")
    )
    manifest$streamline_voxelization_mode <- match.arg(
      as.character(ftq_value(
        analysis$streamline_voxelization_mode,
        manifest$streamline_voxelization_mode,
        "sample_points"
      ))[1L],
      c("sample_points", "segment_interpolation")
    )
    manifest$target_intersection_mode <-
      ftq_normalize_target_intersection_mode(ftq_value(
        analysis$target_intersection_mode,
        manifest$target_intersection_mode,
        "full_streamline"
      ))
    manifest$local_search_radius_mm <- ftq_validate_local_search_radius(
      ftq_value(
        analysis$local_search_radius_mm,
        manifest$local_search_radius_mm,
        2
      )
    )
    manifest$endpoint_analysis <- match.arg(
      as.character(ftq_value(
        analysis$endpoint_analysis,
        manifest$endpoint_analysis,
        "exact_and_dilated1"
      ))[1L],
      c("exact_and_dilated1", "none")
    )
    manifest$background_image <- ftq_value(
      manifest$background_image,
      report$background_image,
      NULL
    )
    if (!is.null(manifest$background_image) &&
        length(manifest$background_image) == 1L &&
        !is.na(manifest$background_image) &&
        nzchar(manifest$background_image) &&
        !is.null(manifest$dataset_root) &&
        !grepl("^(/|[A-Za-z]:[/\\\\])", manifest$background_image)) {
      manifest$background_image <- ftq_manifest_path(
        manifest$background_image,
        normalizePath(
          path.expand(manifest$dataset_root),
          mustWork = FALSE
        )
      )
    }
  }
  prefix <- if (!is.null(manifest$prefix)) manifest$prefix else "fibrtraq"
  report_slices <- ftq_report_slices(ftq_value(report_slices, manifest$report$slices))
  template_candidates <- as.character(c(map_paths$tracks, manifest$fnode))
  template_candidates <- template_candidates[!is.na(template_candidates) & file.exists(template_candidates)]
  template_image <- if (length(template_candidates)) template_candidates[[1L]] else NULL
  if (!is.null(template_image) && file.exists(template_image)) {
    report_slices <- ftq_report_slices(report_slices, dim(RNifti::readNifti(template_image))[3L])
  }
  manifest$report$slices <- report_slices
  report_paths <- ftq_report_output_paths(
    output_dir, prefix, write_endpoint = TRUE, write_local_search = TRUE
  )
  report_html <- report_paths$html
  if (!isTRUE(overwrite) && dir.exists(output_dir)) {
    ftq_check_overwrite(list.files(output_dir, recursive = TRUE, full.names = TRUE), FALSE)
  }

  fig_dir <- file.path(output_dir, "report-figures")
  table_dir <- file.path(output_dir, "tables")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  roi_summary_data <- if (is.null(roi_summary)) {
    attr(result, "fibrtraq_roi_summary")
  } else {
    ftq_read_table(roi_summary)
  }
  connection_summary <- ftq_report_connection_summary(
    result,
    rsc,
    roi_summary = roi_summary_data
  )
  fate_summary <- ftq_report_fate_summary(
    result,
    roi_summary = roi_summary_data
  )
  top_tables <- ftq_write_top_tables(result, table_dir, prefix)
  voxel_exports <- ftq_report_write_voxel_exports(
    result, output_dir, prefix, template_image = template_image, overwrite = overwrite
  )
  top_tables <- c(top_tables, voxel_exports$tables)
  top_tables$connection_summary <- report_paths$tables$connection_summary
  top_tables$fate_summary <- report_paths$tables$fate_summary
  utils::write.csv(connection_summary, top_tables$connection_summary, row.names = FALSE)
  utils::write.csv(
    fate_summary,
    top_tables$fate_summary,
    row.names = FALSE
  )
  endpoint_inputs <- list(
    summary = endpoint_summary,
    penetration = endpoint_penetration,
    context = endpoint_context
  )
  endpoint_supplied <- !vapply(endpoint_inputs, is.null, logical(1))
  endpoint_data <- NULL
  if (all(endpoint_supplied)) {
    endpoint_data <- list(
      summary = ftq_read_table(endpoint_summary),
      penetration = ftq_read_table(endpoint_penetration),
      context = ftq_read_table(endpoint_context),
      source_paths = endpoint_inputs,
      missing = FALSE
    )
    endpoint_data$penetration_summary <- ftq_endpoint_penetration_summary(
      endpoint_data$penetration,
      target_labels = unique(as.character(endpoint_data$summary$target))
    )
    penetration_summary_path <- report_paths$tables$endpoint_penetration_summary
    utils::write.csv(
      endpoint_data$penetration_summary,
      penetration_summary_path,
      row.names = FALSE
    )
    top_tables$endpoint_penetration_summary <- penetration_summary_path
  } else if (any(endpoint_supplied) ||
             (isTRUE(manifest$endpoint_analysis_explicit) &&
              identical(ftq_value(manifest$endpoint_analysis, "none"), "exact_and_dilated1"))) {
    endpoint_data <- list(missing = TRUE)
  }
  local_search_data <- NULL
  if (!is.null(local_search_assignments)) {
    local_search_data <- list(
      assignments = ftq_read_table(local_search_assignments),
      source_path = local_search_assignments
    )
    local_search_data$summary <- ftq_local_search_assignment_summary(
      local_search_data$assignments
    )
    local_search_data$endpoint_status <- ftq_local_search_endpoint_status(
      local_search_data$assignments
    )
    top_tables$local_search_summary <- report_paths$tables$local_search_summary
    top_tables$local_search_endpoint_status <- report_paths$tables$local_search_endpoint_status
    utils::write.csv(
      local_search_data$summary,
      top_tables$local_search_summary,
      row.names = FALSE
    )
    utils::write.csv(
      local_search_data$endpoint_status,
      top_tables$local_search_endpoint_status,
      row.names = FALSE
    )
  }
  figs <- report_paths$figures[c("fate", "rsc", "tracks", "spatial_tracks",
                                "spatial_targets", "spatial_leave_no_target")]
  if (!is.null(endpoint_data) && !isTRUE(endpoint_data$missing)) {
    figs$endpoint_classification <- report_paths$figures$endpoint_classification
    figs$endpoint_comparison <- report_paths$figures$endpoint_comparison
    figs$endpoint_penetration <- report_paths$figures$endpoint_penetration
    ftq_plot_endpoint_classification(endpoint_data$summary, figs$endpoint_classification)
    ftq_plot_endpoint_comparison(endpoint_data$summary, figs$endpoint_comparison)
    ftq_plot_endpoint_penetration(
      endpoint_data$penetration_summary,
      figs$endpoint_penetration
    )
  }
  if (!is.null(local_search_data)) {
    figs$local_search_distances <- report_paths$figures$local_search_distances
    ftq_plot_local_search_distances(
      local_search_data$assignments,
      figs$local_search_distances
    )
  }
  ftq_plot_fate_summary(result, figs$fate, fate = fate_summary)
  ftq_plot_rsc_summary(rsc, figs$rsc, connection_summary = connection_summary)
  ftq_plot_tracks_distribution(result, figs$tracks)
  ftq_plot_spatial_bundle(
    result,
    map_paths,
    figs,
    background_image = manifest$background_image,
    tractogram_context = tractogram_context,
    context_source = ftq_value(manifest$context_source, "auto"),
    target_intersection_mode = ftq_value(
      manifest$target_intersection_mode,
      "full_streamline"
    ),
    local_search_radius_mm = ftq_value(manifest$local_search_radius_mm, 2),
    display_quantile = 0.99,
    report_slices = report_slices
  )
  if (length(voxel_exports$maps)) {
    figs$local_fraction <- report_paths$figures$local_fraction
    arr <- RNifti::readNifti(voxel_exports$maps$local_fraction)
    bg <- ftq_report_background(manifest$background_image, tractogram_context,
                               dim(arr)[1:3], ftq_value(manifest$context_source, "auto"))
    ftq_plot_montage(arr, figs$local_fraction, "Local fraction", background = bg$image,
                     legend_label = "Local / assigned streamlines", context = TRUE,
                     background_label = bg$label, slices = report_slices,
                     validity_mask = as.array(RNifti::readNifti(voxel_exports$maps$local_fraction_valid)) > 0)
  }
  figs <- figs[file.exists(unlist(figs, use.names = FALSE))]

  guide <- ftq_write_output_guide(output_dir)
  source_paths <- list(per_voxel = result_csv, rsc = rsc_csv, roi_summary = roi_summary)
  metadata_path <- NULL
  if (is.character(result_csv) && length(result_csv) == 1L &&
      grepl("_fibrtraq_per-voxel-hits\\.csv$", result_csv)) {
    candidate <- sub("_fibrtraq_per-voxel-hits\\.csv$", "_fibrtraq_run-metadata.json", result_csv)
    if (file.exists(candidate)) metadata_path <- candidate
    candidate <- sub("_fibrtraq_per-voxel-hits\\.csv$", "_fibrtraq_validation.json", result_csv)
    if (file.exists(candidate)) source_paths$validation <- candidate
  }
  ftq_write_html_report(
    report_html,
    result,
    rsc,
    manifest,
    validation,
    preflight,
    figs,
    top_tables,
    tractogram_context,
    endpoint_data = endpoint_data,
    local_search_data = local_search_data,
    connection_summary = connection_summary,
    fate_summary = fate_summary,
    source_paths = source_paths,
    map_paths = c(map_paths, voxel_exports$maps),
    metadata_path = metadata_path,
    guide = guide
  )
  list(html = report_html, figures = figs, tables = top_tables,
       report_maps = voxel_exports$maps, guide = guide)
}

ftq_value <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (!is.null(v)) return(v)
  }
  NULL
}

ftq_manifest_path <- function(path, root) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path)) return(NULL)
  path <- path.expand(path)
  if (ftq_is_absolute_path(path)) return(normalizePath(path, mustWork = FALSE))
  normalizePath(file.path(root, path), mustWork = FALSE)
}

ftq_is_absolute_path <- function(path) {
  grepl("^/", path) ||
    grepl("^[A-Za-z]:[/\\\\]", path) ||
    grepl("^[/\\\\]{2}[^/\\\\]", path)
}

ftq_parse_targets <- function(targets, root) {
  if (is.list(targets) && length(targets) && !is.null(names(targets)) && all(nzchar(names(targets)))) {
    labels <- names(targets)
    paths <- vapply(targets, function(x) {
      if (is.list(x)) x <- x$path
      ftq_manifest_path(as.character(x), root)
    }, character(1))
  } else if (is.list(targets)) {
    labels <- vapply(targets, function(x) as.character(x$label), character(1))
    paths <- vapply(targets, function(x) ftq_manifest_path(as.character(x$path), root), character(1))
  } else {
    stop("targets must be a named mapping or a list of {label, path} objects.")
  }
  if (any(!nzchar(labels)) || anyDuplicated(labels)) stop("Target labels must be non-empty and unique.")
  list(labels = labels, paths = unname(paths))
}

ftq_parse_endpoint_contexts <- function(contexts, root) {
  if (is.null(contexts) || !length(contexts)) {
    return(list(labels = character(), paths = character()))
  }
  parsed <- ftq_parse_targets(contexts, root)
  if ("reference_grid" %in% parsed$labels) {
    stop("endpoint_contexts cannot use the reserved label 'reference_grid'.")
  }
  parsed
}

ftq_parse_expected <- function(expected, root) {
  if (is.null(expected)) return(NULL)
  if (!is.list(expected)) stop("expected must be a mapping of known-good output paths.")
  result_csv <- ftq_value(
    expected$result_csv,
    expected$per_voxel_csv,
    expected$per_voxel_hits,
    expected$result
  )
  rsc_csv <- ftq_value(expected$rsc_csv, expected$rsc_summary, expected$rsc)
  endpoint_summary_csv <- ftq_value(
    expected$endpoint_summary_csv,
    expected$endpoint_summary
  )
  endpoint_penetration_csv <- ftq_value(
    expected$endpoint_penetration_csv,
    expected$endpoint_penetration
  )
  endpoint_context_csv <- ftq_value(
    expected$endpoint_context_csv,
    expected$endpoint_context
  )
  map_paths <- ftq_value(expected$map_paths, expected$maps, NULL)

  maps <- NULL
  if (!is.null(map_paths)) {
    if (!is.list(map_paths) || is.null(names(map_paths)) || any(!nzchar(names(map_paths)))) {
      stop("expected maps must be a named mapping of map keys to paths.")
    }
    maps <- lapply(map_paths, function(x) {
      if (is.list(x)) x <- x$path
      ftq_manifest_path(as.character(x), root)
    })
  }

  out <- list(
    result_csv = ftq_manifest_path(as.character(result_csv), root),
    rsc_csv = ftq_manifest_path(as.character(rsc_csv), root),
    endpoint_summary_csv = ftq_manifest_path(as.character(endpoint_summary_csv), root),
    endpoint_penetration_csv = ftq_manifest_path(as.character(endpoint_penetration_csv), root),
    endpoint_context_csv = ftq_manifest_path(as.character(endpoint_context_csv), root),
    map_paths = maps
  )
  if (is.null(out$result_csv) && is.null(out$rsc_csv) &&
      is.null(out$endpoint_summary_csv) && is.null(out$endpoint_penetration_csv) &&
      is.null(out$endpoint_context_csv) && !length(out$map_paths)) {
    stop("expected must include a result, RSC, endpoint, or map output.")
  }
  out
}

ftq_parse_tractogram_context <- function(context, root) {
  if (is.null(context)) {
    return(list(
      enabled = TRUE,
      density_image = NULL,
      output_dir = NULL,
      prefix = NULL,
      track_chunk_size = NULL,
      density_discretization = "round"
    ))
  }
  if (is.character(context) && length(context) == 1L && nzchar(context)) {
    context <- list(density_image = context)
  }
  if (!is.list(context)) stop("report.tractogram_context must be a mapping or image path.")
  density_image <- ftq_value(context$density_image, context$image, context$path, NULL)
  output_dir <- ftq_value(context$output_dir, NULL)
  list(
    enabled = !identical(context$enabled, FALSE),
    density_image = if (!is.null(density_image)) ftq_manifest_path(as.character(density_image), root) else NULL,
    output_dir = if (!is.null(output_dir)) ftq_manifest_path(as.character(output_dir), root) else NULL,
    prefix = ftq_value(context$prefix, NULL),
    track_chunk_size = ftq_value(context$track_chunk_size, NULL),
    density_discretization = as.character(ftq_value(context$density_discretization, "round"))[1]
  )
}

ftq_run_tractogram_context <- function(m, result, overwrite) {
  cfg <- m$tractogram_context
  if (is.null(cfg)) cfg <- ftq_parse_tractogram_context(NULL, m$dataset_root)
  existing <- ftq_context_density_image(cfg)
  if (!is.null(existing)) {
    if (!file.exists(existing)) {
      if (identical(m$context_source, "tractogram")) {
        stop("Configured tractogram context image does not exist: ", existing)
      }
      warning("Configured tractogram context image does not exist; plotting without fallback context.")
      return(NULL)
    }
    return(list(
      density_image = existing,
      summary_csv = NULL,
      metadata_json = NULL,
      source = "configured tractogram-derived density context"
    ))
  }
  if (!isTRUE(cfg$enabled) || identical(m$context_source, "none")) return(NULL)
  if (!isTRUE(m$write_report) && !identical(m$context_source, "tractogram")) return(NULL)
  has_background <- !is.null(m$background_image) && file.exists(m$background_image)
  if (has_background && m$context_source %in% c("auto", "anatomy")) return(NULL)
  should_generate <- identical(m$context_source, "tractogram") ||
    (identical(m$context_source, "auto") && !has_background && grepl("\\.(trk|tck)$", tolower(m$tractogram)))
  if (!should_generate) return(NULL)
  if (!grepl("\\.(trk|tck)$", tolower(m$tractogram))) {
    if (identical(m$context_source, "tractogram")) {
      stop("Tractogram-derived context requires a .trk or .tck tractogram.")
    }
    return(NULL)
  }
  selected_track_space <- ftq_value(attr(result, "fibrtraq_track_space"), m$track_space, "auto")
  selected_axis_flip <- ftq_value(attr(result, "fibrtraq_axis_flip"), m$axis_flip, c(FALSE, FALSE, FALSE))
  chunk_size <- as.integer(ftq_value(cfg$track_chunk_size, m$track_chunk_size, 50000L))[1]
  density_discretization <- as.character(ftq_value(cfg$density_discretization, "round"))[1]
  tryCatch(
    ftq_build_tractogram_context(
      tractogram = m$tractogram,
      template_image = m$fnode,
      output_dir = ftq_value(cfg$output_dir, m$output_dir),
      prefix = ftq_value(cfg$prefix, m$prefix),
      track_space = selected_track_space,
      axis_flip = selected_axis_flip,
      track_chunk_size = chunk_size,
      density_discretization = density_discretization,
      streamline_voxelization_mode = m$streamline_voxelization_mode,
      ref_image = ftq_value(m$ref_image, m$fnode),
      overwrite = overwrite,
      density_subdir = "maps",
      write_sidecars = FALSE
    ),
    error = function(e) {
      if (identical(m$context_source, "tractogram")) stop(e)
      warning("Could not build tractogram-derived density context: ", conditionMessage(e))
      NULL
    }
  )
}

ftq_compact_context <- function(context) {
  if (is.null(context)) return(NULL)
  if (is.character(context)) context <- list(density_image = context)
  metadata <- context$metadata
  summary <- ftq_value(context$summary, metadata$summary, NULL)
  # Store the occupancy summary once while retaining all density settings.
  if (is.list(metadata)) metadata$summary <- NULL
  list(
    density_image = ftq_context_density_image(context),
    summary_csv = ftq_value(context$summary_csv, NULL),
    metadata_json = ftq_value(context$metadata_json, NULL),
    source = ftq_value(context$source, "tractogram-derived density context"),
    summary = summary,
    metadata = metadata
  )
}

ftq_expected_paths <- function(expected) {
  if (is.null(expected)) return(character())
  unname(c(
    expected$result_csv,
    expected$rsc_csv,
    expected$endpoint_summary_csv,
    expected$endpoint_penetration_csv,
    expected$endpoint_context_csv,
    unlist(expected$map_paths, use.names = FALSE)
  ))
}

ftq_manifest_input_paths <- function(m) {
  endpoint_context_inputs <- if (length(m$endpoint_context_paths)) {
    stats::setNames(
      as.list(m$endpoint_context_paths),
      paste0("endpoint_context:", m$endpoint_context_labels)
    )
  } else {
    list()
  }
  paths <- c(
    list(
      fnode = m$fnode,
      local = m[["local"]],
      tractogram = m$tractogram,
      ref_image = m$ref_image,
      background_image = m$background_image
    ),
    stats::setNames(as.list(m$target_paths), paste0("target:", m$target_labels)),
    endpoint_context_inputs
  )
  paths[!vapply(paths, is.null, logical(1))]
}

ftq_planned_output_paths <- function(output_dir,
                                     prefix,
                                     target_labels = NULL,
                                     write_maps = TRUE,
                                     write_endpoint = FALSE,
                                     write_local_search = FALSE) {
  core <- list(
    result_csv = file.path(output_dir, paste0(prefix, "_fibrtraq_per-voxel-hits.csv")),
    rsc_csv = file.path(output_dir, paste0(prefix, "_fibrtraq_rsc-summary.csv")),
    fate_csv = file.path(output_dir, paste0(prefix, "_fibrtraq_fate-summary.csv")),
    roi_summary_csv = file.path(output_dir, paste0(prefix, "_fibrtraq_roi-summary.csv")),
    metadata_json = file.path(output_dir, paste0(prefix, "_fibrtraq_run-metadata.json")),
    validation_json = file.path(output_dir, paste0(prefix, "_fibrtraq_validation.json"))
  )
  if (isTRUE(write_endpoint)) {
    core$endpoint_summary_csv <- file.path(
      output_dir,
      paste0(prefix, "_fibrtraq_endpoint-summary.csv")
    )
    core$endpoint_penetration_csv <- file.path(
      output_dir,
      paste0(prefix, "_fibrtraq_endpoint-penetration.csv")
    )
    core$endpoint_context_csv <- file.path(
      output_dir,
      paste0(prefix, "_fibrtraq_endpoint-context.csv")
    )
  }
  if (isTRUE(write_local_search)) {
    core$local_search_assignments_csv <- file.path(
      output_dir,
      paste0(prefix, "_fibrtraq_local-search-assignments.csv")
    )
  }
  maps <- list()
  if (isTRUE(write_maps)) {
    labels <- c("tracks", "local", "leave-no-target", "leave-with-target", "tracks-binary")
    if (!is.null(target_labels)) labels <- c(labels, paste0("target-", target_labels))
    maps <- stats::setNames(
      lapply(labels, function(nm) {
        file.path(output_dir, "maps", paste0(prefix, "_fibrtraq_per-voxel-hits_",
                                             ftq_sanitize_label(nm), ".nii.gz"))
      }),
      labels
    )
  }
  list(core = core, maps = maps)
}

ftq_read_volume_for_preflight <- function(path) {
  if (is.null(path) || !file.exists(path)) return(list(data = NULL, error = "File does not exist."))
  tryCatch(
    list(data = ftq_read_volume_data(path), error = NULL),
    error = function(e) list(data = NULL, error = conditionMessage(e))
  )
}

ftq_volume_info <- function(x) {
  full_dims <- dim(x)
  full_dims <- if (is.null(full_dims)) integer() else as.integer(full_dims)
  is_3d <- length(full_dims) == 3L
  dims <- if (length(full_dims) >= 3L) {
    as.integer(full_dims[1:3])
  } else {
    rep(NA_integer_, 3L)
  }
  vals <- as.numeric(x)
  finite <- is.finite(vals)
  finite_vals <- vals[finite]
  list(
    dims = dims,
    full_dims = full_dims,
    dimensionality = length(full_dims),
    is_3d = is_3d,
    finite_voxels = sum(finite),
    nonzero_voxels = sum(finite & vals != 0),
    min = if (length(finite_vals)) min(finite_vals) else NA_real_,
    max = if (length(finite_vals)) max(finite_vals) else NA_real_,
    binary_values = if (length(finite_vals)) all(finite_vals %in% c(0, 1)) else FALSE
  )
}

ftq_dim_label <- function(dims) {
  if (is.null(dims) || !length(dims) || anyNA(dims)) return("unknown")
  paste(as.integer(dims), collapse = "x")
}

ftq_preflight_grid_check <- function(reference_path,
                                     other_path,
                                     reference_label,
                                     other_label,
                                     spatial_units_to_mm = FALSE) {
  fallback <- paste(other_label, "=", other_path, reference_label, "=", reference_path)
  out <- tryCatch({
    reference <- RNifti::readNifti(reference_path, internal = FALSE)
    other <- RNifti::readNifti(other_path, internal = FALSE)
    reference_full_dims <- as.integer(dim(reference))
    other_full_dims <- as.integer(dim(other))
    reference_is_3d <- length(reference_full_dims) == 3L
    other_is_3d <- length(other_full_dims) == 3L
    reference_dims <- if (length(reference_full_dims) >= 3L) {
      reference_full_dims[1:3]
    } else {
      rep(NA_integer_, 3L)
    }
    other_dims <- if (length(other_full_dims) >= 3L) {
      other_full_dims[1:3]
    } else {
      rep(NA_integer_, 3L)
    }
    dims_match <- reference_is_3d && other_is_3d &&
      identical(other_dims, reference_dims)
    affine_available <- dims_match &&
      ftq_has_spatial_grid_metadata(reference) &&
      ftq_has_spatial_grid_metadata(other)
    affine_match <- if (!affine_available) {
      NA
    } else if (isTRUE(spatial_units_to_mm)) {
      ftq_local_search_affine_equal(other, reference)
    } else {
      ftq_affine_equal(other, reference)
    }
    reason <- if (!reference_is_3d || !other_is_3d) {
      paste0(
        "images must be exactly 3-D: ",
        other_label, " = ", ftq_dim_label(other_full_dims),
        " (", length(other_full_dims), " dimensions), ",
        reference_label, " = ", ftq_dim_label(reference_full_dims),
        " (", length(reference_full_dims), " dimensions)"
      )
    } else if (!dims_match) {
      paste0("dimensions differ: ", other_label, " = ", ftq_dim_label(other_dims),
             ", ", reference_label, " = ", ftq_dim_label(reference_dims))
    } else if (!affine_available) {
      paste0(
        "dimensions match; affine/grid comparison is unavailable because ",
        "both images do not declare a valid qform/sform transform"
      )
    } else if (!isTRUE(affine_match)) {
      paste0("dimensions match but affine/grid differs: ", other_label,
             " and ", reference_label, " are not on the same spatial grid")
    } else {
      paste0(other_label, " and ", reference_label, " dimensions and affine/grid match.")
    }
    if ((reference_is_3d && other_is_3d && !dims_match) ||
        (affine_available && !isTRUE(affine_match))) {
      reason <- paste0(
        reason,
        "; use ftq_make_union_mask() when multiple aligned regional masks are available, ",
        "or ftq_regrid_mask() only when the source affine/header is trusted."
      )
    }
    list(
      dims_match = dims_match,
      reference_is_3d = reference_is_3d,
      other_is_3d = other_is_3d,
      affine_available = affine_available,
      affine_match = affine_match,
      grid_match = dims_match && affine_available && isTRUE(affine_match),
      detail = reason
    )
  }, error = function(e) {
    list(
      dims_match = FALSE,
      reference_is_3d = FALSE,
      other_is_3d = FALSE,
      affine_available = FALSE,
      affine_match = FALSE,
      grid_match = FALSE,
      detail = paste(fallback, conditionMessage(e))
    )
  })
  out
}

ftq_preflight_failures <- function(preflight) {
  failed <- preflight$checks[!preflight$checks$passed & preflight$checks$severity == "error", , drop = FALSE]
  if (!nrow(failed)) return(character())
  paste0("- ", failed$check, ": ", failed$detail)
}

ftq_compact_preflight <- function(preflight) {
  if (is.null(preflight)) return(NULL)
  list(
    passed = preflight$passed,
    checks = preflight$checks,
    checked_at = preflight$checked_at
  )
}

ftq_file_inventory <- function(paths) {
  paths <- paths[!vapply(paths, function(x) is.null(x) || is.na(x) || !nzchar(x), logical(1))]
  if (!length(paths)) return(data.frame())
  out <- lapply(names(paths), function(nm) {
    path <- paths[[nm]]
    exists <- file.exists(path)
    info <- if (exists) file.info(path) else NULL
    data.frame(
      name = nm,
      path = path,
      exists = exists,
      size = if (exists) unname(info$size) else NA_real_,
      mtime = if (exists) format(info$mtime, "%Y-%m-%dT%H:%M:%S%z") else NA_character_,
      checksum_md5 = if (exists) unname(tools::md5sum(path)) else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

ftq_result_summary <- function(result, rsc, roi_summary = ftq_roi_summary(result)) {
  leave_with_target <- result$tracks - result$local - result$leave_no_target
  reserved <- c("fNode_Coords", "tracks", "local", "leave_no_target", "none", "cNode_Coords")
  target_cols <- setdiff(names(result), reserved)
  list(
    focal_voxels = nrow(result),
    focal_assignment_mode = ftq_value(attr(result, "fibrtraq_focal_assignment_mode"), NA_character_),
    streamline_voxelization_mode = ftq_value(attr(result, "fibrtraq_streamline_voxelization_mode"), NA_character_),
    target_intersection_mode = ftq_value(attr(result, "fibrtraq_target_intersection_mode"), NA_character_),
    local_search_radius_mm = ftq_value(attr(result, "fibrtraq_local_search_radius_mm"), NA_real_),
    local_search_spatial_units = ftq_value(
      attr(result, "fibrtraq_local_search_spatial_units"),
      NA_character_
    ),
    local_search_spatial_units_scale_to_mm = ftq_value(
      attr(result, "fibrtraq_local_search_spatial_units_scale_to_mm"),
      NA_real_
    ),
    local_search_spatial_units_interpretation = ftq_value(
      attr(result, "fibrtraq_local_search_spatial_units_interpretation"),
      NA_character_
    ),
    reverse_streamlines = ftq_value(attr(result, "fibrtraq_reverse_streamlines"), NA),
    voxelwise_track_associations = sum(result$tracks),
    voxelwise_local_associations = sum(result$local),
    voxelwise_leave_no_target_associations = sum(result$leave_no_target),
    voxelwise_leave_with_target_associations = sum(leave_with_target),
    voxelwise_target_associations = as.list(stats::setNames(vapply(target_cols, function(nm) sum(result[[nm]]), numeric(1)), target_cols)),
    roi_summary = roi_summary,
    rsc = rsc
  )
}

ftq_rsc_association_column <- function(rsc) {
  if ("target_associations" %in% names(rsc)) return(as.numeric(rsc$target_associations))
  if ("focal_voxel_streamline_target_associations" %in% names(rsc)) {
    return(as.numeric(rsc$focal_voxel_streamline_target_associations))
  }
  if ("tracks" %in% names(rsc)) return(as.numeric(rsc$tracks))
  stop("RSC table is missing target association columns.")
}

ftq_rsc_total_voxels <- function(rsc) {
  if ("total_focal_voxels" %in% names(rsc)) return(as.numeric(rsc$total_focal_voxels))
  if ("total_origin_voxels" %in% names(rsc)) return(as.numeric(rsc$total_origin_voxels))
  stop("RSC table is missing total focal voxel columns.")
}

ftq_roi_summary <- function(result) {
  summary <- attr(result, "fibrtraq_roi_summary")
  if (is.data.frame(summary)) return(summary)
  data.frame(
    summary_available = FALSE,
    focal_assignment_mode = ftq_value(attr(result, "fibrtraq_focal_assignment_mode"), NA_character_),
    streamline_voxelization_mode = ftq_value(attr(result, "fibrtraq_streamline_voxelization_mode"), NA_character_),
    target_intersection_mode = ftq_value(attr(result, "fibrtraq_target_intersection_mode"), NA_character_),
    local_search_radius_mm = ftq_value(attr(result, "fibrtraq_local_search_radius_mm"), NA_real_),
    local_search_spatial_units = ftq_value(
      attr(result, "fibrtraq_local_search_spatial_units"),
      NA_character_
    ),
    local_search_spatial_units_scale_to_mm = ftq_value(
      attr(result, "fibrtraq_local_search_spatial_units_scale_to_mm"),
      NA_real_
    ),
    local_search_spatial_units_interpretation = ftq_value(
      attr(result, "fibrtraq_local_search_spatial_units_interpretation"),
      NA_character_
    ),
    local_roi_defaulted_to_focal = ftq_value(attr(result, "fibrtraq_local_roi_defaulted_to_focal"), NA),
    note = "ROI-level unique-streamline summary is only available on freshly computed FibrTraQ results.",
    stringsAsFactors = FALSE
  )
}

ftq_prefix <- function(subject, session = NULL, desc) {
  parts <- c(subject, session)
  parts <- parts[!is.null(parts) & !is.na(parts) & nzchar(parts)]
  paste0(paste(parts, collapse = "_"), "_desc-", desc)
}

ftq_preflight_manifest <- function(m) {
  paths <- c(
    m$fnode, m[["local"]], m$tractogram, m$target_paths,
    m$endpoint_context_paths, m$ref_image, ftq_expected_paths(m$expected)
  )
  missing <- paths[!is.na(paths) & !file.exists(paths)]
  if (length(missing)) stop("Input file(s) do not exist:\n", paste(missing, collapse = "\n"))
  invisible(TRUE)
}

ftq_read_volume_data <- function(path) {
  low <- tolower(path)
  if (grepl("\\.m(gz|gh)$", low)) {
    freesurferformats::read.fs.mgh(path, drop_empty_dims = FALSE)
  } else {
    freesurferformats::read.nifti1.data(path)
  }
}

ftq_read_table <- function(x) {
  if (is.data.frame(x)) return(x)
  df <- utils::read.csv(x, check.names = FALSE, colClasses = "character")
  for (nm in setdiff(names(df), c("fNode_Coords", "node"))) {
    suppressWarnings(num <- as.numeric(df[[nm]]))
    if (!anyNA(num) || all(is.na(df[[nm]]) == is.na(num))) df[[nm]] <- num
  }
  ftq_report_read_voxel_view(x, df)
}

ftq_check_overwrite <- function(paths, overwrite) {
  existing <- paths[file.exists(paths)]
  if (length(existing) && !isTRUE(overwrite)) {
    stop("Output file(s) already exist. Use overwrite=TRUE to replace:\n",
         paste(existing, collapse = "\n"))
  }
}

ftq_coord_matrix <- function(keys) {
  m <- do.call(rbind, strsplit(as.character(keys), "_", fixed = TRUE))
  matrix(as.integer(m), ncol = 3L, dimnames = list(NULL, c("x", "y", "z")))
}

ftq_fate_summary <- function(result) {
  leave_with_target <- result$tracks - result$local - result$leave_no_target
  total <- sum(result$tracks)
  data.frame(
    fate = c("local", "leave_no_target", "leave_with_target"),
    tracks = c(sum(result$local), sum(result$leave_no_target), sum(leave_with_target)),
    percent = if (total > 0) 100 * c(sum(result$local), sum(result$leave_no_target), sum(leave_with_target)) / total else 0,
    stringsAsFactors = FALSE
  )
}

ftq_report_roi_row <- function(roi_summary) {
  if (is.null(roi_summary) || !is.data.frame(roi_summary) || nrow(roi_summary) < 1L) {
    return(NULL)
  }
  if ("summary_available" %in% names(roi_summary)) {
    available <- tolower(as.character(roi_summary$summary_available[[1L]]))
    if (available %in% c("false", "0", "no", "n")) return(NULL)
  }
  roi_summary[1L, , drop = FALSE]
}

ftq_report_connection_summary <- function(result, rsc, roi_summary = NULL) {
  required <- c("node", "voxels_connected", "rsc_percent")
  missing <- setdiff(required, names(rsc))
  if (length(missing)) {
    stop(
      "RSC table is missing columns required for the report: ",
      paste(missing, collapse = ", ")
    )
  }

  roi <- ftq_report_roi_row(roi_summary)
  targets <- as.character(rsc$node)
  target_connections <- rep(NA_real_, length(targets))
  local_connections <- NA_real_
  if (!is.null(roi)) {
    local_column <- "unique_local_streamlines"
    if (local_column %in% names(roi)) {
      local_connections <- as.numeric(roi[[local_column]][[1L]])
    }
    target_connections <- vapply(targets, function(target) {
      column <- paste0("unique_target_streamlines:", target)
      if (!column %in% names(roi)) return(NA_real_)
      as.numeric(roi[[column]][[1L]])
    }, numeric(1))
  }

  total_focal_voxels <- ftq_rsc_total_voxels(rsc)
  target_rows <- data.frame(
    cNode = targets,
    unique_streamline_connections = target_connections,
    connected_fNode_voxels = as.numeric(rsc$voxels_connected),
    rsc_percent = as.numeric(rsc$rsc_percent),
    total_fNode_voxels = total_focal_voxels,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  local_connected <- sum(as.numeric(result$local) > 0, na.rm = TRUE)
  local_total <- nrow(result)
  local_row <- data.frame(
    cNode = "local",
    unique_streamline_connections = local_connections,
    connected_fNode_voxels = local_connected,
    rsc_percent = if (local_total > 0L) 100 * local_connected / local_total else NA_real_,
    total_fNode_voxels = local_total,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out <- rbind(local_row, target_rows)
  attr(out, "unique_connections_available") <- !is.null(roi)
  out
}

ftq_report_fate_summary <- function(result, roi_summary = NULL) {
  roi <- ftq_report_roi_row(roi_summary)
  fates <- c("local", "leave_no_target", "leave_with_target")
  if (!is.null(roi)) {
    required <- c(
      "unique_focal_intersecting_streamlines",
      "unique_local_streamlines",
      "unique_leave_no_target_streamlines",
      "unique_leave_with_target_streamlines"
    )
    if (all(required %in% names(roi))) {
      denominator <- as.numeric(roi$unique_focal_intersecting_streamlines[[1L]])
      counts <- as.numeric(c(
        roi$unique_local_streamlines[[1L]],
        roi$unique_leave_no_target_streamlines[[1L]],
        roi$unique_leave_with_target_streamlines[[1L]]
      ))
      valid <- is.finite(denominator) && denominator >= 0 &&
        all(is.finite(counts) & counts >= 0) &&
        abs(sum(counts) - denominator) <= 1e-8
      if (!valid) {
        stop(
          "ROI-level unique streamline fates do not partition the unique ",
          "focal-intersecting streamline denominator."
        )
      }
      return(data.frame(
        fate = fates,
        count = counts,
        percent_of_denominator = if (denominator > 0) {
          100 * counts / denominator
        } else {
          rep(0, length(counts))
        },
        denominator_count = rep(denominator, length(counts)),
        summary_basis = "unique_focal_intersecting_streamlines",
        stringsAsFactors = FALSE
      ))
    }
  }

  association_fate <- ftq_fate_summary(result)
  denominator <- sum(as.numeric(association_fate$tracks))
  data.frame(
    fate = association_fate$fate,
    count = as.numeric(association_fate$tracks),
    percent_of_denominator = as.numeric(association_fate$percent),
    denominator_count = rep(denominator, nrow(association_fate)),
    summary_basis = "focal_voxel_streamline_associations",
    stringsAsFactors = FALSE
  )
}

ftq_write_maps <- function(result, template_image, map_dir, prefix, target_labels = NULL, overwrite = FALSE) {
  template <- RNifti::readNifti(template_image, internal = FALSE)
  dims <- dim(template)[1:3]
  coords <- ftq_coord_matrix(result$fNode_Coords)
  reserved <- c("fNode_Coords", "tracks", "local", "leave_no_target", "none", "cNode_Coords")
  if (is.null(target_labels)) target_labels <- setdiff(names(result), reserved)
  leave_with_target <- result$tracks - result$local - result$leave_no_target
  values <- c(
    list(
      tracks = result$tracks,
      local = result$local,
      "leave-no-target" = result$leave_no_target,
      "leave-with-target" = leave_with_target,
      "tracks-binary" = as.integer(result$tracks > 0)
    ),
    stats::setNames(as.list(result[target_labels]), paste0("target-", target_labels))
  )
  paths <- list()
  for (nm in names(values)) {
    out <- file.path(map_dir, paste0(prefix, "_fibrtraq_per-voxel-hits_", ftq_sanitize_label(nm), ".nii.gz"))
    ftq_check_overwrite(out, overwrite)
    arr <- array(0, dim = dims)
    arr[cbind(coords[,1], coords[,2], coords[,3])] <- as.numeric(values[[nm]])
    RNifti::writeNifti(arr, out, template = template, datatype = "float")
    paths[[nm]] <- out
  }
  paths
}

ftq_sanitize_label <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "-", x)
}

ftq_validate_maps <- function(result, map_paths) {
  coords <- ftq_coord_matrix(result$fNode_Coords)
  checks <- list()
  add <- function(name, passed, detail = "") {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name,
      passed = isTRUE(passed),
      detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  }
  expected <- list(
    tracks = result$tracks,
    local = result$local,
    "leave-no-target" = result$leave_no_target,
    "leave-with-target" = result$tracks - result$local - result$leave_no_target,
    "tracks-binary" = as.integer(result$tracks > 0)
  )
  reserved <- c("fNode_Coords", "tracks", "local", "leave_no_target", "none", "cNode_Coords")
  for (target in setdiff(names(result), reserved)) {
    expected[[paste0("target-", target)]] <- result[[target]]
  }
  for (nm in intersect(names(expected), names(map_paths))) {
    arr <- RNifti::readNifti(map_paths[[nm]], internal = FALSE)
    vals <- arr[cbind(coords[,1], coords[,2], coords[,3])]
    max_diff <- max(abs(vals - expected[[nm]]))
    sum_diff <- abs(sum(arr) - sum(expected[[nm]]))
    add(paste0("map_matches_csv_at_fnode:", nm), max_diff == 0, paste("max_diff =", max_diff))
    add(paste0("map_sum_matches_csv:", nm), sum_diff == 0, paste("sum_diff =", sum_diff))
  }
  checks
}

ftq_validate_expected <- function(result,
                                  rsc,
                                  map_paths,
                                  expected,
                                  tolerance,
                                  endpoint_summary = NULL,
                                  endpoint_penetration = NULL,
                                  endpoint_context = NULL) {
  checks <- list()
  add <- function(name, passed, detail = "") {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name,
      passed = isTRUE(passed),
      detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  }

  configured <- !is.null(expected$result_csv) || !is.null(expected$rsc_csv) ||
    !is.null(expected$endpoint_summary_csv) ||
    !is.null(expected$endpoint_penetration_csv) ||
    !is.null(expected$endpoint_context_csv) ||
    (!is.null(expected$map_paths) && length(expected$map_paths) > 0L)
  add("expected_outputs_configured", configured)

  if (!is.null(expected$result_csv)) {
    if (!file.exists(expected$result_csv)) {
      add("expected_result_file_exists", FALSE, expected$result_csv)
    } else {
      add("expected_result_file_exists", TRUE, expected$result_csv)
      expected_result <- ftq_read_table(expected$result_csv)
      checks <- c(checks, ftq_compare_tables(result, expected_result, "fNode_Coords", "expected_result", tolerance))
    }
  }

  if (!is.null(expected$rsc_csv)) {
    if (!file.exists(expected$rsc_csv)) {
      add("expected_rsc_file_exists", FALSE, expected$rsc_csv)
    } else {
      add("expected_rsc_file_exists", TRUE, expected$rsc_csv)
      expected_rsc <- ftq_read_table(expected$rsc_csv)
      checks <- c(checks, ftq_compare_tables(rsc, expected_rsc, "node", "expected_rsc", tolerance))
    }
  }

  if (!is.null(expected$map_paths) && length(expected$map_paths)) {
    checks <- c(checks, ftq_compare_expected_maps(map_paths, expected$map_paths, tolerance))
  }

  endpoint_specs <- list(
    endpoint_summary = list(
      actual = endpoint_summary,
      expected = expected$endpoint_summary_csv,
      key = c("target", "tolerance")
    ),
    endpoint_penetration = list(
      actual = endpoint_penetration,
      expected = expected$endpoint_penetration_csv,
      key = c("target", "unique_target_voxels")
    ),
    endpoint_context = list(
      actual = endpoint_context,
      expected = expected$endpoint_context_csv,
      key = "context"
    )
  )
  for (name in names(endpoint_specs)) {
    spec <- endpoint_specs[[name]]
    if (is.null(spec$expected)) next
    if (!file.exists(spec$expected)) {
      add(paste0("expected_", name, "_file_exists"), FALSE, spec$expected)
      next
    }
    add(paste0("expected_", name, "_file_exists"), TRUE, spec$expected)
    if (is.null(spec$actual)) {
      add(paste0("expected_", name, "_actual_supplied"), FALSE)
      next
    }
    actual <- ftq_read_table(spec$actual)
    expected_table <- ftq_read_table(spec$expected)
    if (length(spec$key) > 1L) {
      actual$.ftq_endpoint_key <- do.call(paste, c(actual[spec$key], sep = "\r"))
      expected_table$.ftq_endpoint_key <- do.call(
        paste,
        c(expected_table[spec$key], sep = "\r")
      )
      compare_key <- ".ftq_endpoint_key"
    } else {
      compare_key <- spec$key
    }
    checks <- c(
      checks,
      ftq_compare_tables(
        actual,
        expected_table,
        compare_key,
        paste0("expected_", name),
        tolerance
      )
    )
  }

  checks
}

ftq_compare_tables <- function(actual, expected, key, prefix, tolerance) {
  checks <- list()
  add <- function(name, passed, detail = "") {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name,
      passed = isTRUE(passed),
      detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  }

  column_ok <- setequal(names(actual), names(expected))
  missing_actual <- setdiff(names(expected), names(actual))
  missing_expected <- setdiff(names(actual), names(expected))
  detail <- paste(
    paste("missing_current =", paste(missing_actual, collapse = ",")),
    paste("missing_expected =", paste(missing_expected, collapse = ",")),
    sep = "; "
  )
  add(paste0(prefix, "_columns_match"), column_ok, detail)

  row_ok <- nrow(actual) == nrow(expected)
  add(paste0(prefix, "_row_count_matches"), row_ok,
      paste("current =", nrow(actual), "expected =", nrow(expected)))

  key_ok <- key %in% names(actual) && key %in% names(expected)
  add(paste0(prefix, "_key_column_present"), key_ok, key)
  if (!key_ok) {
    add(paste0(prefix, "_values_match"), FALSE, "Key column missing.")
    return(checks)
  }

  actual_key <- as.character(actual[[key]])
  expected_key <- as.character(expected[[key]])
  unique_ok <- !anyDuplicated(actual_key) && !anyDuplicated(expected_key)
  add(paste0(prefix, "_keys_unique"), unique_ok)
  keys_ok <- setequal(actual_key, expected_key)
  add(paste0(prefix, "_keys_match"), keys_ok,
      paste("missing_current =", paste(setdiff(expected_key, actual_key), collapse = ","),
            "missing_expected =", paste(setdiff(actual_key, expected_key), collapse = ",")))

  if (!(column_ok && row_ok && unique_ok && keys_ok)) {
    add(paste0(prefix, "_values_match"), FALSE, "Column, row, or key mismatch prevents value comparison.")
    return(checks)
  }

  common <- names(actual)
  ordered_keys <- sort(actual_key)
  actual <- actual[match(ordered_keys, actual_key), common, drop = FALSE]
  expected <- expected[match(ordered_keys, expected_key), common, drop = FALSE]

  differing <- character()
  max_diff <- 0
  for (nm in setdiff(common, key)) {
    a <- actual[[nm]]
    e <- expected[[nm]]
    if (ftq_numeric_like(a) && ftq_numeric_like(e)) {
      a_num <- suppressWarnings(as.numeric(a))
      e_num <- suppressWarnings(as.numeric(e))
      na_ok <- identical(is.na(a_num), is.na(e_num))
      diff <- abs(a_num - e_num)
      diff <- diff[is.finite(diff)]
      col_max <- if (length(diff)) max(diff) else 0
      max_diff <- max(max_diff, col_max)
      if (!na_ok || col_max > tolerance) differing <- c(differing, nm)
    } else if (!identical(as.character(a), as.character(e))) {
      differing <- c(differing, nm)
    }
  }
  add(paste0(prefix, "_values_match"), length(differing) == 0L,
      paste("max_abs_diff =", signif(max_diff, 8),
            "differing_columns =", paste(differing, collapse = ",")))
  checks
}

ftq_numeric_like <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(TRUE)
  y <- suppressWarnings(as.numeric(x))
  all(is.na(x) == is.na(y))
}

ftq_compare_expected_maps <- function(current_map_paths, expected_map_paths, tolerance) {
  checks <- list()
  add <- function(name, passed, detail = "") {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name,
      passed = isTRUE(passed),
      detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  }
  if (is.null(current_map_paths)) current_map_paths <- list()
  missing_current <- setdiff(names(expected_map_paths), names(current_map_paths))
  add("expected_map_keys_available", length(missing_current) == 0L,
      paste(missing_current, collapse = ","))
  for (nm in intersect(names(expected_map_paths), names(current_map_paths))) {
    expected_path <- expected_map_paths[[nm]]
    current_path <- current_map_paths[[nm]]
    expected_exists <- file.exists(expected_path)
    current_exists <- file.exists(current_path)
    add(paste0("expected_map_file_exists:", nm), expected_exists, expected_path)
    add(paste0("current_map_file_exists:", nm), current_exists, current_path)
    if (!expected_exists || !current_exists) next
    expected_arr <- as.array(RNifti::readNifti(expected_path, internal = FALSE))
    current_arr <- as.array(RNifti::readNifti(current_path, internal = FALSE))
    dims_ok <- identical(dim(current_arr), dim(expected_arr))
    add(paste0("expected_map_dims_match:", nm), dims_ok,
        paste("current =", paste(dim(current_arr), collapse = "x"),
              "expected =", paste(dim(expected_arr), collapse = "x")))
    if (!dims_ok) next
    diff <- abs(current_arr - expected_arr)
    diff <- diff[is.finite(diff)]
    max_diff <- if (length(diff)) max(diff) else 0
    add(paste0("expected_map_values_match:", nm), max_diff <= tolerance,
        paste("max_abs_diff =", signif(max_diff, 8)))
  }
  checks
}

ftq_write_top_tables <- function(result, table_dir, prefix, n = NULL) {
  if (!is.null(n) && (!is.numeric(n) || length(n) != 1L ||
      is.na(n) || !is.finite(n) || n < 0 || n != floor(n))) {
    stop("n must be NULL for all voxels or a single nonnegative integer.")
  }
  view <- ftq_report_voxel_export(result)
  xyz <- view$data
  coord_names <- view$columns$column[view$columns$role == "coordinate"]
  derived_name <- "leave_with_target"
  while (derived_name %in% names(xyz)) derived_name <- paste0("derived_", derived_name)
  xyz[[derived_name]] <- xyz$tracks - xyz$local - xyz$leave_no_target
  target_cols <- view$columns$column[view$columns$role == "target"]
  out <- list()
  metrics <- c("leave_no_target", target_cols)
  for (metric in metrics) {
    ord <- order(-xyz[[metric]], xyz[[coord_names[1L]]],
                 xyz[[coord_names[2L]]], xyz[[coord_names[3L]]])
    tab <- xyz[ord, unique(c(coord_names, "tracks", "local", "leave_no_target", derived_name, metric)), drop = FALSE]
    if (!is.null(n)) tab <- tab[seq_len(min(n, nrow(tab))), , drop = FALSE]
    path <- ftq_ranked_table_path(table_dir, prefix, metric, n)
    utils::write.csv(tab, path, row.names = FALSE)
    out[[metric]] <- path
  }
  out
}

ftq_plot_fate_summary <- function(result, path, fate = NULL) {
  if (is.null(fate)) fate <- ftq_report_fate_summary(result)
  basis <- unique(as.character(fate$summary_basis))
  unique_basis <- length(basis) == 1L &&
    identical(basis, "unique_focal_intersecting_streamlines")
  grDevices::png(ftq_png_filename(path), width = 2400, height = 1400, res = 240)
  old <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old)
    grDevices::dev.off()
  }, add = TRUE)
  if (!nrow(fate) || all(!is.finite(fate$denominator_count) | fate$denominator_count <= 0)) {
    graphics::plot.new()
    graphics::title(main = "Streamline Fate: Unevaluable")
    graphics::text(0.5, 0.55, "No positive denominator; percentages are undefined.", cex = 1.15)
    graphics::text(0.5, 0.43, "See the complete fate table for counts and denominator.", cex = 1)
    return(invisible(path))
  }
  cols <- c(local = "#1b9e77", leave_no_target = "#d95f02", leave_with_target = "#7570b3")
  labels <- sprintf(
    "%s: %.1f%% (n=%s)",
    gsub("_", " ", fate$fate),
    fate$percent_of_denominator,
    format(fate$count, big.mark = ",")
  )
  graphics::layout(matrix(c(1, 2), nrow = 2L), heights = c(1, 5))
  graphics::par(oma = c(0, 0, 3, 0), mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::legend(
    "center",
    fill = cols[fate$fate],
    legend = labels,
    border = NA,
    horiz = TRUE,
    bty = "n",
    cex = 0.85
  )
  graphics::par(mar = c(5, 6, 1, 2), xaxs = "i")
  bp <- graphics::barplot(
    matrix(fate$percent_of_denominator, ncol = 1L),
    horiz = TRUE,
    col = cols[fate$fate],
    border = NA,
    xlim = c(0, 100),
    names.arg = "",
    xlab = if (unique_basis) {
      "Percent of unique focal-intersecting streamlines"
    } else {
      "Percent of focal-voxel streamline associations"
    }
  )
  left <- c(0, cumsum(utils::head(fate$percent_of_denominator, -1)))
  centers <- left + fate$percent_of_denominator / 2
  show_inside <- fate$percent_of_denominator >= 12
  graphics::text(
    centers[show_inside],
    rep(bp, sum(show_inside)),
    labels = sprintf("%.1f%%\n(n=%s)",
                     fate$percent_of_denominator[show_inside],
                     format(fate$count[show_inside], big.mark = ",")),
    cex = 0.8
  )
  graphics::mtext(
    if (unique_basis) "Unique Streamline Fate Composition" else "Association-Weighted Fate Composition",
    side = 3,
    outer = TRUE,
    line = 0.5,
    cex = 1.25,
    font = 2
  )
}

ftq_plot_rsc_summary <- function(rsc, path, connection_summary = NULL) {
  if (is.null(connection_summary)) {
    connection_summary <- data.frame(
      cNode = as.character(rsc$node),
      unique_streamline_connections = NA_real_,
      connected_fNode_voxels = as.numeric(rsc$voxels_connected),
      rsc_percent = as.numeric(rsc$rsc_percent),
      total_fNode_voxels = ftq_rsc_total_voxels(rsc),
      stringsAsFactors = FALSE
    )
  }
  grDevices::png(ftq_png_filename(path), width = 2000, height = 1400, res = 240)
  old <- graphics::par(mar = c(6, 6, 4, 2))
  on.exit({
    graphics::par(old)
    grDevices::dev.off()
  }, add = TRUE)
  local_idx <- which(connection_summary$cNode == "local")
  target_idx <- which(connection_summary$cNode != "local")
  target_idx <- target_idx[
    order(connection_summary$rsc_percent[target_idx], decreasing = TRUE)
  ]
  ord <- c(local_idx, target_idx)
  rr <- connection_summary[ord, , drop = FALSE]
  cols <- grDevices::hcl.colors(nrow(rr), "Dark 3")
  cols[rr$cNode == "local"] <- "#16827A"
  bp <- graphics::barplot(
    rr$rsc_percent,
    names.arg = rr$cNode,
    col = cols,
    border = NA,
    ylab = "Focal voxels connected (%)",
    main = "Regional Streamline Connectivity (RSC)",
    ylim = c(0, max(1, rr$rsc_percent) * 1.25)
  )
  count_labels <- ifelse(
    is.finite(rr$unique_streamline_connections),
    paste0("n=", format(rr$unique_streamline_connections, big.mark = ",")),
    "unique n unavailable"
  )
  graphics::text(
    bp,
    rr$rsc_percent,
    labels = sprintf("%.2f%%\n%s", rr$rsc_percent, count_labels),
    pos = 3,
    cex = 0.7
  )
}

ftq_plot_tracks_distribution <- function(result, path) {
  x <- result$tracks[result$tracks > 0]
  grDevices::png(ftq_png_filename(path), width = 2400, height = 1400, res = 240)
  old <- graphics::par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
  on.exit({
    graphics::par(old)
    grDevices::dev.off()
  }, add = TRUE)
  if (length(x)) {
    graphics::hist(log1p(x), breaks = 40, col = "grey75", border = "white",
                   main = "Per-Voxel Streamline Count Distribution",
                   xlab = "log1p(streamlines per occupied fNode voxel)")
    graphics::plot(stats::ecdf(x), main = "Per-Voxel Streamline Count ECDF",
                   xlab = "streamlines per occupied fNode voxel", ylab = "ECDF")
  } else {
    graphics::plot.new()
    graphics::title("Per-Voxel Streamline Count Distribution")
    graphics::text(0.5, 0.5, "No occupied fNode voxels")
    graphics::plot.new()
    graphics::title("Per-Voxel Streamline Count ECDF")
    graphics::text(0.5, 0.5, "No occupied fNode voxels")
  }
}

ftq_plot_spatial_bundle <- function(result,
                                    map_paths,
                                    figs,
                                    background_image = NULL,
                                    tractogram_context = NULL,
                                    context_source = "auto",
                                    target_intersection_mode = "full_streamline",
                                    local_search_radius_mm = 2,
                                    display_quantile = 0.99,
                                    report_slices = NULL) {
  if (is.null(map_paths) || !length(map_paths)) return(invisible(NULL))
  assignment_label <- ftq_report_assignment_label(
    target_intersection_mode,
    local_search_radius_mm
  )
  if (!is.null(map_paths$tracks)) {
    arr <- RNifti::readNifti(map_paths$tracks, internal = FALSE)
    bg <- ftq_report_background(background_image, tractogram_context, dim(arr)[1:3], context_source)
    ftq_plot_montage(
      arr,
      figs$spatial_tracks,
      "fNode Streamline Density Heatmap",
      background = bg$image,
      palette = "Inferno",
      legend_label = "Streamline density",
      context = TRUE,
      background_label = bg$label,
      display_quantile = display_quantile,
      slices = report_slices
    )
  }
  if (!is.null(map_paths[["leave-no-target"]])) {
    arr <- RNifti::readNifti(map_paths[["leave-no-target"]], internal = FALSE)
    bg <- ftq_report_background(background_image, tractogram_context, dim(arr)[1:3], context_source)
    ftq_plot_montage(
      arr,
      figs$spatial_leave_no_target,
      paste0("fNode Density: No Assigned cNode\n(", assignment_label, ")"),
      background = bg$image,
      palette = "Rocket",
      legend_label = "Streamline density",
      context = TRUE,
      background_label = bg$label,
      display_quantile = display_quantile,
      slices = report_slices
    )
  }
  target_names <- grep("^target-", names(map_paths), value = TRUE)
  connection_names <- c(
    if (!is.null(map_paths$local)) "local" else character(),
    target_names
  )
  if (length(connection_names)) {
    arrays <- lapply(map_paths[connection_names], RNifti::readNifti, internal = FALSE)
    bg <- ftq_report_background(background_image, tractogram_context, dim(arrays[[1]])[1:3], context_source)
    ftq_plot_multi_map_montage(
      arrays,
      figs$spatial_targets,
      paste0("fNode Connection Densities\n(", assignment_label, ")"),
      background = bg$image,
      background_label = bg$label,
      display_quantile = display_quantile,
      slices = report_slices
    )
  }
}

ftq_report_assignment_label <- function(target_intersection_mode,
                                        local_search_radius_mm = 2) {
  mode <- ftq_normalize_target_intersection_mode(target_intersection_mode)
  switch(
    mode,
    full_streamline = "full-streamline intersections",
    endpoint_voxels = "exact end voxels",
    local_search = paste0(
      "local search, ",
      format(as.numeric(local_search_radius_mm), trim = TRUE),
      " mm"
    )
  )
}

ftq_plot_montage <- function(arr,
                             path,
                             title,
                             nslices = 9L,
                             background = NULL,
                             palette = "Inferno",
                             legend_label = "Count",
                             context = FALSE,
                             background_label = NULL,
                             display_quantile = 1,
                             slices = NULL,
                             validity_mask = NULL) {
  a <- as.array(arr)
  support <- a
  if (!is.null(validity_mask)) {
    if (!identical(dim(validity_mask), dim(a)) || anyNA(validity_mask)) {
      stop("The montage validity mask must match the image grid and contain no missing values.")
    }
    support <- array(as.numeric(validity_mask), dim = dim(a))
    a[!validity_mask] <- NA_real_
  }
  slices <- ftq_report_slices(slices, dim(a)[3L])
  if (is.null(slices)) {
    zsum <- apply(support, 3, sum)
    slices <- which(zsum > 0)
    if (!length(slices)) slices <- seq_len(dim(a)[3])
    if (length(slices) > nslices) {
      slices <- unique(round(seq(min(slices), max(slices), length.out = nslices)))
    }
  }
  crop <- ftq_crop_ranges(list(support), margin = 3L)
  xr <- crop$x
  yr <- crop$y
  use_context <- isTRUE(context) && !is.null(background)
  slices_per_row <- min(3L, length(slices))
  nr <- ceiling(length(slices) / slices_per_row)
  panels_per_row <- if (use_context) slices_per_row * 2L else slices_per_row
  png_width <- if (use_context) 1120L * slices_per_row + 340L else 760L * slices_per_row + 300L
  png_height <- max(900L, 620L * nr)
  grDevices::png(ftq_png_filename(path), width = png_width, height = png_height, res = 240)
  layout_mat <- cbind(
    matrix(seq_len(nr * panels_per_row), nrow = nr, byrow = TRUE),
    rep(nr * panels_per_row + 1L, nr)
  )
  panel_widths <- if (use_context) {
    rep(c(0.78, 1.08), slices_per_row)
  } else {
    rep(1, slices_per_row)
  }
  panel_widths <- c(panel_widths, 0.42)
  graphics::layout(layout_mat, widths = panel_widths, heights = rep(1, nr))
  old <- graphics::par(mar = c(1.5, 1.5, 3, 1), oma = c(2.2, 0, 4.5, 0))
  on.exit({
    graphics::par(old)
    graphics::layout(1)
    grDevices::dev.off()
  }, add = TRUE)
  pal <- ftq_overlay_palette(palette)
  zlim <- ftq_positive_zlim(a, quantile = display_quantile)
  for (z in slices) {
    if (use_context) {
      ftq_draw_overlay_slice(
        a[, , z],
        seq_len(dim(a)[1]),
        seq_len(dim(a)[2]),
        zlim,
        pal,
        paste("z =", z, "full slice"),
        background[, , z],
        crop_box = crop,
        show_zeros = !is.null(validity_mask)
      )
      ftq_draw_overlay_slice(
        a[xr, yr, z],
        xr,
        yr,
        zlim,
        pal,
        paste("z =", z, "focal zoom"),
        background[xr, yr, z],
        show_zeros = !is.null(validity_mask)
      )
    } else {
      mat <- a[xr, yr, z]
      bg_mat <- if (!is.null(background)) background[xr, yr, z] else NULL
      ftq_draw_overlay_slice(mat, xr, yr, zlim, pal, paste("z =", z), bg_mat,
                             show_zeros = !is.null(validity_mask))
    }
  }
  filled_panels <- length(slices) * if (use_context) 2L else 1L
  if (filled_panels < nr * panels_per_row) {
    for (i in seq_len(nr * panels_per_row - filled_panels)) graphics::plot.new()
  }
  ftq_draw_colorbar(
    zlim,
    pal,
    legend_label
  )
  graphics::mtext(title, outer = TRUE, cex = 1.2, font = 2)
  footer <- c(
    if (display_quantile < 1) {
      "Display capped at the 99th percentile of positive voxels; raw map values are unchanged."
    } else {
      character()
    },
    if (!is.null(background_label) && nzchar(background_label)) {
      paste("Background:", background_label)
    } else {
      character()
    }
  )
  if (length(footer)) {
    graphics::mtext(paste(footer, collapse = "  "), side = 1, outer = TRUE, cex = 0.62)
  }
}

ftq_plot_multi_map_montage <- function(arrays,
                                       path,
                                       title,
                                       background = NULL,
                                       background_label = NULL,
                                       display_quantile = 1,
                                       slices = NULL) {
  n <- length(arrays)
  selected_slices <- ftq_report_slices(slices)
  for (arr in arrays) ftq_report_slices(selected_slices, dim(arr)[3L])
  slice_counts <- vapply(arrays, function(x) {
    count <- sum(apply(as.array(x), 3, sum) > 0)
    max(1L, count)
  }, integer(1))
  panel_count <- if (is.null(selected_slices)) min(3L, max(slice_counts)) else length(selected_slices)
  png_width <- max(1460L, 760L * panel_count + 340L)
  grDevices::png(ftq_png_filename(path), width = png_width, height = max(1200, 700 * n), res = 240)
  layout_mat <- matrix(
    seq_len(n * (panel_count + 1L)),
    nrow = n,
    ncol = panel_count + 1L,
    byrow = TRUE
  )
  graphics::layout(
    layout_mat,
    widths = c(rep(1, panel_count), 0.42),
    heights = rep(1, n)
  )
  old <- graphics::par(mar = c(1.5, 1.5, 3, 1), oma = c(2.2, 0, 4.5, 0))
  on.exit({
    graphics::par(old)
    graphics::layout(1)
    grDevices::dev.off()
  }, add = TRUE)
  crop <- ftq_crop_ranges(arrays, margin = 3L)
  palette_names <- c("Inferno", "Viridis", "Plasma", "Mako", "Rocket", "Cividis")
  for (i in seq_along(arrays)) {
    nm <- names(arrays)[[i]]
    a <- as.array(arrays[[nm]])
    pal <- ftq_overlay_palette(palette_names[[1L + ((i - 1L) %% length(palette_names))]])
    zlim <- ftq_positive_zlim(a, quantile = display_quantile)
    slices <- selected_slices
    if (is.null(slices)) {
      zsum <- apply(a, 3, sum)
      slices <- which(zsum > 0)
      if (!length(slices)) slices <- max(1L, round(dim(a)[3] / 2))
      slices <- unique(round(seq(
        min(slices), max(slices), length.out = panel_count
      )))
    }
    for (z in slices) {
      mat <- a[crop$x, crop$y, z]
      bg_mat <- if (!is.null(background)) background[crop$x, crop$y, z] else NULL
      label <- if (identical(nm, "local")) "local" else gsub("^target-", "", nm)
      ftq_draw_overlay_slice(
        mat,
        crop$x,
        crop$y,
        zlim,
        pal,
        paste(label, "z =", z),
        bg_mat
      )
    }
    if (length(slices) < panel_count) {
      for (unused in seq_len(panel_count - length(slices))) graphics::plot.new()
    }
    ftq_draw_colorbar(
      zlim,
      pal,
      "Streamline density"
    )
  }
  graphics::mtext(title, outer = TRUE, cex = 1.2, font = 2)
  footer <- c(
    if (display_quantile < 1) {
      "Independent 99th-percentile scales; compare patterns within rows. Raw values are unchanged."
    } else {
      character()
    },
    if (!is.null(background_label) && nzchar(background_label)) {
      paste("Background:", background_label)
    } else {
      character()
    }
  )
  if (length(footer)) {
    graphics::mtext(paste(footer, collapse = "  "), side = 1, outer = TRUE, cex = 0.62)
  }
}

ftq_report_background <- function(background_image, tractogram_context, dims, context_source = "auto") {
  context_image <- ftq_context_density_image(tractogram_context)
  if (identical(context_source, "tractogram")) {
    bg <- ftq_read_context_background(context_image, dims)
    if (!is.null(bg)) {
      return(list(image = bg, label = "Tractogram-derived density context"))
    }
  }
  bg <- ftq_read_background(background_image, dims)
  if (!is.null(bg)) {
    return(list(image = bg, label = "anatomical image"))
  }
  bg <- ftq_read_context_background(context_image, dims)
  if (!is.null(bg)) {
    return(list(image = bg, label = "Tractogram-derived density context"))
  }
  list(image = NULL, label = NULL)
}

ftq_read_background <- function(path, dims) {
  if (is.null(path)) return(NULL)
  if (!file.exists(path)) {
    warning("Background image does not exist; plotting without background.")
    return(NULL)
  }
  bg <- as.array(ftq_read_volume_data(path))
  if (length(dim(bg)) > 3L) bg <- bg[, , , 1]
  if (!identical(dim(bg)[1:3], as.integer(dims))) {
    warning("Background image dimensions do not match map dimensions; plotting without background.")
    return(NULL)
  }
  bg
}

ftq_read_context_background <- function(path, dims) {
  bg <- ftq_read_background(path, dims)
  if (is.null(bg)) return(NULL)
  bg[!is.finite(bg) | bg < 0] <- 0
  density <- log1p(bg)
  vals <- density[density > 0 & is.finite(density)]
  if (!length(vals)) return(density)
  upper <- stats::quantile(vals, 0.98, names = FALSE, na.rm = TRUE)
  if (!is.finite(upper) || upper <= 0) upper <- max(vals)
  density <- pmin(density / upper, 1)
  coverage <- ftq_dilate_mask(bg > 0, iterations = 2L)
  pmax(density, 0.18 * coverage)
}

ftq_dilate_mask <- function(mask, iterations = 1L) {
  out <- array(as.logical(mask), dim = dim(mask))
  dims <- dim(out)
  iterations <- as.integer(iterations)[1]
  if (is.na(iterations) || iterations <= 0L) return(out)
  for (i in seq_len(iterations)) {
    nxt <- out
    shifted <- array(FALSE, dim = dims)
    shifted[-1L, , ] <- out[-dims[1], , ]
    nxt <- nxt | shifted
    shifted[] <- FALSE
    shifted[-dims[1], , ] <- out[-1L, , ]
    nxt <- nxt | shifted
    shifted[] <- FALSE
    shifted[, -1L, ] <- out[, -dims[2], ]
    nxt <- nxt | shifted
    shifted[] <- FALSE
    shifted[, -dims[2], ] <- out[, -1L, ]
    nxt <- nxt | shifted
    shifted[] <- FALSE
    shifted[, , -1L] <- out[, , -dims[3]]
    nxt <- nxt | shifted
    shifted[] <- FALSE
    shifted[, , -dims[3]] <- out[, , -1L]
    out <- nxt | shifted
  }
  out
}

ftq_crop_ranges <- function(arrays, margin = 3L) {
  dims <- dim(as.array(arrays[[1]]))[1:3]
  nz_list <- lapply(arrays, function(x) which(as.array(x) > 0, arr.ind = TRUE))
  nz_list <- nz_list[vapply(nz_list, nrow, integer(1)) > 0L]
  if (!length(nz_list)) {
    return(list(x = seq_len(dims[1]), y = seq_len(dims[2])))
  }
  nz <- do.call(rbind, nz_list)
  xr <- range(nz[, 1])
  yr <- range(nz[, 2])
  list(
    x = max(1L, xr[1] - margin):min(dims[1], xr[2] + margin),
    y = max(1L, yr[1] - margin):min(dims[2], yr[2] + margin)
  )
}

ftq_positive_zlim <- function(arr, quantile = 1) {
  if (is.list(arr) && !is.array(arr)) {
    vals <- unlist(lapply(arr, function(x) as.numeric(as.array(x))), use.names = FALSE)
  } else {
    vals <- as.numeric(arr)
  }
  vals <- vals[is.finite(vals) & vals > 0]
  if (!length(vals)) return(c(0, 1))
  quantile <- as.numeric(quantile)[1L]
  if (!is.finite(quantile) || quantile <= 0 || quantile > 1) {
    stop("quantile must be one finite number in (0, 1].")
  }
  upper <- if (quantile < 1) {
    stats::quantile(vals, quantile, names = FALSE, na.rm = TRUE, type = 8)
  } else {
    max(vals)
  }
  c(0, max(1, upper))
}

ftq_overlay_palette <- function(palette, n = 64L) {
  raw <- grDevices::hcl.colors(n + 12L, palette)
  grDevices::adjustcolor(raw[seq.int(13L, n + 12L)], alpha.f = 0.88)
}

ftq_draw_overlay_slice <- function(map_mat, x, y, zlim, pal, main, bg_mat = NULL, crop_box = NULL,
                                  show_zeros = FALSE) {
  if (!is.null(bg_mat)) {
    graphics::image(x, y, ftq_normalize_background(bg_mat), col = grDevices::gray.colors(128, 0, 0.95),
                    zlim = c(0, 1), axes = FALSE, asp = 1, main = main, xlab = "", ylab = "")
  } else {
    blank <- matrix(0, nrow = nrow(map_mat), ncol = ncol(map_mat))
    graphics::image(x, y, blank, col = "#f8fafc", zlim = c(0, 1), axes = FALSE, asp = 1,
                    main = main, xlab = "", ylab = "")
  }
  overlay <- map_mat
  overlay[!is.finite(overlay) | if (show_zeros) overlay < 0 else overlay <= 0] <- NA
  if (any(!is.na(overlay))) {
    graphics::image(x, y, overlay, col = pal, zlim = zlim, axes = FALSE, asp = 1,
                    add = TRUE, xlab = "", ylab = "")
  }
  if (!is.null(crop_box)) {
    graphics::rect(min(crop_box$x), min(crop_box$y), max(crop_box$x), max(crop_box$y),
                   border = "#1f9bd1", lwd = 2)
  }
  graphics::box()
}

ftq_normalize_background <- function(mat) {
  vals <- as.numeric(mat)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(matrix(0, nrow = nrow(mat), ncol = ncol(mat)))
  q <- stats::quantile(vals, c(0.02, 0.98), names = FALSE, na.rm = TRUE)
  if (q[1] == q[2]) q <- range(vals)
  if (q[1] == q[2]) return(matrix(0, nrow = nrow(mat), ncol = ncol(mat)))
  out <- pmax(pmin(mat, q[2]), q[1])
  (out - q[1]) / (q[2] - q[1])
}

ftq_draw_colorbar <- function(zlim, pal, label) {
  old <- graphics::par(mar = c(1.2, 0.1, 2.4, 3.6), xaxs = "i", yaxs = "i")
  on.exit(graphics::par(old), add = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = zlim)
  bar <- grDevices::as.raster(matrix(rev(pal), ncol = 1L))
  graphics::rasterImage(bar, 0.15, zlim[1], 0.55, zlim[2])
  graphics::rect(0.15, zlim[1], 0.55, zlim[2], border = "black", lwd = 0.8)
  ticks <- pretty(zlim, n = 5)
  ticks <- sort(unique(c(zlim[1], ticks[ticks >= zlim[1] & ticks <= zlim[2]], zlim[2])))
  graphics::axis(
    4,
    at = ticks,
    labels = format(signif(ticks, 3L), big.mark = ",", trim = TRUE, scientific = FALSE),
    pos = 0.55,
    las = 1,
    cex.axis = 0.65,
    tcl = -0.25,
    lwd = 0,
    lwd.ticks = 0.8
  )
  graphics::mtext(label, side = 4, line = 1.7, cex = 0.7)
}

ftq_report_preflight_summary <- function(preflight) {
  if (is.null(preflight) || is.null(preflight$checks) || !nrow(preflight$checks)) {
    return(data.frame())
  }
  checks <- preflight$checks
  groups <- list(
    "Required files" =
      "^(fnode|tractogram|local|ref_image|target:|endpoint_context:).*_exists$",
    "Node masks" =
      "^(fnode|local|target:).*(readable|is_3d|nonempty|binary_values)$",
    "Image grids" =
      "(dims_match_fnode|affine_matches_fnode)$",
    "Coordinate and assignment settings" = paste0(
      "^(track_space_supported|discretization_supported|focal_assignment_mode_supported|",
      "streamline_voxelization_mode_supported|target_intersection_mode_supported|",
      "local_search_radius_mm_positive|local_search_geometry_available|",
      "local_search_spatial_units_known)$"
    ),
    "Tractogram reader" =
      "^(tractogram_format_supported|tck_header_readable|tck_world_coordinate_mapping|tck_spatial_geometry:.*|track_reader_supported|track_chunk_size_positive)$"
  )
  rows <- lapply(names(groups), function(label) {
    selected <- grepl(groups[[label]], checks$check, perl = TRUE)
    if (!any(selected)) return(NULL)
    part <- checks[selected, , drop = FALSE]
    failed <- part[!part$passed, , drop = FALSE]
    failed_errors <- failed[failed$severity == "error", , drop = FALSE]
    status <- if (nrow(failed_errors)) {
      "FAIL"
    } else if (nrow(failed)) {
      "REVIEW"
    } else {
      "PASS"
    }
    detail <- paste0(sum(part$passed), "/", nrow(part), " checks passed")
    if (nrow(failed)) {
      detail <- paste0(
        detail,
        "; review: ",
        paste(failed$check, collapse = ", ")
      )
    }
    data.frame(section = label, status = status, detail = detail,
               stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

ftq_write_html_report <- function(path,
                                  result,
                                  rsc,
                                  manifest,
                                  validation,
                                  preflight,
                                  figs,
                                  top_tables,
                                  tractogram_context = NULL,
                                  endpoint_data = NULL,
                                  local_search_data = NULL,
                                  connection_summary = NULL,
                                  fate_summary = NULL,
                                  source_paths = NULL,
                                  map_paths = NULL,
                                  metadata_path = NULL,
                                  guide = NULL) {
  if (is.null(connection_summary)) {
    connection_summary <- ftq_report_connection_summary(result, rsc)
  }
  if (is.null(fate_summary)) fate_summary <- ftq_report_fate_summary(result)
  esc <- function(x) {
    x <- as.character(x)
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x <- gsub('"', "&quot;", x, fixed = TRUE)
    x <- gsub("'", "&#39;", x, fixed = TRUE)
    x
  }
  table_html <- function(df) {
    if (is.null(df) || !nrow(df)) return("<p class='note'>No rows available.</p>")
    header <- paste0("<tr>", paste0("<th>", esc(names(df)), "</th>", collapse = ""), "</tr>")
    rows <- apply(df, 1, function(r) paste0("<tr>", paste0("<td>", esc(r), "</td>", collapse = ""), "</tr>"))
    paste0(
      "<div class='table-scroll' role='region' tabindex='0'>",
      "<table>", header, paste(rows, collapse = "\n"), "</table>",
      "</div>"
    )
  }
  fig_tag <- function(p, title, caption = title) {
    rp <- ftq_report_href(p, dirname(path))
    paste0(
      "<figure><img src='", esc(rp), "' alt='", esc(title), "'>",
      "<figcaption><b>", esc(title), ".</b> ", esc(caption),
      " <a href='", esc(rp), "'>Open or save image (PNG)</a>",
      "</figcaption></figure>"
    )
  }
  css <- paste0(
    "<style>",
    "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;",
    "margin:32px;line-height:1.45;color:#222;overflow-wrap:anywhere}",
    "h1,h2{letter-spacing:0}",
    "h2{margin-top:36px}h3{margin-top:28px}",
    ".table-scroll{width:100%;max-width:100%;overflow-x:auto;",
    "-webkit-overflow-scrolling:touch;margin:12px 0 24px}",
    ".table-scroll:focus{outline:2px solid #2563eb;outline-offset:2px}",
    "table{border-collapse:collapse;min-width:100%;margin:0}",
    "th,td{border:1px solid #ddd;padding:6px 10px;text-align:right;",
    "vertical-align:top;min-width:72px}",
    "th{background:#f3f4f6;min-width:96px;max-width:260px;",
    "overflow-wrap:anywhere}",
    "th:first-child,td:first-child{text-align:left}",
    "img{max-width:100%;height:auto;border:1px solid #eee}",
    "figure{margin:20px 0 36px}",
    "figcaption{max-width:1050px;margin-top:8px;color:#475569;font-size:.94rem}",
    "details{border-top:1px solid #d8dee6;border-bottom:1px solid #d8dee6;",
    "padding:10px 0;margin:18px 0}",
    "summary{cursor:pointer;font-weight:700;color:#334155}",
    "details .table-scroll{margin-top:14px}",
    ".pass{color:#176b35;font-weight:700}",
    ".fail{color:#9f1239;font-weight:700}",
    ".warn{color:#92400e;font-weight:700}",
    ".note,.meta{background:#f8fafc;padding:12px 16px;",
    "border-left:4px solid #64748b}",
    ".note{border-left-color:#94a3b8}",
    "@media(max-width:700px){body{margin:16px}",
    "h1{font-size:1.6rem}h2{font-size:1.3rem}",
    "th,td{padding:5px 7px;font-size:.85rem}}",
    "</style>"
  )
  identity <- unlist(list(manifest$subject, manifest$focal_region), use.names = FALSE)
  identity <- identity[!is.na(identity) & nzchar(identity)]
  title <- if (length(identity)) paste("FibrTraQ:", paste(identity, collapse = " / ")) else "FibrTraQ Report"
  preflight_html <- ""
  if (!is.null(preflight)) {
    cls <- if (isTRUE(preflight$passed)) "pass" else "fail"
    failed_warnings <- sum(!preflight$checks$passed & preflight$checks$severity %in% c("warning", "recommendation"))
    preflight_summary <- ftq_report_preflight_summary(preflight)
    preflight_html <- paste0(
      "<h2>Preflight</h2><p class='", cls, "'>",
      if (isTRUE(preflight$passed)) "PASS" else "FAIL",
      if (failed_warnings > 0L) {
        paste0(" <span class='warn'>(", failed_warnings, " non-blocking note(s))</span>")
      } else {
        ""
      },
      "</p>",
      if (any(preflight_summary$status != "PASS")) table_html(preflight_summary) else "",
      "<details><summary>Full preflight details (", nrow(preflight$checks),
      " checks)</summary>", table_html(preflight$checks), "</details>"
    )
  }
  validation_html <- ""
  if (!is.null(validation)) {
    cls <- if (isTRUE(validation$passed)) "pass" else "fail"
    biological_status <- ftq_value(validation$biological_qc_status, "unavailable")
    biological_class <- if (identical(biological_status, "pass")) "pass" else "warn"
    biological_warnings <- ftq_value(validation$biological_qc_warnings, character())
    warning_html <- if (length(biological_warnings)) {
      paste0(
        "<ul>",
        paste0("<li>", esc(biological_warnings), "</li>", collapse = ""),
        "</ul>"
      )
    } else {
      ""
    }
    checks <- validation$checks
    passed_checks <- if (is.null(checks)) 0L else sum(checks$passed)
    total_checks <- if (is.null(checks)) 0L else nrow(checks)
    validation_html <- paste0(
      "<h2>Validation</h2><p class='", cls, "'>",
      passed_checks, " of ", total_checks,
      " automated computational checks passed. Overall status: ",
      if (isTRUE(validation$passed)) "PASS" else "FAIL",
      ".</p><p><b>Biological QC:</b> <span class='",
      biological_class, "'>", esc(toupper(biological_status)),
      "</span></p>", warning_html,
      if (total_checks > 0L) {
        paste0(
          "<details><summary>Full validation details (", total_checks,
          " checks)</summary>", table_html(checks), "</details>"
        )
      } else {
        ""
      }
    )
  }
  target_assignment_mode <- ftq_value(
    manifest$target_intersection_mode,
    "full_streamline"
  )
  local_search_radius <- ftq_value(manifest$local_search_radius_mm, 2)
  assignment_label <- ftq_report_assignment_label(
    target_assignment_mode,
    local_search_radius
  )
  manifest_html <- ""
  if (!is.null(manifest)) {
    local_search_radius_detail <- if (identical(
      target_assignment_mode,
      "local_search"
    )) {
      paste0(
        "<br><b>Local search radius (mm):</b> ",
        esc(local_search_radius)
      )
    } else {
      ""
    }
    context_image <- ftq_context_density_image(tractogram_context)
    background_detail <- if (!is.null(manifest$background_image)) {
      paste0("<br><b>Background image:</b> ", esc(manifest$background_image))
    } else {
      "<br><b>Background image:</b> not provided"
    }
    context_detail <- if (!is.null(context_image)) {
      paste0("<br><b>Tractogram-derived density context:</b> ", esc(context_image))
    } else {
      "<br><b>Fallback spatial context:</b> not used"
    }
    manifest_html <- paste0(
      "<details><summary>Methods and run metadata</summary><div class='meta'><p><b>Subject:</b> ", esc(ftq_value(manifest$subject, "Not recorded")),
      "<br><b>Session:</b> ", esc(ftq_value(manifest$session, "")),
	      "<br><b>Description:</b> ", esc(manifest$description),
	      "<br><b>Focal assignment mode:</b> ", esc(ftq_value(manifest$focal_assignment_mode, "all_focal_voxels")),
	      "<br><b>Streamline voxelization mode:</b> ", esc(ftq_value(manifest$streamline_voxelization_mode, "sample_points")),
	      "<br><b>Target assignment mode:</b> ", esc(target_assignment_mode),
	      " (", esc(assignment_label), ")",
	      local_search_radius_detail,
      "<br><b>End-voxel audit:</b> ", esc(ftq_value(manifest$endpoint_analysis, "none")),
      "<br><b>End-voxel contexts:</b> ", esc(if (length(manifest$endpoint_context_labels)) {
        paste(manifest$endpoint_context_labels, collapse = ", ")
      } else {
        "none supplied"
      }),
      "<br><b>Reverse streamlines:</b> ", esc(ftq_value(manifest$reverse_streamlines, FALSE)),
      "<br><b>Discretization:</b> ", esc(paste(manifest$discretization, collapse = ",")),
      "<br><b>Tractogram format:</b> ", esc(ftq_value(manifest$tractogram_format, "Not recorded")),
      "<br><b>Track space:</b> ", esc(manifest$track_space),
      if (identical(manifest$tractogram_format, "tck")) {
        "<br><b>Coordinate mapping:</b> native TCK RAS+ world millimeters mapped through the declared NIfTI affine; no axis flips or overlap-based orientation selection."
      } else "",
      "<br><b>Output directory:</b> ", esc(manifest$output_dir),
      "<br><b>Report axial slices (one-based indices):</b> ",
      esc(if (is.null(manifest$report$slices)) "Automatic selection" else paste(manifest$report$slices, collapse = ", ")),
      background_detail,
      context_detail, "</p></div>", ftq_report_metadata_html(metadata_path, dirname(path)), "</details>"
    )
  }
  local_search_unit_interpretation <- ftq_value(
    validation$local_search_spatial_units_interpretation,
    NA_character_
  )
  local_search_units_known <- !is.na(local_search_unit_interpretation) &&
    !grepl("unknown", local_search_unit_interpretation, ignore.case = TRUE)
  assignment_explanation <- switch(
    target_assignment_mode,
    full_streamline = paste0(
      "A focal-positive streamline is connected to a cNode when any selected ",
      "voxel along its full path intersects that cNode."
    ),
    endpoint_voxels = paste0(
      "A focal-positive streamline is connected only when its exact end voxels ",
      "form an unambiguous opposite fNode-cNode pair."
    ),
    local_search = paste0(
      "A focal-positive streamline is connected only when its two continuous ",
      "end coordinates resolve to an unambiguous opposite fNode-cNode pair after ",
      "nearest labeled voxel-center searches within ",
      format(local_search_radius, trim = TRUE), " mm."
    )
  )
  assignment_warning <- if (identical(target_assignment_mode, "local_search") &&
                            !local_search_units_known) {
    if (is.na(local_search_unit_interpretation)) {
      paste0("<p class='warn'>Spatial-unit validation was not supplied with this report. ",
        "Check the original validation record and reference-image units before interpreting physical distances.</p>")
    } else paste0(
      "<p class='warn'>The reference header does not declare known spatial units. ",
      "This run recorded an assumed millimeter interpretation; publication use ",
      "requires verified scanner-space units or a trusted, non-resampled header correction.</p>"
    )
  } else {
    ""
  }
  assignment_html <- paste0(
    "<h2>Primary Connection Definition</h2><div class='note'><b>",
    esc(assignment_label), ":</b> ", esc(assignment_explanation),
    " Target-density and no-cNode maps use this same definition. The separate ",
    "end-voxel audit below is sensitivity QC and does not replace the primary ",
    "assignment method.</div><p><b>Focal voxel allocation:</b> ",
    esc(ftq_value(manifest$focal_assignment_mode, "Not recorded")), "</p>",
    if (!is.null(manifest$provenance_warning)) paste0("<p class='warn'>", esc(manifest$provenance_warning), "</p>") else "",
    assignment_warning
  )
  build_local_search_html <- function() {
    if (is.null(local_search_data)) return("")
    local_summary <- local_search_data$summary
    local_status <- local_search_data$endpoint_status
    local_display <- local_summary
    if (nrow(local_display)) {
      count_columns <- intersect(c(
        "unique_connections",
        "connections_with_full_path_intersection",
        "connections_without_full_path_intersection"
      ), names(local_display))
      for (column in count_columns) {
        local_display[[column]] <- format_count(local_display[[column]])
      }
      numeric_columns <- intersect(c(
        "full_path_intersection_percent_of_connections",
        "focal_distance_median_mm", "target_distance_median_mm",
        "paired_max_distance_q25_mm", "paired_max_distance_median_mm",
        "paired_max_distance_q75_mm", "paired_max_distance_q95_mm",
        "paired_max_distance_max_mm",
        "full_path_target_distance_median_mm",
        "full_path_max_distance_median_mm",
        "full_path_max_distance_q95_mm",
        "proximity_only_target_distance_median_mm",
        "proximity_only_max_distance_median_mm",
        "proximity_only_max_distance_q95_mm",
        "search_radius_mm"
      ), names(local_display))
      for (column in numeric_columns) {
        local_display[[column]] <- ifelse(
          is.finite(as.numeric(local_summary[[column]])) &
            if (identical(column, "full_path_intersection_percent_of_connections"))
              local_summary$unique_connections > 0 else TRUE,
          sprintf("%.3f", as.numeric(local_summary[[column]])),
          "Unavailable"
        )
      }
      local_labels <- c(
        target = "cNode",
        unique_connections = "Unique local-search connections",
        connections_with_full_path_intersection =
          "Also enters cNode along full path",
        connections_without_full_path_intersection =
          "Proximity-only connection",
        full_path_intersection_percent_of_connections =
          "Also enters cNode (% of connections)",
        focal_distance_median_mm = "Median fNode-end distance (mm)",
        target_distance_median_mm = "Median cNode-end distance (mm)",
        paired_max_distance_q25_mm = "Paired max distance P25 (mm)",
        paired_max_distance_median_mm = "Paired max distance median (mm)",
        paired_max_distance_q75_mm = "Paired max distance P75 (mm)",
        paired_max_distance_q95_mm = "Paired max distance P95 (mm)",
        paired_max_distance_max_mm = "Paired max distance maximum (mm)",
        full_path_target_distance_median_mm =
          "Full-path cNode-end distance median (mm)",
        full_path_max_distance_median_mm =
          "Full-path paired max distance median (mm)",
        full_path_max_distance_q95_mm =
          "Full-path paired max distance P95 (mm)",
        proximity_only_target_distance_median_mm =
          "Proximity-only cNode-end distance median (mm)",
        proximity_only_max_distance_median_mm =
          "Proximity-only paired max distance median (mm)",
        proximity_only_max_distance_q95_mm =
          "Proximity-only paired max distance P95 (mm)",
        search_radius_mm = "Search radius (mm)"
      )
      names(local_display) <- unname(local_labels[names(local_display)])
    }
    status_display <- local_status
    if (nrow(status_display)) {
      status_display$endpoint <- ifelse(
        status_display$endpoint == "start",
        "Stored start",
        "Stored end"
      )
      status_display$status <- gsub(
        "_", " ", status_display$status, fixed = TRUE
      )
      status_display$streamline_count <- format_count(
        status_display$streamline_count
      )
      status_display$percent_of_focal_intersecting_streamlines <- ifelse(
        local_status$denominator_focal_intersecting_streamlines > 0,
        sprintf("%.3f", as.numeric(local_status$percent_of_focal_intersecting_streamlines)),
        "Unavailable")
      status_display$denominator_focal_intersecting_streamlines <-
        format_count(
          status_display$denominator_focal_intersecting_streamlines
        )
      names(status_display) <- c(
        "Stored endpoint", "Nearest-node status", "Streamline count",
        "Percent of focal-intersecting streamlines",
        "Focal-intersecting denominator"
      )
    }
    local_figure <- figure_if_present(
      "local_search_distances",
      "Local-search end-voxel distances and full-path reconciliation",
      paste0(
        "Left: among assigned fNode-cNode pairs, the cumulative distribution ",
        "of the larger of the two endpoint-to-nearest-node distances. Right: ",
        "the same unique primary connections partitioned by whether the full ",
        "voxelized streamline also enters the assigned cNode. Proximity-only ",
        "connections satisfy the stated local-search rule but do not constitute ",
        "full-path mask intersections."
      )
    )
    source_link <- ""
    source <- local_search_data$source_path
    if (is.character(source) && length(source) == 1L && file.exists(source)) {
      href <- ftq_report_href(source, dirname(path))
      source_link <- paste0(
        "<p><a href='", esc(href),
        "'>Download streamline-level local-search audit CSV</a></p>"
      )
    }
    paste0(
      "<h2>Local-Search Connection Audit</h2>",
      "<div class='note'><b>Primary rule and independent reconciliation:</b> ",
      "the primary local-search count requires opposite streamline ends to be ",
      "nearest to the fNode and cNode within the stated physical radius. The ",
      "table below then asks an independent question: did that streamline also ",
      "enter the cNode anywhere along its full voxelized path? These categories ",
      "partition the primary connection count and do not alter it.</div>",
      table_html(local_display[seq_len(min(5L, ncol(local_display)))]),
      "<details><summary>Endpoint distances and full local-search audit</summary>",
      table_html(local_display), local_figure, source_link, "</details>",
      "<details><summary>All endpoint nearest-node statuses</summary>",
      "<p>Each focal-intersecting streamline contributes one stored start and ",
      "one stored end status. Streamlines are treated as unoriented for paired ",
      "connection counting; start and end are retained here only for auditability.</p>",
      table_html(status_display), "</details>"
    )
  }
  background_note <- ""
  context_image <- ftq_context_density_image(tractogram_context)
  if (!is.null(context_image)) {
    reason <- if (identical(ftq_value(manifest$context_source, "auto"), "tractogram")) {
      "because report.context_source requested tractogram context."
    } else {
      "because no anatomical image was used."
    }
    background_note <- paste0(
      "<div class='note'><b>Spatial QC note:</b> Anatomical background images are optional. ",
      "These spatial figures use <b>Tractogram-derived density context</b> ", reason, " ",
      "This context is not anatomical MRI and is not brain, tissue, or lesion volume.</div>"
    )
  } else if (is.null(manifest$background_image)) {
    background_note <- paste0(
      "<div class='note'><b>Spatial QC note:</b> Anatomical background images are optional. ",
      "This report uses map-only spatial figures because no report.background_image or tractogram-derived context was supplied.</div>"
    )
  }
  link_names <- gsub("_", " ", names(top_tables), fixed = TRUE)
  ranked <- vapply(top_tables, function(p) {
    is.character(p) && length(p) == 1L && !is.na(p) &&
      grepl("_ranked_.*_voxels[.]csv$", basename(p))
  }, logical(1))
  link_names[ranked] <- paste0("All voxels, ranked by ", link_names[ranked], " (CSV)")
  link_names[link_names == "endpoint penetration summary"] <-
    "end-voxel penetration summary"
  primary_table_keys <- intersect(c("all_voxels", "all_voxel_columns", "local_fraction", "connection_summary", "fate_summary"), names(top_tables))
  primary_tables <- ftq_report_downloads(top_tables[primary_table_keys], dirname(path),
    c(all_voxels = "All focal voxels (CSV)", all_voxel_columns = "Column roles and coordinates (CSV)",
      local_fraction = "Local fraction and validity (CSV)", connection_summary = "RSC summary (CSV)",
      fate_summary = "Streamline fate summary (CSV)")[primary_table_keys])
  advanced_keys <- setdiff(names(top_tables), primary_table_keys)
  links <- ftq_report_downloads(top_tables[advanced_keys], dirname(path),
                               link_names[match(advanced_keys, names(top_tables))])
  original_links <- ftq_report_downloads(source_paths, dirname(path),
    c(per_voxel = "Original per-voxel results (CSV)", rsc = "Complete RSC fields and compatibility aliases (CSV)",
      roi_summary = "ROI-level unique counts (CSV)", validation = "Validation record (JSON)")[names(source_paths)])
  map_links <- ftq_report_downloads(c(map_paths, list(context_density = context_image)), dirname(path),
    paste0(gsub("[-_]", " ", c(names(map_paths), "context density")), " (NIfTI)"))
  guide_link <- ftq_report_downloads(list(guide), dirname(path), "Output guide: fields, units and interpretation")

  format_count <- function(x) {
    ifelse(is.finite(as.numeric(x)), format(as.numeric(x), big.mark = ",", trim = TRUE), "Unavailable")
  }
  connection_display <- connection_summary
  connection_display$unique_streamline_connections <- format_count(
    connection_display$unique_streamline_connections
  )
  connection_display$connected_fNode_voxels <- format_count(
    connection_display$connected_fNode_voxels
  )
  connection_display$rsc_percent <- ifelse(
    is.finite(connection_summary$rsc_percent),
    sprintf("%.2f", connection_summary$rsc_percent),
    "Unavailable"
  )
  connection_display$total_fNode_voxels <- NULL
  names(connection_display) <- c(
    "Destination",
    "Unique streamline connections",
    "Connected fNode voxels",
    "Focal voxels connected (%) - RSC"
  )
  unique_connections_available <- any(is.finite(
    connection_summary$unique_streamline_connections
  ))

  fate_basis <- unique(as.character(fate_summary$summary_basis))
  unique_fate <- length(fate_basis) == 1L &&
    identical(fate_basis, "unique_focal_intersecting_streamlines")
  fate_display <- fate_summary[c(
    "fate", "count", "percent_of_denominator"
  )]
  fate_display$fate <- gsub("_", " ", fate_display$fate, fixed = TRUE)
  fate_display$count <- format_count(fate_display$count)
  fate_display$percent_of_denominator <- ifelse(fate_summary$denominator_count > 0,
    sprintf("%.2f", as.numeric(fate_summary$percent_of_denominator)), "Unavailable")
  names(fate_display) <- c(
    "Fate", "Count", "Percent of denominator"
  )

  figure_if_present <- function(name, title, caption) {
    if (!name %in% names(figs) || !file.exists(figs[[name]])) return("")
    fig_tag(figs[[name]], title, caption)
  }
  rsc_figure_html <- figure_if_present(
    "rsc",
    "Regional streamline connectivity",
    paste0(
      "RSC is the percentage of fNode voxels with at least one connection in ",
      "the indicated row. A focal voxel can contribute to several rows, so ",
      "these percentages do not need to sum to 100%. Labels report ROI-level ",
      "unique streamline counts; the local row is included for comparison ",
      "and is not a cNode."
    )
  )
  fate_figure_html <- figure_if_present(
    "fate",
    if (unique_fate) "Unique streamline fate" else "Association-weighted fate",
    if (unique_fate) {
      paste0(
        "Each unique focal-intersecting streamline appears once. Local has ",
        "priority; the remaining nonlocal streamlines are separated by whether ",
        "they connect to any cNode under the selected assignment method."
      )
    } else {
      paste0(
        "ROI-level unique-streamline data were unavailable, so this legacy ",
        "display weights a streamline once per occupied fNode voxel."
      )
    }
  )
  tracks_figure_html <- figure_if_present(
    "spatial_tracks",
    "fNode streamline density heatmap",
    paste0(
      "Each focal voxel shows its assigned streamline count under the recorded ",
      "focal assignment rule, regardless of local or cNode status. The color ",
      "scale is capped at the 99th percentile for display only; source CSV and ",
      "NIfTI values are unchanged."
    )
  )
  target_figure_html <- figure_if_present(
    "spatial_targets",
    "fNode connection density heatmaps",
    paste0(
      "The local row and each cNode row show the number of qualifying ",
      "streamlines at each fNode voxel using ", assignment_label, ". Each row ",
      "has a distinct gradient and its own 99th-percentile display scale, so ",
      "compare spatial patterns rather than color intensity between rows."
    )
  )
  no_target_figure_html <- figure_if_present(
    "spatial_leave_no_target",
    "Nonlocal streamlines without an assigned cNode",
    paste0(
      "Density of focal-intersecting streamlines that are not local and do not ",
      "connect to any cNode under ", assignment_label, ". This category changes ",
      "when the primary assignment method or search radius changes."
    )
  )
  local_search_html <- build_local_search_html()

  endpoint_html <- ""
  endpoint_figure_names <- intersect(
    c("endpoint_classification", "endpoint_comparison", "endpoint_penetration"),
    names(figs)
  )
  if (!is.null(endpoint_data)) {
    if (isTRUE(endpoint_data$missing)) {
      endpoint_html <- paste0(
        "<h2>Advanced End-Voxel QC</h2>",
        "<div class='note fail'><b>End-voxel audit unavailable:</b> this run requested the audit, ",
        "but one or more endpoint tables are missing. Existing full-path metrics remain available; ",
        "the endpoint section must not be interpreted as validated.</div>"
      )
    } else {
      summary_columns <- intersect(c(
        "target", "tolerance", "target_intersection_streamlines",
        "paired_focal_target_endpoints", "shared_mask_ambiguous",
        "target_endpoint_only", "focal_endpoint_only", "neither_relevant_endpoint",
        "focal_target_overlap_voxels"
      ), names(endpoint_data$summary))
      context_columns <- intersect(c(
        "context", "context_type", "endpoint_count", "endpoints_inside",
        "endpoints_outside", "endpoints_inside_percent",
        "endpoints_one_voxel_outside", "endpoints_more_than_one_voxel_outside",
        "streamlines_with_both_endpoints_inside_percent"
      ), names(endpoint_data$context))
      summary_display <- endpoint_data$summary[summary_columns]
      summary_display$unambiguous_paired_percent <- rep(NA_real_, nrow(summary_display))
      if (all(c("target_intersection_streamlines", "paired_focal_target_endpoints") %in%
              names(endpoint_data$summary))) {
        summary_display$unambiguous_paired_percent <- ifelse(
          endpoint_data$summary$target_intersection_streamlines > 0,
          100 * endpoint_data$summary$paired_focal_target_endpoints /
            endpoint_data$summary$target_intersection_streamlines, NA_real_)
      }
      summary_display$unambiguous_paired_percent <- ifelse(
        is.finite(summary_display$unambiguous_paired_percent),
        sprintf("%.2f", summary_display$unambiguous_paired_percent), "Unavailable")
      summary_labels <- c(
        target = "cNode",
        tolerance = "End-voxel definition",
        target_intersection_streamlines = "Full-path intersections",
        paired_focal_target_endpoints = "Unambiguous paired fNode-cNode ends",
        unambiguous_paired_percent = "Unambiguous paired (% of full-path intersections)",
        shared_mask_ambiguous = "Shared-mask ambiguous",
        target_endpoint_only = "cNode end only",
        focal_endpoint_only = "fNode end only",
        neither_relevant_endpoint = "Neither relevant end",
        focal_target_overlap_voxels = "fNode-cNode overlap voxels"
      )
      names(summary_display) <- unname(summary_labels[names(summary_display)])
      context_display <- endpoint_data$context[context_columns]
      if ("endpoints_inside_percent" %in% names(context_display)) {
        context_display$endpoints_inside_percent <- ifelse(
          endpoint_data$context$endpoint_count > 0,
          sprintf("%.2f", endpoint_data$context$endpoints_inside_percent), "Unavailable")
      }
      if ("streamlines_with_both_endpoints_inside_percent" %in% names(context_display)) {
        context_display$streamlines_with_both_endpoints_inside_percent <- ifelse(
          endpoint_data$context$valid_streamlines > 0,
          sprintf("%.2f", endpoint_data$context$streamlines_with_both_endpoints_inside_percent), "Unavailable")
      }
      context_labels <- c(
        context = "Context",
        context_type = "Context type",
        endpoint_count = "Total end voxels",
        endpoints_inside = "End voxels inside",
        endpoints_outside = "End voxels outside",
        endpoints_inside_percent = "End voxels inside (%)",
        endpoints_one_voxel_outside = "One voxel outside",
        endpoints_more_than_one_voxel_outside = "More than one voxel outside",
        streamlines_with_both_endpoints_inside_percent =
          "Streamlines with both ends inside (%)"
      )
      names(context_display) <- unname(context_labels[names(context_display)])
      endpoint_figure_specs <- list(
        endpoint_classification = c(
          "Exact end-voxel classes",
          paste0(
            "Every full-path target intersection is assigned to exactly one ",
            "end-voxel class. The denominator is the full-path intersection ",
            "count for that cNode, not the number of all focal streamlines."
          )
        ),
        endpoint_comparison = c(
          "Exact versus one-voxel end-voxel sensitivity",
          paste0(
            "Compares exact fNode-cNode end-voxel pairs with a 26-neighbor ",
            "one-voxel dilation. This sensitivity is voxel-grid based and is ",
            "not the same as the primary physical-distance local search. ",
            "These curves use geometric pairs, including shared-mask membership, ",
            "whereas the classification table excludes shared-mask ambiguity ",
            "from its paired class. Each panel states its own denominator."
          )
        ),
        endpoint_penetration = c(
          "Target penetration of full-path intersections",
          paste0(
            "Intersection-positive streamlines are placed into mutually ",
            "exclusive bins based on the number of unique exact cNode voxels ",
            "visited along the selected full path."
          )
        )
      )
      endpoint_figures <- paste(vapply(endpoint_figure_names, function(nm) {
        spec <- endpoint_figure_specs[[nm]]
        fig_tag(figs[[nm]], spec[[1L]], spec[[2L]])
      }, character(1)), collapse = "\n")
      source_links <- character()
      if (length(endpoint_data$source_paths)) {
        for (nm in names(endpoint_data$source_paths)) {
          source <- endpoint_data$source_paths[[nm]]
          if (is.character(source) && length(source) == 1L && file.exists(source)) {
            href <- ftq_report_href(source, dirname(path))
            source_links <- c(
              source_links,
              paste0("<li><a href='", esc(href), "'>", esc(nm), " CSV</a></li>")
            )
          }
        }
      }
      source_link_html <- if (length(source_links)) {
        paste0("<h3>Machine-Readable End-Voxel Tables</h3><ul>", paste(source_links, collapse = ""), "</ul>")
      } else {
        ""
      }
      penetration_note <- if (any(
        endpoint_data$penetration$streamline_voxelization_mode == "sample_points"
      )) {
        paste0(
          "<div class='note'><b>Penetration interpretation:</b> these values count unique cNode voxels ",
          "represented by stored streamline points. Sparse stored points can undercount continuously traversed voxels.</div>"
        )
      } else {
        ""
      }
      penetration_columns <- intersect(c(
        "target", "target_intersection_streamlines",
        "one_voxel_touch_percent_of_target_intersections",
        "two_voxel_only_percent_of_target_intersections",
        "three_to_four_voxels_percent_of_target_intersections",
        "five_to_nine_voxels_percent_of_target_intersections",
        "ten_or_more_voxels_percent_of_target_intersections",
        "penetration_p25_unique_target_voxels",
        "penetration_median_unique_target_voxels",
        "penetration_p75_unique_target_voxels",
        "penetration_p95_unique_target_voxels"
      ), names(endpoint_data$penetration_summary))
      penetration_display <- endpoint_data$penetration_summary[
        penetration_columns
      ]
      percent_columns <- grep("percent", names(penetration_display), value = TRUE)
      for (column in percent_columns) {
        penetration_display[[column]] <- ifelse(
          penetration_display$target_intersection_streamlines > 0,
          sprintf("%.2f", penetration_display[[column]]), "Unavailable")
      }
      penetration_labels <- c(
        target = "cNode",
        target_intersection_streamlines = "Full-path intersections",
        one_voxel_touch_percent_of_target_intersections = "1 voxel (%)",
        two_voxel_only_percent_of_target_intersections = "2 voxels (%)",
        three_to_four_voxels_percent_of_target_intersections = "3-4 voxels (%)",
        five_to_nine_voxels_percent_of_target_intersections = "5-9 voxels (%)",
        ten_or_more_voxels_percent_of_target_intersections = "10+ voxels (%)",
        penetration_p25_unique_target_voxels = "25th percentile voxels",
        penetration_median_unique_target_voxels = "Median voxels",
        penetration_p75_unique_target_voxels = "75th percentile voxels",
        penetration_p95_unique_target_voxels = "95th percentile voxels"
      )
      names(penetration_display) <- unname(
        penetration_labels[names(penetration_display)]
      )
      endpoint_html <- paste0(
        "<h2>Advanced End-Voxel QC</h2><details><summary>",
        "Full-Path Intersections vs End-Voxel Connectivity</summary>",
        "<div class='note'><b>Different estimands:</b> a full-path intersection ",
        "means any selected streamline voxel entered a cNode. Paired end voxels ",
        "ask whether opposite ends lie in the fNode and cNode. The five classes ",
        "are mutually exclusive after shared-mask ambiguity receives priority. ",
        "The one-voxel result is a 26-neighbor grid sensitivity analysis, not a ",
        "uniform physical-distance tolerance and not the primary local-search method.</div>",
        "<h3>End-Voxel Classification</h3>",
        "<p><b>Paired:</b> opposite unambiguous fNode and cNode ends. ",
        "<b>Shared-mask ambiguity:</b> an end belongs to both masks. ",
        "<b>cNode end only, fNode end only, and neither:</b> only the named relevant ",
        "end membership is present.</p>",
        table_html(summary_display),
        endpoint_figures,
        "<h3>cNode Penetration</h3>", penetration_note,
        "<p>Penetration counts distinct exact cNode voxels visited per full-path ",
        "target-positive streamline; it is not distance in millimetres or tissue volume. ",
        "The displayed percentages use mutually exclusive bins (1, 2, 3-4, ",
        "5-9, or 10+ unique cNode voxels), so each nonempty row sums to 100% before rounding. The ",
        "quartiles summarize the same per-streamline voxel-count distribution.</p>",
        table_html(penetration_display),
        "<h3>Endpoint coverage</h3>",
        "<p>Reference-grid rows test whether end coordinates fall inside the ",
        "image grid. Optional context rows test containment in a supplied, ",
        "trusted anatomical mask; they are unavailable when no such mask was supplied. ",
        "These rows describe all processed input endpoints, not just the focal-intersecting ",
        "population. Even 100% image-grid coverage does not imply full-path containment ",
        "in localROI, anatomical coverage, or a valid focal-target connection.</p>",
        table_html(context_display),
        source_link_html, "</details>"
      )
    }
  }
  distribution_html <- figure_if_present(
    "tracks",
    "Per-fNode-voxel streamline-count distribution",
    paste0(
      "The histogram shows the distribution of log-transformed per-fNode-voxel ",
      "streamline counts. The ECDF shows the fraction of occupied fNode voxels ",
      "at or below each raw count. This is a density diagnostic, not a ",
      "streamline-fate or connection summary."
    )
  )
  if (nzchar(distribution_html)) {
    distribution_html <- paste0(
      "<details><summary>Advanced density-distribution diagnostic</summary>",
      distribution_html, "</details>"
    )
  }
  connection_note <- paste0(
    "<div class='note'><b>How to read this table:</b> unique streamline ",
    "connections count each qualifying streamline once per row. Connected ",
    "fNode voxels count fNode voxels with at least one such connection. RSC is ",
    "100 * connected fNode voxels / all focal voxels, including voxels with ",
    "no assigned streamlines. A focal voxel can contain streamlines connected ",
    "to several targets, even when each streamline is assigned to only one ",
    "target. The local row uses local streamline status and can also overlap ",
    "target-positive rows. These are overlapping coverage percentages, not ",
    "parts of a whole: rows should not be summed and do not need to add to ",
    "100%. For a partition of unique streamlines, use the Fate Summary below.</div>",
    if (!unique_connections_available) {
      paste0(
        "<p class='warn'>Unique connection counts are unavailable because this ",
        "legacy report call did not supply an ROI-summary table.</p>"
      )
    } else {
      ""
    }
  )
  fate_note <- if (unique_fate) {
    paste0(
      "<div class='note'><b>Mutually exclusive unique-streamline fates:</b> ",
      "local means the selected voxelized path remains inside localROI and has ",
      "priority in this display, even if that streamline also reaches a cNode. ",
      "Among the remaining nonlocal streamlines, leave no target means no cNode ",
      "was assigned and leave with target means at least one cNode was assigned. ",
      "The three counts sum exactly to all unique focal-intersecting streamlines.</div>"
    )
  } else {
    paste0(
      "<div class='note warn'><b>Legacy association-weighted fate:</b> an ",
      "ROI-summary table was not supplied, so a streamline can contribute once ",
      "for every occupied fNode voxel. Rerun or call ftq_report() with ",
      "roi_summary for unique-streamline publication counts.</div>"
    )
  }
  voxel_count <- nrow(result)
  zero_count <- sum(result$tracks == 0)
  coverage_html <- paste0(
    "<p><b>Focal voxel coverage:</b> ", format_count(voxel_count - zero_count), " of ",
    format_count(voxel_count), " focal voxels have assigned streamlines; ", format_count(zero_count),
    " have zero assignments. This is a voxel count, separate from the three streamline fates. ",
    if (identical(manifest$focal_assignment_mode, "first_intersection")) {
      "Under first_intersection, zero assignments do not establish that no streamline traversed the voxel."
    } else "", "</p>"
  )
  endpoint_warning <- ""
  if (!is.null(endpoint_data) && !isTRUE(endpoint_data$missing)) {
    grid <- endpoint_data$context[endpoint_data$context$context_type == "reference_grid", , drop = FALSE]
    if (nrow(grid) && any(grid$endpoints_outside > 0 | grid$endpoints_nonfinite > 0, na.rm = TRUE)) {
      endpoint_warning <- paste0("<p class='warn'>Endpoint coverage needs review: ",
        format_count(grid$endpoints_outside[1L]), " endpoints are outside the reference grid and ",
        format_count(grid$endpoints_nonfinite[1L]), " are nonfinite. See Endpoint coverage in advanced QC.</p>")
    }
  }
  fraction_figure_html <- figure_if_present("local_fraction", "Local fraction",
    paste0("Local assigned streamlines divided by all assigned streamlines at each focal voxel (0 to 1). ",
      "Colored zero values are defined; voxels without assigned streamlines are uncolored and undefined. ",
      "The CSV uses NA for undefined values. Read the NIfTI together with its validity mask. ",
      "Local refers to the supplied localROI, which may be a regional union."))
  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width,initial-scale=1'>",
    "<title>", esc(title), "</title>", css, "</head><body>",
    "<h1>", esc(title), "</h1>",
    "<!-- FIBRTRAQ_NAVIGATION --><!-- /FIBRTRAQ_NAVIGATION -->",
    guide_link,
    assignment_html,
    endpoint_warning,
    preflight_html, validation_html,
    "<h2>RSC summary</h2><p>Total number of fNode voxels: ", format_count(nrow(result)),
    ".</p>", connection_note,
    table_html(connection_display), rsc_figure_html,
    "<h2>Spatial Density Heatmaps</h2>", background_note,
    tracks_figure_html, target_figure_html, no_target_figure_html, fraction_figure_html,
    "<details><summary>Download spatial maps (NIfTI)</summary>", map_links, "</details>",
    "<h2>Fate Summary</h2><p>Denominator: ",
    paste(unique(format_count(fate_summary$denominator_count)), collapse = ", "),
    if (unique_fate) " unique focal-intersecting streamlines." else " focal voxel-streamline associations.",
    "</p>", fate_note,
    table_html(fate_display), coverage_html, fate_figure_html,
    local_search_html,
    endpoint_html,
    distribution_html,
    "<h2>Download Tables</h2><p>The complete coordinate-aware view includes every focal voxel. ",
    "The column-role file preserves the meaning of coordinate and target fields.</p>", primary_tables,
    "<details><summary>Original results and advanced tables</summary>", original_links,
    "<p>Ranked CSVs include every focal voxel, including zero counts, sorted by the stated metric.</p>", links, "</details>",
    manifest_html,
    "</body></html>"
  )
  writeLines(html, path)
}
