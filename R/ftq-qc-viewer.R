#' Build an interactive 3D streamline QC viewer
#'
#' Sample real streamlines from a FibrTraQ manifest and write a standalone HTML
#' viewer for spatial quality control. The viewer renders a boundary
#' representation of the focal ROI, a target-balanced quota subset of positive
#' streamlines, target ROI extents, end voxels, and context points from an
#' anatomical image or tractogram-derived density context. Positive streamlines
#' are accepted in tractogram scan order until each target quota or the scan
#' limit is reached; the subset is not random or prevalence weighted. Sampling
#' follows the manifest's primary `target_intersection_mode`. Endpoint-enabled
#' manifests also add exact versus one-voxel sensitivity controls,
#' endpoint-class filters, and optional focal/target hit markers without
#' changing the selected-method target quotas. The initial display uses an
#' anatomical RAS basis with superior up and an orientation compass. Endpoint
#' markers use continuous streamline coordinates, with colours determined by
#' voxel-based classifications. Circle/square shapes identify the first/last
#' stored point, not biological direction. Sampling and method details are
#' available below the canvas; the nonrandom-subset warning remains visible.
#' Drag to rotate freely, drag near an edge to roll, and scroll or pinch with
#' two fingers to zoom. Reset Orientation restores the initial view and zoom.
#'
#' @param manifest Manifest path/list, or a completed object returned by
#'   [ftq_run()].
#' @param output_html Optional output HTML path. If omitted, the viewer is
#'   written under the run's report figures directory.
#' @param target_labels Optional subset of manifest target labels to sample.
#' @param sample_per_target Number of streamlines to sample per target. Supply
#'   a single value for all targets, or a named integer vector with per-target
#'   quotas.
#' @param max_streamlines Maximum streamlines to scan while sampling.
#' @param max_points_per_streamline Maximum points retained per rendered
#'   streamline segment.
#' @param brain_points Maximum anatomical context points sampled from
#'   `background_image` or tractogram-derived density context.
#' @param target_points Maximum target ROI boundary points rendered per target.
#' @param focal_window Number of points on each side of a focal hit to include
#'   in the focal segment mode.
#' @param tractogram_context Optional tractogram-derived context output or
#'   density image path used when no anatomical background is available.
#' @param overwrite Whether to replace an existing viewer file.
#' @param open Whether to open the viewer after writing.
#'
#' @return A list with viewer path and sampling summary.
#' @export
ftq_interactive_qc_viewer <- function(manifest,
                                      output_html = NULL,
                                      target_labels = NULL,
                                      sample_per_target = 25L,
                                      max_streamlines = 3000000L,
                                      max_points_per_streamline = 140L,
                                      brain_points = 12000L,
                                      target_points = 2000L,
                                      focal_window = 3L,
                                      overwrite = FALSE,
                                      open = FALSE,
                                      tractogram_context = NULL) {
  run_result <- NULL
  run_context <- NULL
  is_completed_run <- is.list(manifest) &&
    "manifest" %in% names(manifest) &&
    "outputs" %in% names(manifest)
  m <- if (is_completed_run) {
    run_result <- manifest[["result"]]
    run_context <- manifest[["outputs"]][["context"]]
    manifest[["manifest"]]
  } else {
    ftq_read_manifest(manifest)
  }

  max_streamlines <- as.integer(max_streamlines)[1]
  max_points_per_streamline <- as.integer(max_points_per_streamline)[1]
  brain_points <- as.integer(brain_points)[1]
  target_points <- as.integer(target_points)[1]
  focal_window <- as.integer(focal_window)[1]
  if (is.na(max_streamlines) || max_streamlines <= 0L) stop("max_streamlines must be positive.")
  if (is.na(max_points_per_streamline) || max_points_per_streamline <= 1L) {
    stop("max_points_per_streamline must be greater than 1.")
  }
  if (is.na(target_points) || target_points < 0L) stop("target_points must be non-negative.")
  if (is.na(focal_window) || focal_window < 0L) stop("focal_window must be non-negative.")

  if (is.null(target_labels)) {
    target_labels <- m$target_labels
  } else {
    target_labels <- as.character(target_labels)
    missing <- setdiff(target_labels, m$target_labels)
    if (length(missing)) {
      stop("target_labels are not present in the manifest: ", paste(missing, collapse = ", "))
    }
  }
  sample_quota <- ftq_qc_sample_quota(sample_per_target, target_labels)
  if (!length(target_labels) || !any(sample_quota > 0L)) {
    stop("At least one target must have a positive sample quota.")
  }

  if (is.null(output_html)) {
    output_html <- file.path(
      m$output_dir,
      "reports",
      m$prefix,
      "report-figures",
      paste0(m$prefix, "_fig_interactive-3d-qc.html")
    )
  }
  output_html <- normalizePath(path.expand(output_html), mustWork = FALSE)
  if (file.exists(output_html) && !isTRUE(overwrite)) {
    stop("Interactive QC viewer already exists. Use overwrite = TRUE to replace: ", output_html)
  }
  dir.create(dirname(output_html), recursive = TRUE, showWarnings = FALSE)

  target_intersection_mode <- ftq_normalize_target_intersection_mode(
    m$target_intersection_mode
  )
  local_search_radius_mm <- ftq_validate_local_search_radius(
    ftq_value(m[["local_search_radius_mm"]], 2)
  )
  assignment_target_labels <- if (identical(
    target_intersection_mode,
    "full_streamline"
  )) {
    target_labels
  } else {
    m$target_labels
  }

  focal <- ftq_qc_load_roi(m$fnode)
  assignment_targets <- lapply(
    match(assignment_target_labels, m$target_labels),
    function(i) ftq_qc_load_roi(m$target_paths[[i]])
  )
  names(assignment_targets) <- assignment_target_labels
  targets <- assignment_targets[target_labels]
  local_context <- NULL
  if (!is.null(m[["local"]]) && file.exists(m[["local"]])) {
    local_context <- ftq_qc_load_roi(m[["local"]])
  }
  vol_dim <- dim(focal$mask)[1:3]
  ref_image <- if (!is.null(m$ref_image)) RNifti::readNifti(m$ref_image, internal = FALSE) else focal$image

  focal_lut <- ftq_qc_make_lut(focal$coords, vol_dim)
  assignment_target_luts <- lapply(
    assignment_targets,
    function(x) ftq_qc_make_lut(x$coords, vol_dim)
  )
  target_luts <- assignment_target_luts[target_labels]
  endpoint_enabled <- identical(
    ftq_value(m$endpoint_analysis, "exact_and_dilated1"),
    "exact_and_dilated1"
  )
  focal_lut_dilated <- if (endpoint_enabled) {
    ftq_endpoint_dilate_lut(focal_lut, vol_dim)
  } else {
    rep(FALSE, length(focal_lut))
  }
  target_luts_dilated <- if (endpoint_enabled) {
    lapply(target_luts, ftq_endpoint_dilate_lut, vol_dim = vol_dim)
  } else {
    lapply(target_luts, function(x) rep(FALSE, length(x)))
  }
  endpoint_context_luts <- list()
  if (endpoint_enabled && length(m$endpoint_context_paths)) {
    endpoint_context_luts <- lapply(m$endpoint_context_paths, function(path) {
      mask <- ftq_qc_load_roi(path)
      ftq_qc_make_lut(mask$coords, vol_dim)
    })
    names(endpoint_context_luts) <- m$endpoint_context_labels
  }
  target_lut_matrix <- if (length(assignment_target_luts)) {
    do.call(cbind, assignment_target_luts)
  } else {
    matrix(FALSE, nrow = prod(vol_dim), ncol = 0L)
  }
  methods <- ftq_qc_discretization_methods(m$discretization)
  local_search_geometry <- if (identical(
    target_intersection_mode,
    "local_search"
  )) {
    ftq_local_search_geometry(focal$image, ref_image = ref_image)
  } else {
    list(
      voxel_to_world = diag(3L),
      voxel_radius_per_mm = rep(1, 3L)
    )
  }
  local_search_voxel_radius <- if (identical(
    target_intersection_mode,
    "local_search"
  )) {
    ftq_local_search_voxel_radius(
      local_search_radius_mm,
      local_search_geometry
    )
  } else {
    rep(0, 3L)
  }

  tractogram <- ftq_read_tractogram(
    m$tractogram, reader = m$track_reader, chunk_size = m$track_chunk_size
  )
  tractogram_format <- ftq_track_format(tractogram)
  stream <- if (ftq_is_track_stream(tractogram)) tractogram else NULL
  track_list <- if (is.null(stream)) tractogram$tracks else NULL
  track_affine <- NULL
  track_voxel_size <- NULL
  selected_track_space <- ftq_value(attr(run_result, "fibrtraq_track_space"), m$track_space)
  axis_flip <- ftq_value(attr(run_result, "fibrtraq_axis_flip"), m$axis_flip, c(FALSE, FALSE, FALSE))
  if (identical(tractogram_format, "tck")) {
    settings <- ftq_tck_coordinate_settings(selected_track_space, axis_flip, ref_image)
    spatial_images <- c(
      list(focal$image), lapply(assignment_targets, `[[`, "image"),
      if (!is.null(local_context)) list(local_context$image),
      lapply(m$endpoint_context_paths, RNifti::readNifti, internal = FALSE)
    )
    for (spatial_image in spatial_images) {
      ftq_tck_validate_reference(spatial_image)
      ftq_validate_direct_reference_grid(focal$image, spatial_image, vol_dim)
    }
    ftq_validate_direct_reference_grid(focal$image, settings$reference, vol_dim)
    track_space <- settings$track_space
    axis_flip <- settings$axis_flip
    ref_image <- settings$reference
  } else {
    track_affine <- ftq_trackvis_affine(tractogram$header)
    track_voxel_size <- ftq_trackvis_voxel_size(tractogram$header)
    track_space <- ftq_qc_resolve_track_space(
      track_space = selected_track_space,
      track_affine = track_affine,
      ref_image = ref_image,
      track_voxel_size = track_voxel_size
    )
  }

  samples <- stats::setNames(vector("list", length(target_labels)), target_labels)
  for (target in target_labels) samples[[target]] <- list()
  hits_seen <- stats::setNames(integer(length(target_labels)), target_labels)
  scanned <- 0L
  focal_hits <- 0L
  done <- FALSE

  map_track <- function(track) {
    mapped <- ftq_qc_map_track(
      track = track,
      track_space = track_space,
      ref_image = ref_image,
      track_affine = track_affine,
      track_voxel_size = track_voxel_size,
      streamline_voxelization_mode = m$streamline_voxelization_mode
    )
    mapped$render_world <- ftq_qc_render_world(
      voxel = mapped$voxel,
      ref_image = ref_image,
      axis_flip = axis_flip,
      vol_dim = vol_dim
    )
    mapped
  }

  target_assignment <- function(mapped, focal_hit) {
    chosen_keys <- focal_hit$keys
    valid_keys <- unique(chosen_keys[!is.na(chosen_keys)])
    full_path_target_hit <- if (length(valid_keys) && ncol(target_lut_matrix)) {
      colSums(target_lut_matrix[valid_keys, , drop = FALSE]) > 0L
    } else {
      rep(FALSE, ncol(target_lut_matrix))
    }
    assignment <- ftq_primary_target_hits_r(
      mode = target_intersection_mode,
      full_path_target_hit = full_path_target_hit,
      chosen_keys = chosen_keys,
      mapped = mapped$voxel,
      two_endpoints = nrow(mapped$voxel) > 1L,
      lut_f = focal_lut,
      target_lut = target_lut_matrix,
      axis_flip = axis_flip,
      vol_dim = vol_dim,
      local_search_radius_mm = local_search_radius_mm,
      local_search_geometry = local_search_geometry
    )
    details <- stats::setNames(
      vector("list", length(assignment_target_labels)),
      assignment_target_labels
    )
    for (target_index in which(assignment$target_hit)) {
      if (identical(target_intersection_mode, "full_streamline")) {
        indices <- which(
          !is.na(chosen_keys) &
            target_lut_matrix[chosen_keys, target_index]
        )
        marker_voxels <- focal_hit$discrete[indices, , drop = FALSE]
      } else if (identical(
        target_intersection_mode,
        "endpoint_voxels"
      )) {
        endpoint_indices <- c(1L, length(chosen_keys))
        endpoint_labels <- vapply(
          chosen_keys[endpoint_indices],
          ftq_unique_node_label_at_key,
          integer(1),
          lut_f = focal_lut,
          target_lut = target_lut_matrix
        )
        selected <- which(endpoint_labels == target_index + 1L)
        indices <- endpoint_indices[selected]
        marker_voxels <- focal_hit$discrete[indices, , drop = FALSE]
      } else {
        endpoint_indices <- c(1L, nrow(mapped$voxel))
        endpoint_labels <- vapply(
          assignment$endpoints,
          function(x) x$label,
          integer(1)
        )
        selected <- which(endpoint_labels == target_index + 1L)
        indices <- endpoint_indices[selected]
        marker_voxels <- do.call(
          rbind,
          lapply(assignment$endpoints[selected], function(x) x$voxel)
        )
      }
      marker_world <- ftq_voxel_to_world_mm(marker_voxels, ref_image)
      details[[target_index]] <- list(
        hit = TRUE,
        indices = as.integer(indices),
        world = marker_world
      )
    }
    list(target_hit = assignment$target_hit, details = details)
  }

  record_sample <- function(track, chosen) {
    mapped <- tryCatch(map_track(track), error = function(e) NULL)
    if (is.null(mapped)) return(FALSE)

    focal_hit <- ftq_qc_hit_detail(
      mapped$voxel,
      focal_lut,
      vol_dim = vol_dim,
      methods = methods,
      axis_flip = axis_flip
    )
    if (!isTRUE(focal_hit$hit)) return(FALSE)

    assignment <- target_assignment(mapped, focal_hit)
    target_detail <- assignment$details[[chosen]]
    if (is.null(target_detail) || !isTRUE(target_detail$hit)) return(FALSE)
    endpoint_annotation <- if (endpoint_enabled) {
      ftq_qc_endpoint_annotation(
        mapped = mapped,
        focal_method = focal_hit$method,
        focal_lut = focal_lut,
        focal_lut_dilated = focal_lut_dilated,
        target_lut = target_luts[[chosen]],
        target_lut_dilated = target_luts_dilated[[chosen]],
        context_luts = endpoint_context_luts,
        axis_flip = axis_flip,
        vol_dim = vol_dim
      )
    } else {
      NULL
    }

    samples[[chosen]][[length(samples[[chosen]]) + 1L]] <<- ftq_qc_sample_record(
      target = chosen,
      world = mapped$render_world,
      focal_indices = focal_hit$indices,
      target_indices = target_detail$indices,
      focal_hit_world = ftq_qc_hit_world(focal_hit, ref_image),
      target_hit_world = target_detail$world,
      endpoint_annotation = endpoint_annotation,
      max_points = max_points_per_streamline,
      focal_window = focal_window
    )
    TRUE
  }

  process_track <- function(track) {
    scanned <<- scanned + 1L
    if (scanned > max_streamlines) {
      done <<- TRUE
      return(invisible(NULL))
    }
    mapped <- tryCatch(map_track(track), error = function(e) NULL)
    if (is.null(mapped)) return(invisible(NULL))

    focal_hit <- ftq_qc_hit_detail(
      mapped$voxel,
      focal_lut,
      vol_dim = vol_dim,
      methods = methods,
      axis_flip = axis_flip
    )
    if (!isTRUE(focal_hit$hit)) return(invisible(NULL))

    focal_hits <<- focal_hits + 1L
    assignment <- target_assignment(mapped, focal_hit)
    selected_target_hit <- assignment$target_hit[
      match(target_labels, assignment_target_labels)
    ]
    available <- target_labels[
      selected_target_hit &
        lengths(samples) < sample_quota[target_labels]
    ]
    if (!length(available)) return(invisible(NULL))

    chosen <- available[which.min(lengths(samples)[available])]
    hits_seen[chosen] <<- hits_seen[chosen] + 1L
    endpoint_annotation <- if (endpoint_enabled) {
      ftq_qc_endpoint_annotation(
        mapped = mapped,
        focal_method = focal_hit$method,
        focal_lut = focal_lut,
        focal_lut_dilated = focal_lut_dilated,
        target_lut = target_luts[[chosen]],
        target_lut_dilated = target_luts_dilated[[chosen]],
        context_luts = endpoint_context_luts,
        axis_flip = axis_flip,
        vol_dim = vol_dim
      )
    } else {
      NULL
    }
    samples[[chosen]][[length(samples[[chosen]]) + 1L]] <<- ftq_qc_sample_record(
      target = chosen,
      world = mapped$render_world,
      focal_indices = focal_hit$indices,
      target_indices = assignment$details[[chosen]]$indices,
      focal_hit_world = ftq_qc_hit_world(focal_hit, ref_image),
      target_hit_world = assignment$details[[chosen]]$world,
      endpoint_annotation = endpoint_annotation,
      max_points = max_points_per_streamline,
      focal_window = focal_window
    )
    if (all(lengths(samples) >= sample_quota)) done <<- TRUE
    invisible(NULL)
  }

  use_cpp_sampler <- !is.null(stream) && ftq_qc_compiled_sampler_available()
  if (use_cpp_sampler) {
    raw_to_voxel <- ftq_qc_raw_to_voxel_affine(
      track_space = track_space,
      ref_image = ref_image,
      track_affine = track_affine,
      track_voxel_size = track_voxel_size
    )
    cpp_quota <- stats::setNames(
      integer(length(assignment_target_labels)),
      assignment_target_labels
    )
    cpp_quota[target_labels] <- sample_quota[target_labels]
    cpp_counts <- stats::setNames(
      integer(length(assignment_target_labels)),
      assignment_target_labels
    )
    ftq_track_stream_iterate(
      stream,
      callback = function(chunk, start_index, end_index) {
        if (done) return(FALSE)
        res <- ftq_qc_sample_chunk_cpp(
          tracks = chunk,
          raw_to_voxel = raw_to_voxel,
          methods = methods,
          axis_flip = as.logical(axis_flip),
          vol_dim = as.integer(vol_dim),
          lut_f = as.logical(focal_lut),
          target_lut = target_lut_matrix,
          sample_quota = as.integer(cpp_quota),
          current_counts = as.integer(cpp_counts),
          streamline_voxelization_mode = m$streamline_voxelization_mode,
          target_intersection_mode = target_intersection_mode,
          local_search_radius_mm = local_search_radius_mm,
          local_search_voxel_to_world =
            local_search_geometry$voxel_to_world,
          local_search_voxel_radius = local_search_voxel_radius
        )
        scanned <<- scanned + as.integer(res$processed)
        focal_hits <<- focal_hits + as.integer(res$focal_hits)
        seen <- as.integer(res$target_seen)
        names(seen) <- assignment_target_labels
        hits_seen <<- hits_seen + seen[target_labels]

        track_indices <- as.integer(res$track_indices)
        target_indices <- as.integer(res$target_indices)
        if (length(track_indices)) {
          for (i in seq_along(track_indices)) {
            target <- assignment_target_labels[[target_indices[[i]]]]
            record_sample(chunk[[track_indices[[i]]]], target)
          }
        }
        cpp_counts[] <<- 0L
        cpp_counts[target_labels] <<- lengths(samples)
        done <<- all(lengths(samples) >= sample_quota) || scanned >= max_streamlines
        if (done) return(FALSE)
        TRUE
      },
      max_tracks = max_streamlines,
      chunk_size = stream$chunk_size
    )
  } else if (!is.null(stream)) {
    ftq_track_stream_iterate(
      stream,
      callback = function(chunk, start_index, end_index) {
        for (track in chunk) {
          if (done) break
          process_track(track)
        }
        if (done) return(FALSE)
        TRUE
      },
      max_tracks = max_streamlines,
      chunk_size = stream$chunk_size
    )
  } else {
    for (track in track_list) {
      if (done) break
      process_track(track)
    }
  }

  sample_counts <- lengths(samples)
  if (!any(sample_counts)) {
    stop(
      "No sampled focal-positive streamlines were connected to the selected ",
      "targets under target_intersection_mode = '",
      target_intersection_mode,
      "'."
    )
  }

  focal_render <- ftq_qc_focal_boundary_coords(focal$mask, focal$coords, max_points = 8000L)
  focal_world <- ftq_voxel_to_world_mm(focal_render$coords, focal$image)
  target_centers <- ftq_qc_target_centers(targets)
  target_boundaries <- ftq_qc_target_boundaries(targets, max_points = target_points)
  context_image <- ftq_context_density_image(ftq_value(tractogram_context, run_context, m$tractogram_context, NULL))
  background <- ftq_qc_context_cloud(
    background_image = m$background_image,
    tractogram_context_image = context_image,
    context_mask = if (!is.null(local_context)) local_context$mask else NULL,
    max_points = brain_points
  )
  colors <- stats::setNames(grDevices::hcl.colors(length(target_labels), "Dark 3"), target_labels)

  ftq_qc_write_html(
    path = output_html,
    streamlines = unlist(samples, recursive = FALSE),
    focal_world = focal_world,
    focal_total_count = nrow(focal$coords),
    focal_rendered_count = nrow(focal_render$coords),
    focal_render_mode = focal_render$mode,
    target_centers = target_centers,
    target_boundaries = target_boundaries,
    background_world = background$points,
    context_label = background$label,
    colors = colors,
    scanned = scanned,
    focal_hits = focal_hits,
    sample_counts = sample_counts[target_labels],
    endpoint_analysis = endpoint_enabled,
    target_intersection_mode = target_intersection_mode,
    local_search_radius_mm = local_search_radius_mm,
    title = paste0("Interactive 3D QC: ", m$prefix)
  )

  if (isTRUE(open)) ftq_qc_open_output(output_html)
  summary <- data.frame(
    target = target_labels,
    sampled_streamlines = as.integer(sample_counts[target_labels]),
    target_hits_seen_while_sampling = as.integer(hits_seen[target_labels]),
    stringsAsFactors = FALSE
  )
  summary <- rbind(
    data.frame(
      target = c("streamlines_scanned", "focal_hits_seen"),
      sampled_streamlines = c(scanned, focal_hits),
      target_hits_seen_while_sampling = NA_integer_,
      stringsAsFactors = FALSE
    ),
    summary
  )
  list(
    html = output_html,
    summary = summary,
    target_labels = target_labels,
    track_space = track_space,
    target_intersection_mode = target_intersection_mode,
    local_search_radius_mm = if (identical(
      target_intersection_mode,
      "local_search"
    )) {
      local_search_radius_mm
    } else {
      NA_real_
    },
    focal_voxels = nrow(focal$coords)
  )
}

ftq_qc_load_roi <- function(path) {
  img <- RNifti::readNifti(path, internal = FALSE)
  list(image = img, mask = as.array(img) > 0, coords = iso_coords(img))
}

ftq_qc_sample_quota <- function(sample_per_target, target_labels) {
  quota <- as.integer(sample_per_target)
  if (!length(quota) || any(is.na(quota))) {
    stop("sample_per_target must contain finite integer quotas.")
  }
  if (length(quota) == 1L) {
    quota <- rep(quota, length(target_labels))
    names(quota) <- target_labels
  } else {
    quota_names <- names(quota)
    if (is.null(quota_names) || any(!nzchar(quota_names))) {
      if (length(quota) != length(target_labels)) {
        stop("Unnamed vector sample_per_target values must match target_labels length.")
      }
      names(quota) <- target_labels
    } else {
      missing <- setdiff(target_labels, quota_names)
      if (length(missing)) {
        stop("sample_per_target is missing quota(s) for target label(s): ",
             paste(missing, collapse = ", "))
      }
      quota <- quota[target_labels]
    }
  }
  if (any(quota < 0L)) stop("sample_per_target quotas must be non-negative.")
  quota
}

ftq_qc_discretization_methods <- function(discretization) {
  if (identical(discretization, "any") || any(discretization == "any")) {
    c("round", "floor", "ceil")
  } else {
    unique(match.arg(discretization, c("round", "floor", "ceil"), several.ok = TRUE))
  }
}

ftq_qc_make_lut <- function(coords, vol_dim) {
  lut <- rep(FALSE, prod(vol_dim))
  key <- ftq_qc_encode_voxels(coords, vol_dim)
  lut[key[!is.na(key)]] <- TRUE
  lut
}

ftq_qc_encode_voxels <- function(coords, vol_dim) {
  mm <- as.matrix(coords)[, 1:3, drop = FALSE]
  storage.mode(mm) <- "integer"
  x <- mm[, 1]
  y <- mm[, 2]
  z <- mm[, 3]
  in_bounds <- x >= 1L & x <= vol_dim[1] &
    y >= 1L & y <= vol_dim[2] &
    z >= 1L & z <= vol_dim[3]
  key <- x + (y - 1L) * vol_dim[1] + (z - 1L) * (vol_dim[1] * vol_dim[2])
  key[!in_bounds] <- NA_integer_
  key
}

ftq_qc_resolve_track_space <- function(track_space,
                                       track_affine,
                                       ref_image,
                                       track_voxel_size = NULL) {
  if (!identical(track_space, "auto")) return(track_space)
  if (!is.null(track_voxel_size)) return("voxelmm")
  if (!is.null(track_affine) && !is.null(ref_image)) return("trackvis")
  if (!is.null(ref_image)) return("world")
  "voxel"
}

ftq_qc_compiled_sampler_available <- function() {
  is.function(get0("ftq_qc_sample_chunk_cpp", mode = "function"))
}

ftq_qc_raw_to_voxel_affine <- function(track_space,
                                       ref_image,
                                       track_affine = NULL,
                                       track_voxel_size = NULL) {
  basis <- rbind(
    c(0, 0, 0),
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1)
  )
  mapped <- ftq_qc_map_track(
    track = basis,
    track_space = track_space,
    ref_image = ref_image,
    track_affine = track_affine,
    track_voxel_size = track_voxel_size
  )$voxel
  origin <- mapped[1, ]
  aff <- matrix(0, nrow = 4L, ncol = 4L)
  aff[1:3, 1] <- mapped[2, ] - origin
  aff[1:3, 2] <- mapped[3, ] - origin
  aff[1:3, 3] <- mapped[4, ] - origin
  aff[1:3, 4] <- origin
  aff[4, 4] <- 1
  aff
}

ftq_qc_map_track <- function(track,
                             track_space,
                             ref_image,
                             track_affine = NULL,
                             track_voxel_size = NULL,
                             streamline_voxelization_mode = "sample_points") {
  trk <- as_track(track)
  raw <- as.matrix(trk$coords)[, 1:3, drop = FALSE]
  if (identical(track_space, "voxelmm")) {
    voxel <- ftq_trackvis_voxelmm_to_voxel(
      raw,
      track_voxel_size,
      caller = "ftq_interactive_qc_viewer"
    )
    world <- ftq_voxel_to_world_mm(voxel, ref_image)
  } else if (identical(track_space, "vox10")) {
    voxel <- raw * 10
    world <- ftq_voxel_to_world_mm(voxel, ref_image)
  } else if (identical(track_space, "voxel")) {
    voxel <- raw
    world <- ftq_voxel_to_world_mm(voxel, ref_image)
  } else if (identical(track_space, "world")) {
    world <- raw
    voxel <- ftq_world_mm_to_voxel(world, ref_image)
  } else if (identical(track_space, "trackvis")) {
    if (is.null(track_affine)) stop("track_space = 'trackvis' requires a TrackVis affine.")
    world <- cbind(raw, 1) %*% t(track_affine)
    world <- world[, 1:3, drop = FALSE]
    voxel <- ftq_world_mm_to_voxel(world, ref_image)
  } else {
    stop("Unsupported track_space for QC viewer: ", track_space)
  }
  if (identical(streamline_voxelization_mode, "segment_interpolation")) {
    voxel <- ftq_densify_voxel_path(voxel)
    world <- ftq_voxel_to_world_mm(voxel, ref_image)
  } else if (!identical(streamline_voxelization_mode, "sample_points")) {
    stop("Unsupported streamline voxelization mode for QC viewer: ",
         streamline_voxelization_mode)
  }
  list(world = world, voxel = voxel)
}

ftq_qc_render_world <- function(voxel, ref_image, axis_flip, vol_dim) {
  render_voxel <- if (any(axis_flip)) {
    ftq_qc_apply_flip(voxel, axis_flip, vol_dim)
  } else {
    voxel
  }
  ftq_voxel_to_world_mm(render_voxel, ref_image)
}

ftq_qc_hit_detail <- function(mapped, lut, vol_dim, methods, axis_flip) {
  for (method in methods) {
    discrete <- switch(method, round = round(mapped), floor = floor(mapped), ceil = ceiling(mapped))
    discrete <- ftq_qc_apply_flip(discrete, axis_flip, vol_dim)
    key <- ftq_qc_encode_voxels(discrete, vol_dim)
    hit_indices <- which(!is.na(key) & lut[key])
    if (length(hit_indices)) {
      return(list(
        hit = TRUE,
        method = method,
        indices = hit_indices,
        discrete = discrete,
        keys = key
      ))
    }
  }
  list(
    hit = FALSE,
    method = NA_character_,
    indices = integer(),
    discrete = NULL,
    keys = integer()
  )
}

ftq_qc_hit_world <- function(hit, ref_image) {
  if (!isTRUE(hit$hit) || is.null(hit$discrete) || !length(hit$indices)) {
    return(matrix(numeric(), ncol = 3L))
  }
  ftq_voxel_to_world_mm(
    hit$discrete[hit$indices, , drop = FALSE],
    ref_image
  )
}

ftq_qc_apply_flip <- function(coords, axis_flip, vol_dim) {
  out <- as.matrix(coords)[, 1:3, drop = FALSE]
  if (axis_flip[1]) out[, 1] <- vol_dim[1] - out[, 1] + 1
  if (axis_flip[2]) out[, 2] <- vol_dim[2] - out[, 2] + 1
  if (axis_flip[3]) out[, 3] <- vol_dim[3] - out[, 3] + 1
  out
}

ftq_qc_endpoint_annotation <- function(mapped,
                                       focal_method,
                                       focal_lut,
                                       focal_lut_dilated,
                                       target_lut,
                                       target_lut_dilated,
                                       context_luts,
                                       axis_flip,
                                       vol_dim) {
  discrete <- ftq_endpoint_discretize(
    mapped$voxel,
    focal_method,
    axis_flip,
    vol_dim
  )
  keys <- ftq_qc_encode_voxels(discrete, vol_dim)
  two_endpoints <- nrow(discrete) > 1L
  endpoint_indices <- if (two_endpoints) c(1L, nrow(discrete)) else 1L
  endpoint_keys <- keys[endpoint_indices]
  start_key <- endpoint_keys[[1L]]
  end_key <- if (two_endpoints) endpoint_keys[[2L]] else start_key

  classification <- function(focal, target) {
    ftq_endpoint_classification(
      ftq_endpoint_membership(start_key, focal),
      ftq_endpoint_membership(end_key, focal),
      ftq_endpoint_membership(start_key, target),
      ftq_endpoint_membership(end_key, target),
      two_endpoints
    )
  }
  exact <- classification(focal_lut, target_lut)
  dilated1 <- classification(focal_lut_dilated, target_lut_dilated)
  class_names <- c(
    "shared_mask_ambiguous",
    "paired_focal_target_endpoints",
    "target_endpoint_only",
    "focal_endpoint_only",
    "neither_relevant_endpoint"
  )

  endpoint_status <- function(key, focal, target) {
    if (is.na(key)) return(list(status = "grid_outside", contexts = character()))
    in_focal <- ftq_endpoint_membership(key, focal)
    in_target <- ftq_endpoint_membership(key, target)
    contexts <- if (length(context_luts)) {
      names(context_luts)[vapply(
        context_luts,
        function(lut) ftq_endpoint_membership(key, lut),
        logical(1)
      )]
    } else {
      character()
    }
    status <- if (in_focal && in_target) {
      "shared"
    } else if (in_focal) {
      "focal"
    } else if (in_target) {
      "target"
    } else if (length(context_luts) && !length(contexts)) {
      "context_outside"
    } else {
      "other"
    }
    list(status = status, contexts = contexts)
  }

  endpoint_records <- lapply(seq_along(endpoint_indices), function(i) {
    index <- endpoint_indices[[i]]
    key <- endpoint_keys[[i]]
    exact_status <- endpoint_status(key, focal_lut, target_lut)
    d1_status <- endpoint_status(key, focal_lut_dilated, target_lut_dilated)
    list(
      end = if (!two_endpoints) "single" else if (i == 1L) "start" else "end",
      point = unname(as.numeric(mapped$render_world[index, 1:3])),
      statusExact = exact_status$status,
      statusDilated1 = d1_status$status,
      contexts = exact_status$contexts
    )
  })
  valid_keys <- unique(keys[!is.na(keys)])
  penetration <- if (length(valid_keys)) sum(target_lut[valid_keys]) else 0
  list(
    endpointClassExact = class_names[[exact$class]],
    endpointClassDilated1 = class_names[[dilated1$class]],
    pairedDilated1Only = !exact$paired && dilated1$paired,
    endpoints = endpoint_records,
    uniqueTargetVoxels = as.integer(penetration),
    discretization = focal_method
  )
}

ftq_qc_sample_record <- function(target,
                                 world,
                                 focal_indices,
                                 target_indices,
                                 focal_hit_world = NULL,
                                 target_hit_world = NULL,
                                 endpoint_annotation = NULL,
                                 max_points,
                                 focal_window) {
  n <- nrow(world)
  focal_indices <- sort(unique(focal_indices))
  target_indices <- sort(unique(target_indices))
  if (is.null(focal_hit_world)) {
    focal_hit_world <- world[focal_indices, , drop = FALSE]
  }
  if (is.null(target_hit_world)) {
    target_hit_world <- world[target_indices, , drop = FALSE]
  }
  focal_range <- unique(unlist(lapply(focal_indices, function(i) {
    seq.int(max(1L, i - focal_window), min(n, i + focal_window))
  }), use.names = FALSE))
  target_range <- unique(unlist(lapply(target_indices, function(i) {
    seq.int(max(1L, i - focal_window), min(n, i + focal_window))
  }), use.names = FALSE))
  after_range <- seq.int(min(focal_indices), n)
  list(
    target = target,
    points = ftq_qc_rows(ftq_qc_decimate(world, max_points)),
    focal = ftq_qc_rows(ftq_qc_decimate(world[focal_range, , drop = FALSE], max_points)),
    target_segment = ftq_qc_rows(ftq_qc_decimate(world[target_range, , drop = FALSE], max_points)),
    after = ftq_qc_rows(ftq_qc_decimate(world[after_range, , drop = FALSE], max_points)),
    hits = ftq_qc_rows(focal_hit_world),
    target_hits = ftq_qc_rows(target_hit_world),
    endpoint_annotation = endpoint_annotation
  )
}

ftq_qc_decimate <- function(points, n) {
  if (nrow(points) <= n) return(points)
  keep <- unique(round(seq(1, nrow(points), length.out = n)))
  points[keep, , drop = FALSE]
}

ftq_qc_rows <- function(x) {
  m <- as.matrix(x)[, 1:3, drop = FALSE]
  lapply(seq_len(nrow(m)), function(i) unname(as.numeric(m[i, ])))
}

ftq_qc_focal_boundary_coords <- function(mask, coords, max_points = 8000L) {
  coords <- as.matrix(coords)[, 1:3, drop = FALSE]
  storage.mode(coords) <- "integer"
  if (!nrow(coords)) {
    return(list(coords = coords, mode = "empty"))
  }
  dims <- dim(mask)[1:3]
  boundary_mask <- mask
  if (all(dims > 2L)) {
    interior <- array(FALSE, dim = dims)
    x_mid <- 2L:(dims[1] - 1L)
    y_mid <- 2L:(dims[2] - 1L)
    z_mid <- 2L:(dims[3] - 1L)
    interior[x_mid, y_mid, z_mid] <-
      mask[x_mid, y_mid, z_mid] &
      mask[x_mid - 1L, y_mid, z_mid] &
      mask[x_mid + 1L, y_mid, z_mid] &
      mask[x_mid, y_mid - 1L, z_mid] &
      mask[x_mid, y_mid + 1L, z_mid] &
      mask[x_mid, y_mid, z_mid - 1L] &
      mask[x_mid, y_mid, z_mid + 1L]
    boundary_mask <- mask & !interior
  }
  boundary <- which(boundary_mask, arr.ind = TRUE)
  mode <- "boundary"
  max_points <- as.integer(max_points)[1]
  if (!is.na(max_points) && max_points > 0L && nrow(boundary) > max_points) {
    keep <- unique(round(seq(1, nrow(boundary), length.out = max_points)))
    boundary <- boundary[keep, , drop = FALSE]
    mode <- "boundary_sample"
  }
  if (!nrow(boundary)) {
    boundary <- coords
    mode <- "all_voxels"
  }
  list(coords = boundary, mode = mode)
}

ftq_qc_target_centers <- function(targets) {
  do.call(rbind, lapply(names(targets), function(target) {
    pts <- ftq_voxel_to_world_mm(
      targets[[target]]$coords,
      targets[[target]]$image
    )
    center <- colMeans(pts)
    data.frame(target = target, x = center[1], y = center[2], z = center[3], stringsAsFactors = FALSE)
  }))
}

ftq_qc_target_boundaries <- function(targets, max_points = 2000L) {
  stats::setNames(lapply(names(targets), function(target) {
    roi <- targets[[target]]
    if (is.na(max_points) || max_points <= 0L) {
      points <- data.frame(x = numeric(), y = numeric(), z = numeric())
      return(list(
        target = target,
        voxel_count = nrow(roi$coords),
        rendered_count = 0L,
        render_mode = "hidden",
        points = points
      ))
    }
    boundary <- ftq_qc_focal_boundary_coords(roi$mask, roi$coords, max_points = max_points)
    world <- ftq_voxel_to_world_mm(boundary$coords, roi$image)
    list(
      target = target,
      voxel_count = nrow(roi$coords),
      rendered_count = nrow(boundary$coords),
      render_mode = boundary$mode,
      points = data.frame(x = world[, 1], y = world[, 2], z = world[, 3])
    )
  }), names(targets))
}

ftq_qc_background_cloud <- function(path, max_points = 12000L, context_mask = NULL) {
  if (is.null(path) || !file.exists(path) || is.na(max_points) || max_points <= 0L) {
    return(data.frame(x = numeric(), y = numeric(), z = numeric()))
  }
  img <- RNifti::readNifti(path, internal = FALSE)
  arr <- as.array(img)
  candidate <- is.finite(arr)
  if (!is.null(context_mask) && identical(dim(context_mask)[1:3], dim(arr)[1:3])) {
    seed <- context_mask & is.finite(arr) & arr > 0
    vals <- arr[seed]
    if (length(vals)) {
      threshold <- stats::quantile(vals, 0.05, na.rm = TRUE)
      thresholded <- context_mask & is.finite(arr) & arr > threshold
      candidate <- if (any(thresholded)) thresholded else seed
    } else {
      candidate <- context_mask & is.finite(arr)
    }
  } else {
    vals <- arr[candidate & arr > 0]
    if (!length(vals)) return(data.frame(x = numeric(), y = numeric(), z = numeric()))
    threshold <- stats::quantile(vals, 0.20, na.rm = TRUE)
    thresholded <- candidate & arr > threshold
    if (any(thresholded)) candidate <- thresholded
  }
  ijk <- which(candidate, arr.ind = TRUE)
  if (nrow(ijk) > max_points) {
    keep <- ftq_qc_evenly_spaced_indices(nrow(ijk), max_points)
    ijk <- ijk[keep, , drop = FALSE]
  }
  world <- ftq_voxel_to_world_mm(ijk, img)
  data.frame(x = world[, 1], y = world[, 2], z = world[, 3])
}

ftq_qc_context_cloud <- function(background_image = NULL,
                                 tractogram_context_image = NULL,
                                 context_mask = NULL,
                                 max_points = 12000L) {
  if (!is.null(background_image) && file.exists(background_image)) {
    return(list(
      points = ftq_qc_background_cloud(
        background_image,
        max_points = max_points,
        context_mask = context_mask
      ),
      label = if (!is.null(context_mask)) {
        "anatomical image foreground calibrated from local/QC mask"
      } else {
        "anatomical image"
      }
    ))
  }
  if (!is.null(tractogram_context_image) && file.exists(tractogram_context_image)) {
    return(list(
      points = ftq_qc_density_cloud(tractogram_context_image, max_points = max_points),
      label = "Tractogram-derived density context"
    ))
  }
  list(
    points = data.frame(x = numeric(), y = numeric(), z = numeric()),
    label = "none"
  )
}

ftq_qc_density_cloud <- function(path, max_points = 12000L) {
  if (is.null(path) || !file.exists(path) || is.na(max_points) || max_points <= 0L) {
    return(data.frame(x = numeric(), y = numeric(), z = numeric()))
  }
  img <- RNifti::readNifti(path, internal = FALSE)
  arr <- as.array(img)
  vals <- arr[is.finite(arr) & arr > 0]
  if (!length(vals)) return(data.frame(x = numeric(), y = numeric(), z = numeric()))
  ijk <- which(is.finite(arr) & arr > 0, arr.ind = TRUE)
  if (nrow(ijk) > max_points) {
    keep <- ftq_qc_evenly_spaced_indices(nrow(ijk), max_points)
    ijk <- ijk[keep, , drop = FALSE]
  }
  world <- ftq_voxel_to_world_mm(ijk, img)
  data.frame(x = world[, 1], y = world[, 2], z = world[, 3])
}

ftq_qc_evenly_spaced_indices <- function(n, max_points) {
  n <- as.integer(n)[1]
  max_points <- as.integer(max_points)[1]
  if (is.na(n) || n < 1L || is.na(max_points) || max_points < 1L) {
    return(integer())
  }
  if (n <= max_points) return(seq_len(n))
  as.integer(round(seq.int(1, n, length.out = max_points)))
}

ftq_qc_write_html <- function(path,
                              streamlines,
                              focal_world,
                              focal_total_count,
                              focal_rendered_count,
                              focal_render_mode,
                              target_centers,
                              target_boundaries,
                              background_world,
                              context_label,
                              colors,
                              scanned,
                              focal_hits,
                              sample_counts,
                              endpoint_analysis,
                              title,
                              target_intersection_mode = "full_streamline",
                              local_search_radius_mm = 2) {
  target_intersection_mode <- ftq_normalize_target_intersection_mode(
    target_intersection_mode
  )
  assignment_sample_label <- switch(
    target_intersection_mode,
    full_streamline = "intersection-positive",
    endpoint_voxels = "endpoint-voxel-assigned",
    local_search = "local-search-assigned"
  )
  assignment_detail <- switch(
    target_intersection_mode,
    full_streamline = "any selected path voxel intersects the target ROI",
    endpoint_voxels =
      "the endpoint voxels form an unambiguous focal-target pair",
    local_search = paste0(
      "the endpoints form an unambiguous focal-target pair after nearest-",
      "label voxel-centre searches within ", local_search_radius_mm, " mm"
    )
  )
  label_colors <- stats::setNames(unname(colors), gsub("_", " ", names(colors)))
  streamline_payload <- lapply(seq_along(streamlines), function(i) {
    x <- streamlines[[i]]
    label <- gsub("_", " ", x$target)
    endpoint <- x$endpoint_annotation
    list(
      id = paste0(label, " ", i),
      target = label,
      color = unname(label_colors[[label]]),
      points = x$points,
      focal = x$focal,
      targetSegment = x$target_segment,
      after = x$after,
      hits = x$hits,
      targetHits = x$target_hits,
      endpointClassExact = if (is.null(endpoint)) "unclassified" else endpoint$endpointClassExact,
      endpointClassDilated1 = if (is.null(endpoint)) "unclassified" else endpoint$endpointClassDilated1,
      pairedDilated1Only = if (is.null(endpoint)) FALSE else endpoint$pairedDilated1Only,
      endpoints = if (is.null(endpoint)) list() else endpoint$endpoints,
      uniqueTargetVoxels = if (is.null(endpoint)) NA_integer_ else endpoint$uniqueTargetVoxels,
      endpointDiscretization = if (is.null(endpoint)) NA_character_ else endpoint$discretization
    )
  })
  focal_df <- data.frame(x = focal_world[, 1], y = focal_world[, 2], z = focal_world[, 3])
  target_payload <- lapply(seq_len(nrow(target_centers)), function(i) {
    target <- target_centers$target[i]
    label <- gsub("_", " ", target)
    boundary <- target_boundaries[[target]]
    list(
      target = label,
      color = unname(label_colors[[label]]),
      point = unname(as.numeric(target_centers[i, c("x", "y", "z")])),
      boundary = ftq_qc_rows(boundary$points),
      voxelCount = boundary$voxel_count,
      renderedCount = boundary$rendered_count,
      renderMode = boundary$render_mode
    )
  })
  # Each display mode is decimated separately, so include its points and hit
  # markers when framing the scene, not only the whole-streamline sample.
  stream_world <- do.call(rbind, lapply(streamlines, function(x) {
    points <- unlist(x[c("points", "focal", "target_segment", "after",
                         "hits", "target_hits")], recursive = FALSE)
    do.call(rbind, points)
  }))
  stream_world <- data.frame(x = stream_world[, 1], y = stream_world[, 2], z = stream_world[, 3])
  target_extent_world <- do.call(rbind, lapply(target_boundaries, function(x) {
    x$points[, c("x", "y", "z"), drop = FALSE]
  }))
  if (is.null(target_extent_world)) {
    target_extent_world <- data.frame(x = numeric(), y = numeric(), z = numeric())
  }
  all_world <- rbind(
    stream_world,
    focal_df,
    data.frame(target_centers[, c("x", "y", "z")]),
    target_extent_world,
    data.frame(background_world[, c("x", "y", "z"), drop = FALSE])
  )
  bounds <- apply(all_world, 2, range, finite = TRUE)
  center <- colMeans(bounds)
  # A sphere encloses the scene after any rotation. Half the longest box side
  # can crop diagonal paths, especially when the canvas is narrow.
  radius <- max(sqrt(rowSums(sweep(as.matrix(all_world), 2, center, "-")^2)))
  if (!is.finite(radius) || radius <= 0) radius <- 1
  total_sampled <- length(streamline_payload)
  sampled_by_target <- paste(
    paste0(names(sample_counts), " n=", format(as.integer(sample_counts), big.mark = ",")),
    collapse = "; "
  )
  payload <- list(
    metadata = list(
      title = title,
      subtitle = paste0(
        "Target-balanced quota QC subset: ", format(total_sampled, big.mark = ","),
        " ", assignment_sample_label, " streamline examples. ",
        "This is not a random or prevalence-weighted sample. ",
        "Use the full-run tables and maps for quantitative results."
      ),
      method = paste0("Target assignment: ", target_intersection_mode,
                      ": ", assignment_detail, "."),
      context = if (identical(context_label, "Tractogram-derived density context")) {
        "Context: tractogram-derived density; not anatomical MRI or brain volume."
      } else if (grepl("^anatomical image", context_label)) {
        "Context: supplied anatomical image."
      } else {
        "Context: unavailable; anatomical alignment cannot be assessed from a background image."
      },
      interaction = "Drag to rotate freely; drag near an edge to roll. Scroll or pinch to zoom. Reset restores the initial view. The compass follows the scene.",
      note = paste0(
        "The initial anatomical view places superior upward; the R/L, A/P, and S/I compass rotates with the scene. Drag to rotate freely, drag near an edge to roll, and scroll or use two fingers to pinch to zoom. Reset restores the initial orientation and zoom. All ",
        format(focal_total_count, big.mark = ","),
        " focal ROI voxels are represented by a translucent boundary",
        if (identical(focal_render_mode, "boundary_sample")) {
          paste0(" sample of ", format(focal_rendered_count, big.mark = ","), " points")
        } else {
          ""
        },
        ". ",
        "Displayed streamlines are a target-balanced quota subset for spatial QC; quantitative tables and maps come from the full FibrTraQ counting run. ",
        "Examples were accepted in tractogram order until target quotas were filled or scanning stopped. ",
        "The sampler scanned ", format(scanned, big.mark = ","), " streamlines and observed ",
        format(focal_hits, big.mark = ","), " focal hits. ",
        if (nzchar(sampled_by_target)) paste0("Displayed per target: ", sampled_by_target, ". ") else "",
        "Endpoint markers are drawn at continuous streamline endpoints; their colours use voxel-based classification. ",
        "Circle and square mark the first and last stored path points, not biological direction; a diamond marks a single-point path. ",
        "After Focal Hit follows stored path order and does not imply signal direction. ",
        "The selected target-assignment method is ",
        target_intersection_mode,
        ": ", assignment_detail, ". ",
        if (isTRUE(endpoint_analysis)) {
          paste0(
            "End-voxel classes describe this same ",
            assignment_sample_label,
            " QC sample and do not rebalance sample selection; the audit itself retains full-streamline target-intersection denominators. ",
            "Exact end voxels are primary for this audit; 1-voxel sensitivity uses a 26-neighbor voxel dilation and is distinct from physical-distance local search. "
          )
        } else {
          ""
        },
        "Whole Streamline mode is intentionally unclipped and can include raw tractography points outside the displayed foreground context. ",
        "Focal-hit markers show where sampled streamlines intersected the focal ROI. ",
        "Target ROI extents and target-hit markers show the path intersection or assigned endpoint voxel centre used by the selected method, ",
        "which is important for large target regions. Target centroid labels are orientation aids only. ",
        if (identical(context_label, "Tractogram-derived density context")) {
          "Background points are Tractogram-derived density context, not anatomical MRI or brain volume."
        } else if (grepl("^anatomical image", context_label)) {
          "Background points are sampled from the supplied anatomical image."
        } else {
          "No background context is shown."
        }
      )
    ),
    bounds = list(center = unname(as.numeric(center)), radius = unname(as.numeric(radius))),
    brain = ftq_qc_rows(background_world),
    contextLabel = context_label,
    focal = ftq_qc_rows(focal_df),
    focalCount = focal_total_count,
    focalRenderedCount = focal_rendered_count,
    focalRenderMode = focal_render_mode,
    endpointAnalysis = isTRUE(endpoint_analysis),
    targetIntersectionMode = target_intersection_mode,
    localSearchRadiusMm = if (identical(
      target_intersection_mode,
      "local_search"
    )) {
      local_search_radius_mm
    } else {
      NULL
    },
    assignmentSampleLabel = assignment_sample_label,
    samplingStrategy = "target_quota_scan_order",
    endpointClasses = c(
      "paired_focal_target_endpoints",
      "shared_mask_ambiguous",
      "target_endpoint_only",
      "focal_endpoint_only",
      "neither_relevant_endpoint"
    ),
    targets = target_payload,
    streamlines = streamline_payload,
    legend = lapply(names(label_colors), function(label) list(target = label, color = unname(label_colors[[label]])))
  )
  json <- jsonlite::toJSON(payload, auto_unbox = TRUE, digits = 7, pretty = FALSE)
  # A script element is terminated by the HTML parser before JSON parsing. Encode
  # markup-significant characters so an untrusted target label cannot close the
  # application/json element and inject executable markup into the QC viewer.
  json <- as.character(json)
  escape_script_character <- function(value, literal, unicode_code) {
    paste(
      strsplit(value, literal, fixed = TRUE)[[1L]],
      collapse = paste0(intToUtf8(92L), "u", unicode_code)
    )
  }
  json <- escape_script_character(json, "&", "0026")
  json <- escape_script_character(json, "<", "003c")
  json <- escape_script_character(json, ">", "003e")
  guide <- ftq_write_output_guide(dirname(path))
  writeLines(ftq_qc_html_template(
    json, guide_href = ftq_report_href(guide, dirname(path))
  ), path)
  invisible(path)
}

ftq_qc_html_template <- function(json, guide_href = "OUTPUTS_README.html") {
  paste0('<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Interactive 3D FibrTraQ QC</title>
<style>
[hidden]{display:none!important}
body{margin:0;background:#f8fafc;color:#111827;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
main{width:min(1500px,calc(100vw - 28px));margin:14px auto}
h1{margin:0 0 4px;font-size:24px;line-height:1.15}p{margin:0;color:#475569;font-size:14px}
.toolbar{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin:12px 0}
button{border:1px solid #cbd5e1;background:white;border-radius:6px;padding:7px 10px;color:#111827;cursor:pointer}
button[aria-pressed="true"]{background:#e0f2fe;border-color:#0284c7}
.advanced-controls{width:100%;border-top:1px solid #e2e8f0;padding-top:8px;margin-top:2px}
.advanced-controls summary{cursor:pointer;font-size:13px;font-weight:700;color:#334155;width:max-content}
.advanced-panel,.endpoint-controls{display:flex;flex-wrap:wrap;gap:8px;align-items:center;padding-top:8px}
.endpoint-controls{width:100%;border-top:1px solid #e2e8f0;margin-top:2px}
.control-label{font-size:12px;font-weight:700;color:#475569;margin-right:2px}.class-filters{display:flex;flex-wrap:wrap;gap:6px;align-items:center}
.class-summary{width:100%;min-height:20px;margin:0;color:#475569;font-size:12px}
.viewer{position:relative;height:min(78vh,900px);min-height:560px;border:1px solid #cbd5e1;border-radius:8px;overflow:hidden;background:white}
canvas{display:block;width:100%;height:100%;cursor:grab;touch-action:none}.dragging{cursor:grabbing}
.hud{position:absolute;left:12px;top:12px;max-width:380px;padding:9px 10px;border-radius:7px;background:rgba(255,255,255,.88);border:1px solid rgba(203,213,225,.9);font-size:12px;color:#334155}
.viewer-details{margin-top:14px;padding:12px;border:1px solid #cbd5e1;border-radius:7px;background:white}.viewer-details summary{cursor:pointer;font-weight:700}.viewer-details p{margin-top:8px}.method-context{margin-top:7px}.endpoint-legend{margin-top:12px;padding:10px;background:#fff;border:1px solid #cbd5e1;border-radius:7px}.endpoint-legend p{font-size:12px;margin-top:5px}.endpoint-status{font-weight:700}.marker-key{display:flex;flex-wrap:wrap;gap:5px 16px;font-size:12px}.marker-key span{display:inline-flex;align-items:center;gap:5px}.marker-dot{display:inline-block;width:10px;height:10px;border-radius:50%;background:var(--marker-color)}
.legend{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:6px 14px;margin-top:10px;font-size:13px;color:#334155}
.legend span{display:inline-flex;align-items:center;gap:7px;min-width:0}.swatch{width:23px;height:4px;border-radius:99px;display:inline-block;flex:0 0 auto}
@media (max-width:600px){.viewer{display:flex;flex-direction:column;height:auto;min-height:0;overflow:hidden}.hud{position:static;order:-1;max-width:none;box-sizing:border-box;border-width:0 0 1px;border-radius:0;background:#fff}.viewer canvas{height:520px;flex:0 0 520px}}
</style>
</head>
<body>
<main>
<h1 id="title"></h1><p id="subtitle"></p>
<p class="method-context" id="methodOverview"></p><p id="contextOverview"></p>
<div class="toolbar">
<button id="reset" type="button">Reset Orientation</button>
<button id="toggleBrain" type="button" aria-pressed="true">Spatial Context</button>
<button id="toggleFocal" type="button" aria-pressed="true">fNode Boundary</button>
<button id="toggleTargets" type="button" aria-pressed="true">cNode Extents</button>
<button id="toggleEndpoints" type="button" aria-pressed="true">Endpoint Markers</button>
<details class="advanced-controls" id="advancedControls">
<summary>Advanced controls</summary>
<div class="advanced-panel">
<button id="toggleHits" type="button" aria-pressed="false">fNode-Hit Markers</button>
<button id="toggleTargetHits" type="button" aria-pressed="false">cNode-Hit Markers</button>
<button id="toggleLabels" type="button" aria-pressed="false">cNode Labels</button>
<span class="control-label">Path display</span>
<button class="segment" data-mode="points" type="button" aria-pressed="true">Whole Streamline</button>
<button class="segment" data-mode="focal" type="button" aria-pressed="false">Focal Segment</button>
<button class="segment" data-mode="targetSegment" type="button" aria-pressed="false">Target Segment</button>
<button class="segment" data-mode="after" type="button" aria-pressed="false">After Focal Hit</button>
<div class="endpoint-controls" id="endpointControls">
<span class="control-label">Endpoint classification</span>
<button class="endpoint-mode" data-endpoint-mode="Exact" type="button" aria-pressed="true">Exact</button>
<button class="endpoint-mode" data-endpoint-mode="Dilated1" type="button" aria-pressed="false">1-voxel sensitivity</button>
<span class="control-label">Classes</span><span class="class-filters" id="classFilters"></span>
</div>
<div class="class-summary" id="classSummary"></div>
</div>
</details>
</div>
<section class="viewer" aria-label="Interactive 3D streamline quality-control viewer"><canvas id="scene"></canvas><div class="hud" id="hud"></div></section>
<div class="legend" id="legend" aria-label="Target colours"></div>
<section class="endpoint-legend" id="endpointLegend" aria-label="Endpoint marker legend">
<p class="endpoint-status" id="endpointStatus" aria-live="polite"></p>
<p>Markers sit at continuous streamline endpoints; colours describe their voxel-based classification.</p>
<p>Shapes: circle = first stored point; square = last stored point; diamond = single-point path. Stored order does not show biological direction.</p>
<div class="marker-key">
<span><i class="marker-dot" style="--marker-color:#f97316"></i>Orange: focal ROI</span>
<span>Target colour: target ROI</span>
<span><i class="marker-dot" style="--marker-color:#c24156"></i>Rose: shared focal/target mask</span>
<span><i class="marker-dot" style="--marker-color:#111827"></i>Black: outside reference grid</span>
<span><i class="marker-dot" style="--marker-color:#dc2626"></i>Red: outside supplied endpoint contexts, after focal/target checks</span>
<span><i class="marker-dot" style="--marker-color:#64748b"></i>Grey: other or unavailable context</span>
</div>
</section>
<p class="method-context"><a href="', ftq_report_html_escape(guide_href), '">Output definitions and guide</a></p>
<details class="viewer-details" id="viewerDetails"><summary>Sampling, methods and interpretation</summary><p id="viewerMethods"></p></details>
</main>
<script id="scene-data" type="application/json">', json, '</script>
<script>
const data=JSON.parse(document.getElementById("scene-data").textContent);
const canvas=document.getElementById("scene"),ctx=canvas.getContext("2d");
const elTitle=document.getElementById("title"),elSubtitle=document.getElementById("subtitle"),elHud=document.getElementById("hud"),elLegend=document.getElementById("legend"),elEndpointControls=document.getElementById("endpointControls"),elClassFilters=document.getElementById("classFilters"),elClassSummary=document.getElementById("classSummary");
const btnReset=document.getElementById("reset"),btnBrain=document.getElementById("toggleBrain"),btnFocal=document.getElementById("toggleFocal"),btnTargets=document.getElementById("toggleTargets"),btnHits=document.getElementById("toggleHits"),btnTargetHits=document.getElementById("toggleTargetHits"),btnLabels=document.getElementById("toggleLabels"),btnEndpoints=document.getElementById("toggleEndpoints");
const classLabels={paired_focal_target_endpoints:"Paired fNode-cNode ends",shared_mask_ambiguous:"Shared-mask ambiguity",target_endpoint_only:"cNode end only",focal_endpoint_only:"fNode end only",neither_relevant_endpoint:"Neither"};
// Unit quaternions allow continuous rotation through every anatomical view.
function multiplyQuaternion(a,b){const [x,y,z,w]=a,[u,v,t,s]=b;return[w*u+x*s+y*t-z*v,w*v-x*t+y*s+z*u,w*t+x*v-y*u+z*s,w*s-x*u-y*v-z*t]}
function normalizeQuaternion(q){const n=Math.hypot(...q);return n?q.map(v=>v/n):[0,0,0,1]}
function initialOrientation(){return multiplyQuaternion([Math.sin(.12),0,0,Math.cos(.12)],[0,Math.sin(-.31),0,Math.cos(-.31)])}
function rotationMatrix(q){const [x,y,z,w]=q;return[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w),2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w),2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]}
const state={orientation:initialOrientation(),zoom:1,showBrain:true,showFocal:true,showTargets:true,showHits:false,showTargetHits:false,showLabels:false,showEndpoints:true,segment:"points",endpointMode:"Exact",activeClasses:new Set(data.endpointClasses||[])};
const pointers=new Map();
let gesture=null,frame=null,renderPending=false;
elTitle.textContent=data.metadata.title;elSubtitle.textContent=data.metadata.subtitle;elHud.textContent=data.metadata.interaction;document.getElementById("methodOverview").textContent=data.metadata.method;document.getElementById("contextOverview").textContent=data.metadata.context;document.getElementById("viewerMethods").textContent=data.metadata.note;
for(const d of data.legend){const row=document.createElement("span"),swatch=document.createElement("i");swatch.className="swatch";swatch.style.backgroundColor=typeof d.color==="string"?d.color:"";row.appendChild(swatch);row.appendChild(document.createTextNode(String(d.target)));elLegend.appendChild(row)}
function pressed(id,v){document.getElementById(id).setAttribute("aria-pressed",v?"true":"false")}
function updateEndpointStatus(){const legend=document.getElementById("endpointLegend"),status=document.getElementById("endpointStatus");if(!data.endpointAnalysis){status.textContent="Endpoint classification was not generated for this run.";legend.querySelectorAll("p:not(.endpoint-status),.marker-key").forEach(e=>e.hidden=true);return}status.textContent=`Endpoint markers: ${state.showEndpoints?"shown":"hidden"}. Classification: ${state.endpointMode==="Exact"?"exact endpoint voxels":"1-voxel sensitivity (26-neighbour dilation)"}. Showing ${data.streamlines.filter(visible).length} of ${data.streamlines.length} sampled streamlines.`}
function endpointClass(s){return s[`endpointClass${state.endpointMode}`]||"unclassified"}
function visible(s){return !data.endpointAnalysis||state.activeClasses.has(endpointClass(s))}
	function updateClassSummary(){if(!data.endpointAnalysis){elEndpointControls.hidden=true;elClassSummary.hidden=true;btnEndpoints.hidden=true;return}const counts={};for(const c of data.endpointClasses)counts[c]=0;for(const s of data.streamlines){const c=endpointClass(s);counts[c]=(counts[c]||0)+1}const shown=data.streamlines.filter(visible).length;elClassSummary.textContent=`${state.endpointMode==="Exact"?"Exact":"1-voxel sensitivity"} end-voxel classes in this ${data.assignmentSampleLabel||"target-positive"} quota subset: ${data.endpointClasses.map(c=>`${classLabels[c]} ${counts[c]||0}`).join("; ")}. Showing ${shown} of ${data.streamlines.length}.`}
if(data.endpointAnalysis){for(const c of data.endpointClasses){const b=document.createElement("button");b.type="button";b.dataset.endpointClass=c;b.textContent=classLabels[c];b.setAttribute("aria-pressed","true");b.onclick=()=>{if(state.activeClasses.has(c)){state.activeClasses.delete(c);b.setAttribute("aria-pressed","false")}else{state.activeClasses.add(c);b.setAttribute("aria-pressed","true")}updateClassSummary();updateEndpointStatus();render()};elClassFilters.appendChild(b)}}else{elEndpointControls.hidden=true;elClassSummary.hidden=true}
btnReset.onclick=()=>{clearGesture();state.orientation=initialOrientation();state.zoom=1;render()};
btnBrain.onclick=()=>{state.showBrain=!state.showBrain;pressed("toggleBrain",state.showBrain);render()};
btnFocal.onclick=()=>{state.showFocal=!state.showFocal;pressed("toggleFocal",state.showFocal);render()};
btnTargets.onclick=()=>{state.showTargets=!state.showTargets;pressed("toggleTargets",state.showTargets);render()};
btnHits.onclick=()=>{state.showHits=!state.showHits;pressed("toggleHits",state.showHits);render()};
btnTargetHits.onclick=()=>{state.showTargetHits=!state.showTargetHits;pressed("toggleTargetHits",state.showTargetHits);render()};
btnLabels.onclick=()=>{state.showLabels=!state.showLabels;pressed("toggleLabels",state.showLabels);render()};
document.querySelectorAll(".segment").forEach(b=>b.onclick=()=>{state.segment=b.dataset.mode;document.querySelectorAll(".segment").forEach(x=>x.setAttribute("aria-pressed",x===b?"true":"false"));render()});
btnEndpoints.onclick=()=>{state.showEndpoints=!state.showEndpoints;pressed("toggleEndpoints",state.showEndpoints);updateEndpointStatus();render()};
document.querySelectorAll(".endpoint-mode").forEach(b=>b.onclick=()=>{state.endpointMode=b.dataset.endpointMode;document.querySelectorAll(".endpoint-mode").forEach(x=>x.setAttribute("aria-pressed",x===b?"true":"false"));updateClassSummary();updateEndpointStatus();render()});
function clearGesture(){const ids=[...pointers.keys()];pointers.clear();gesture=null;canvas.classList.remove("dragging");for(const id of ids)if(canvas.hasPointerCapture(id))canvas.releasePointerCapture(id)}
function resize(){clearGesture();const r=canvas.getBoundingClientRect(),d=window.devicePixelRatio||1;canvas.width=Math.max(1,Math.round(r.width*d));canvas.height=Math.max(1,Math.round(r.height*d));ctx.setTransform(d,0,0,d,0,0);render()}
window.onresize=resize;
// Anchor each drag to its starting pose so event frequency cannot change the result.
function trackball(p,r){const radius=Math.max(1,Math.min(r.width,r.height)/2.25),x=(p.x-r.left-r.width/2)/radius,y=(r.top+r.height/2-p.y)/radius,d=x*x+y*y,z=d<=.5?Math.sqrt(1-d):.5/Math.sqrt(d),n=Math.hypot(x,y,z);return[x/n,y/n,z/n]}
function rotationBetween(a,b){const cross=[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],w=1+a[0]*b[0]+a[1]*b[1]+a[2]*b[2];if(w<1e-10){const axis=Math.abs(a[0])<.9?[0,a[2],-a[1]]:[-a[2],0,a[0]];return normalizeQuaternion([...axis,0])}return normalizeQuaternion([...cross,w])}
function pointerSpan(){const p=[...pointers.values()];return Math.hypot(p[1].x-p[0].x,p[1].y-p[0].y)}
function beginGesture(){if(!pointers.size){gesture=null;canvas.classList.remove("dragging");return}const r=canvas.getBoundingClientRect();gesture={rect:r,orientation:state.orientation.slice(),start:trackball(pointers.values().next().value,r),span:pointers.size===2?pointerSpan():0,zoom:state.zoom};canvas.classList.add("dragging")}
function boundedZoom(value){return Math.max(.35,Math.min(5,value))}
canvas.onpointerdown=e=>{if(e.button!==0||pointers.has(e.pointerId)||pointers.size>=2)return;if(pointers.size&&(e.pointerType!=="touch"||[...pointers.values()].some(p=>p.type!=="touch")))return;pointers.set(e.pointerId,{x:e.clientX,y:e.clientY,type:e.pointerType});canvas.setPointerCapture(e.pointerId);beginGesture()};
canvas.onpointermove=e=>{if(!pointers.has(e.pointerId)||!gesture)return;pointers.set(e.pointerId,{x:e.clientX,y:e.clientY,type:e.pointerType});if(pointers.size===2){if(gesture.span>0)state.zoom=boundedZoom(gesture.zoom*pointerSpan()/gesture.span);else beginGesture()}else{const end=trackball(pointers.get(e.pointerId),gesture.rect);state.orientation=normalizeQuaternion(multiplyQuaternion(rotationBetween(gesture.start,end),gesture.orientation))}render()};
function endPointer(e){if(!pointers.delete(e.pointerId))return;if(canvas.hasPointerCapture(e.pointerId))canvas.releasePointerCapture(e.pointerId);beginGesture()}
canvas.onpointerup=endPointer;canvas.onpointercancel=endPointer;canvas.onlostpointercapture=endPointer;
canvas.addEventListener("wheel",e=>{e.preventDefault();const unit=e.deltaMode===1?16:e.deltaMode===2?canvas.getBoundingClientRect().height:1;state.zoom=boundedZoom(state.zoom*Math.exp(-e.deltaY*unit*.001));if(pointers.size)beginGesture();render()},{passive:false});
// Reuse one viewport measurement and rotation matrix throughout each frame.
function project(pt){const c=data.bounds.center;let x=pt[0]-c[0],y=pt[2]-c[2],z=-(pt[1]-c[1]);const m=frame.matrix,r=frame.rect,s=frame.scale;return{x:r.width/2+(m[0]*x+m[1]*y+m[2]*z)*s,y:r.height/2-(m[3]*x+m[4]*y+m[5]*z)*s,z:m[6]*x+m[7]*y+m[8]*z}}
function rgba(hex,a){const h=hex.replace("#","");return`rgba(${parseInt(h.slice(0,2),16)},${parseInt(h.slice(2,4),16)},${parseInt(h.slice(4,6),16)},${a})`}
function cloud(pts,color,a,rad){ctx.fillStyle=rgba(color,a);for(const pt of pts){const p=project(pt);ctx.beginPath();ctx.arc(p.x,p.y,rad,0,Math.PI*2);ctx.fill()}}
function focalBoundary(){cloud(data.focal,"#f97316",.30,1.65)}
function targetRegions(){for(const t of data.targets){if(t.boundary&&t.boundary.length)cloud(t.boundary,t.color,.20,1.35)}}
function lines(){for(const s of data.streamlines){if(!visible(s))continue;const pts=s[state.segment]||s.points;if(!pts.length)continue;ctx.strokeStyle=rgba(s.color,.74);ctx.lineWidth=state.segment==="focal"?2.6:1.8;ctx.lineJoin="round";ctx.lineCap="round";ctx.beginPath();pts.forEach((pt,i)=>{const p=project(pt);if(i===0)ctx.moveTo(p.x,p.y);else ctx.lineTo(p.x,p.y)});ctx.stroke()}}
function hits(){for(const s of data.streamlines){if(!visible(s))continue;ctx.fillStyle="white";ctx.strokeStyle=s.color;ctx.lineWidth=2.3;for(const h of s.hits){const p=project(h);ctx.beginPath();ctx.arc(p.x,p.y,4.4,0,Math.PI*2);ctx.fill();ctx.stroke()}}}
function diamond(p,r){ctx.beginPath();ctx.moveTo(p.x,p.y-r);ctx.lineTo(p.x+r,p.y);ctx.lineTo(p.x,p.y+r);ctx.lineTo(p.x-r,p.y);ctx.closePath()}
function targetHits(){for(const s of data.streamlines){if(!visible(s)||!s.targetHits||!s.targetHits.length)continue;ctx.fillStyle=rgba(s.color,.94);ctx.strokeStyle="#0f172a";ctx.lineWidth=1.2;for(const h of s.targetHits){const p=project(h);diamond(p,5.2);ctx.fill();ctx.stroke()}}}
function endpointMarkers(){const colors={focal:"#f97316",shared:"#c24156",grid_outside:"#111827",context_outside:"#dc2626",other:"#64748b"};for(const s of data.streamlines){if(!visible(s)||!s.endpoints)continue;for(const e of s.endpoints){const p=project(e.point),status=e[`status${state.endpointMode}`]||"other";ctx.fillStyle=status==="target"?s.color:(colors[status]||colors.other);ctx.strokeStyle="#ffffff";ctx.lineWidth=1.8;if(e.end==="end"){ctx.beginPath();ctx.rect(p.x-5,p.y-5,10,10)}else if(e.end==="single"){diamond(p,5.8)}else{ctx.beginPath();ctx.arc(p.x,p.y,5.2,0,Math.PI*2)}ctx.fill();ctx.stroke()}}}
function targets(){ctx.font="12px -apple-system,BlinkMacSystemFont,Segoe UI,sans-serif";ctx.textBaseline="middle";const projected=data.targets.map(t=>({t,p:project(t.point)}));for(const d of projected){ctx.fillStyle="white";ctx.strokeStyle=d.t.color;ctx.lineWidth=3;ctx.beginPath();ctx.arc(d.p.x,d.p.y,7,0,Math.PI*2);ctx.fill();ctx.stroke()}if(!state.showLabels)return;const r=frame.rect,boxes=[];for(const d of projected){const label=`${d.t.target} centroid`,w=ctx.measureText(label).width+10,h=20;let x=d.p.x+10;if(x+w>r.width-4)x=d.p.x-w-10;x=Math.max(4,Math.min(r.width-w-4,x));const offsets=[-10,12,-32,34,-54,56];let y=Math.max(2,Math.min(r.height-h-2,d.p.y-10));for(const offset of offsets){const candidate=Math.max(2,Math.min(r.height-h-2,d.p.y+offset));const box={left:x,right:x+w,top:candidate,bottom:candidate+h};if(boxes.every(b=>box.right<=b.left||box.left>=b.right||box.bottom<=b.top||box.top>=b.bottom)){y=candidate;break}}boxes.push({left:x,right:x+w,top:y,bottom:y+h});ctx.fillStyle="rgba(255,255,255,.88)";ctx.fillRect(x,y,w,h);ctx.fillStyle="#111827";ctx.fillText(label,x+5,y+h/2)}}
function axes(){const c=data.bounds.center,r=data.bounds.radius*.32,o=project(c),box=frame.rect,origin={x:64,y:box.height-64},dirs=[["R",[c[0]+r,c[1],c[2]],"#b91c1c"],["L",[c[0]-r,c[1],c[2]],"#b91c1c"],["A",[c[0],c[1]+r,c[2]],"#15803d"],["P",[c[0],c[1]-r,c[2]],"#15803d"],["S",[c[0],c[1],c[2]+r],"#0369a1"],["I",[c[0],c[1],c[2]-r],"#0369a1"]];ctx.font="bold 11px -apple-system,BlinkMacSystemFont,Segoe UI,sans-serif";ctx.textAlign="center";ctx.textBaseline="middle";ctx.fillStyle="rgba(255,255,255,.88)";ctx.beginPath();ctx.arc(origin.x,origin.y,50,0,Math.PI*2);ctx.fill();for(const d of dirs){const q=project(d[1]),dx=q.x-o.x,dy=q.y-o.y,n=Math.hypot(dx,dy)||1,ux=dx/n,uy=dy/n,end={x:origin.x+ux*34,y:origin.y+uy*34};ctx.strokeStyle=d[2];ctx.fillStyle=d[2];ctx.lineWidth=1.4;ctx.beginPath();ctx.moveTo(origin.x,origin.y);ctx.lineTo(end.x,end.y);ctx.stroke();ctx.fillText(d[0],origin.x+ux*44,origin.y+uy*44)}ctx.textAlign="start";ctx.textBaseline="alphabetic"}
function render(){if(renderPending)return;renderPending=true;window.requestAnimationFrame(()=>{renderPending=false;drawScene()})}
function drawScene(){const r=canvas.getBoundingClientRect();frame={rect:r,matrix:rotationMatrix(state.orientation),scale:Math.min(r.width,r.height)/(data.bounds.radius*2.25)*state.zoom};ctx.clearRect(0,0,r.width,r.height);if(state.showBrain)cloud(data.brain,"#475569",.22,1.25);if(state.showTargets)targetRegions();if(state.showFocal)focalBoundary();lines();if(state.showHits)hits();if(state.showTargetHits)targetHits();if(state.showEndpoints&&data.endpointAnalysis)endpointMarkers();if(state.showTargets)targets();axes()}
updateClassSummary();
updateEndpointStatus();
resize();
if(document.fonts&&document.fonts.ready)document.fonts.ready.then(resize);
</script>
</body>
</html>')
}

ftq_qc_open_output <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  utils::browseURL(ftq_file_url(path))
  invisible(path)
}
