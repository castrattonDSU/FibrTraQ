#' Build tractogram-derived spatial context for reports
#'
#' Stream a TrackVis `.trk` file and build a discretized streamline-occupancy
#' density image on a reference voxel grid. This is a visualization and
#' quality-control context only: it is not an anatomical image and must not be
#' interpreted as brain, tissue, or lesion volume.
#'
#' Each streamline contributes at most one count to a voxel under the selected
#' stored-point or segment-interpolation sampling and discretization contract.
#'
#' @param tractogram TrackVis `.trk` file.
#' @param template_image NIfTI/MGH image defining the output voxel grid.
#' @param output_dir Directory where a `context/` subdirectory is written. If
#'   omitted, the template image directory is used.
#' @param prefix Output filename prefix.
#' @param track_space Track coordinate space. Use the same setting as the
#'   FibrTraQ run whenever possible. `"voxelmm"` uses TrackVis header voxel
#'   sizes and the TrackVis corner-to-center convention.
#' @param axis_flip Logical length-3 vector applied after discretization.
#' @param track_chunk_size Number of streamlines read per chunk.
#' @param density_discretization Discretization used for the visualization
#'   density. This does not affect FibrTraQ counting outputs.
#' @param streamline_voxelization_mode Whether density uses stored sample
#'   points only or interpolates segments at a maximum 0.5-voxel axis step.
#' @param ref_image Optional reference image for world/TrackVis conversion.
#'   Defaults to `template_image`.
#' @param overwrite Whether existing context outputs may be replaced.
#'
#' @return A list with `density_image`, `summary_csv`, `metadata_json`,
#'   `summary`, and `metadata`.
#' @export
ftq_tractogram_context <- function(tractogram,
                                   template_image,
                                   output_dir = NULL,
                                   prefix = "fibrtraq",
                                   track_space = c("auto", "voxelmm", "vox10", "voxel", "world", "trackvis"),
                                   axis_flip = c(FALSE, FALSE, FALSE),
                                   track_chunk_size = 50000L,
                                   density_discretization = c("round", "floor", "ceil"),
                                   streamline_voxelization_mode = c("sample_points", "segment_interpolation"),
                                   ref_image = NULL,
                                   overwrite = FALSE) {
  ftq_build_tractogram_context(
    tractogram = tractogram,
    template_image = template_image,
    output_dir = output_dir,
    prefix = prefix,
    track_space = track_space,
    axis_flip = axis_flip,
    track_chunk_size = track_chunk_size,
    density_discretization = density_discretization,
    streamline_voxelization_mode = streamline_voxelization_mode,
    ref_image = ref_image,
    overwrite = overwrite
  )
}

# Keep the standalone research export intact. Standard reports use maps/ and
# retain the returned summary and metadata in their existing run-metadata JSON.
ftq_build_tractogram_context <- function(tractogram,
                                         template_image,
                                         output_dir = NULL,
                                         prefix = "fibrtraq",
                                         track_space = c("auto", "voxelmm", "vox10", "voxel", "world", "trackvis"),
                                         axis_flip = c(FALSE, FALSE, FALSE),
                                         track_chunk_size = 50000L,
                                         density_discretization = c("round", "floor", "ceil"),
                                         streamline_voxelization_mode = c("sample_points", "segment_interpolation"),
                                         ref_image = NULL,
                                         overwrite = FALSE,
                                         density_subdir = "context",
                                         write_sidecars = TRUE) {
  tractogram <- normalizePath(path.expand(tractogram), mustWork = TRUE)
  template_image <- normalizePath(path.expand(template_image), mustWork = TRUE)
  if (!grepl("\\.trk$", tolower(tractogram))) {
    stop("ftq_tractogram_context() currently requires a TrackVis .trk file.")
  }
  track_space <- match.arg(track_space)
  density_discretization <- match.arg(density_discretization)
  streamline_voxelization_mode <- match.arg(streamline_voxelization_mode)
  axis_flip <- as.logical(axis_flip)[1:3]
  if (length(axis_flip) != 3L || anyNA(axis_flip)) {
    stop("axis_flip must be a logical vector of length 3.")
  }
  track_chunk_size <- as.integer(track_chunk_size)[1]
  if (is.na(track_chunk_size) || track_chunk_size <= 0L) {
    stop("track_chunk_size must be a positive integer.")
  }

  if (is.null(output_dir)) output_dir <- dirname(template_image)
  output_dir <- normalizePath(path.expand(output_dir), mustWork = FALSE)
  context_dir <- file.path(output_dir, density_subdir)
  dir.create(context_dir, recursive = TRUE, showWarnings = FALSE)

  paths <- list(
    density_image = file.path(context_dir, paste0(prefix, "_fibrtraq_tractogram-density.nii.gz")),
    summary_csv = if (isTRUE(write_sidecars)) file.path(context_dir, paste0(prefix, "_fibrtraq_tractogram-context-summary.csv")) else NULL,
    metadata_json = if (isTRUE(write_sidecars)) file.path(context_dir, paste0(prefix, "_fibrtraq_tractogram-context-metadata.json")) else NULL
  )
  ftq_check_overwrite(unlist(paths, use.names = FALSE), overwrite)

  template <- RNifti::readNifti(template_image, internal = FALSE)
  dims <- as.integer(dim(template)[1:3])
  if (length(dims) != 3L || anyNA(dims)) stop("template_image must be a 3D volume.")
  ref <- if (is.null(ref_image)) template else RNifti::readNifti(ref_image, internal = FALSE)
  density <- array(0, dim = dims)
  stream <- ftq_trk_stream(tractogram, chunk_size = track_chunk_size)
  track_affine <- ftq_trackvis_affine(stream$header)
  track_voxel_size <- ftq_trackvis_voxel_size(stream$header)
  selected_track_space <- ftq_qc_resolve_track_space(
    track_space,
    track_affine,
    ref,
    track_voxel_size = track_voxel_size
  )

  processed <- 0L
  mapped_errors <- 0L
  stored_points <- 0
  total_points <- 0
  out_of_bounds_points <- 0
  occupancy_events <- 0
  voxel_min <- rep(Inf, 3)
  voxel_max <- rep(-Inf, 3)

  map_track <- function(track) {
    ftq_qc_map_track(
      track = track,
      track_space = selected_track_space,
      ref_image = ref,
      track_affine = track_affine,
      track_voxel_size = track_voxel_size
    )
  }

  use_cpp_density <- is.function(get0("ftq_density_chunk_cpp", mode = "function"))
  raw_to_voxel <- if (use_cpp_density) {
    ftq_qc_raw_to_voxel_affine(
      track_space = selected_track_space,
      ref_image = ref,
      track_affine = track_affine,
      track_voxel_size = track_voxel_size
    )
  } else {
    NULL
  }

  ftq_trk_stream_iterate(
    stream,
    callback = function(chunk, start_index, end_index) {
      if (use_cpp_density) {
        res <- ftq_density_chunk_cpp(
          tracks = chunk,
          raw_to_voxel = raw_to_voxel,
          method = density_discretization,
          axis_flip = axis_flip,
          vol_dim = dims,
          streamline_voxelization_mode = streamline_voxelization_mode
        )
        density[] <<- density[] + as.integer(res$density)
        processed <<- processed + as.integer(res$processed)
        stored_points <<- stored_points + as.numeric(res$stored_points)
        total_points <<- total_points + as.numeric(res$voxelization_points)
        out_of_bounds_points <<- out_of_bounds_points + as.numeric(res$out_of_bounds_points)
        occupancy_events <<- occupancy_events + as.numeric(res$occupancy_events)
        chunk_min <- as.numeric(res$voxel_min)
        chunk_max <- as.numeric(res$voxel_max)
        if (all(chunk_max > 0) && all(chunk_min <= dims)) {
          voxel_min <<- pmin(voxel_min, chunk_min)
          voxel_max <<- pmax(voxel_max, chunk_max)
        }
        return(invisible(NULL))
      }
      for (track in chunk) {
        processed <<- processed + 1L
        mapped <- tryCatch(map_track(track), error = function(e) NULL)
        if (is.null(mapped)) {
          mapped_errors <<- mapped_errors + 1L
          next
        }
        stored_points <<- stored_points + nrow(mapped$voxel)
        if (identical(streamline_voxelization_mode, "segment_interpolation")) {
          mapped$voxel <- ftq_densify_voxel_path(mapped$voxel)
        }
        total_points <<- total_points + nrow(mapped$voxel)
        discrete <- switch(
          density_discretization,
          round = round(mapped$voxel),
          floor = floor(mapped$voxel),
          ceil = ceiling(mapped$voxel)
        )
        discrete <- ftq_qc_apply_flip(discrete, axis_flip, dims)
        key <- ftq_qc_encode_voxels(discrete, dims)
        out_of_bounds_points <<- out_of_bounds_points + sum(is.na(key))
        keep <- !is.na(key)
        if (!any(keep)) next
        kept_key <- key[keep]
        kept_vox <- discrete[keep, , drop = FALSE]
        first <- !duplicated(kept_key)
        unique_key <- kept_key[first]
        unique_vox <- kept_vox[first, , drop = FALSE]
        density[unique_key] <<- density[unique_key] + 1
        occupancy_events <<- occupancy_events + length(unique_key)
        voxel_min <<- pmin(voxel_min, apply(unique_vox, 2, min))
        voxel_max <<- pmax(voxel_max, apply(unique_vox, 2, max))
      }
    },
    chunk_size = track_chunk_size
  )

  voxel_volume <- ftq_voxel_volume_mm3(template)
  occupied_voxels <- sum(density > 0)
  bbox_available <- all(is.finite(voxel_min)) && all(is.finite(voxel_max))
  bbox_world <- if (bbox_available) ftq_bbox_world(template, voxel_min, voxel_max) else rep(NA_real_, 6L)
  declared <- as.integer(stream$header$n_count)
  if (is.na(declared) || declared <= 0L) declared <- NA_integer_

  summary <- data.frame(
    processed_streamlines = as.numeric(processed),
    declared_streamlines = as.numeric(declared),
    nonzero_tractogram_occupied_voxels = as.numeric(occupied_voxels),
    tractogram_occupancy_volume_mm3 = as.numeric(occupied_voxels * voxel_volume),
    total_streamline_voxel_occupancy_events = as.numeric(occupancy_events),
    max_density = if (length(density)) max(density) else 0,
    total_sample_points = as.numeric(stored_points),
    total_voxelization_points = as.numeric(total_points),
    out_of_bounds_sample_points = as.numeric(out_of_bounds_points),
    mapping_error_streamlines = as.numeric(mapped_errors),
    voxel_bbox_x_min = if (bbox_available) voxel_min[1] else NA_real_,
    voxel_bbox_x_max = if (bbox_available) voxel_max[1] else NA_real_,
    voxel_bbox_y_min = if (bbox_available) voxel_min[2] else NA_real_,
    voxel_bbox_y_max = if (bbox_available) voxel_max[2] else NA_real_,
    voxel_bbox_z_min = if (bbox_available) voxel_min[3] else NA_real_,
    voxel_bbox_z_max = if (bbox_available) voxel_max[3] else NA_real_,
    world_bbox_x_min = bbox_world[1],
    world_bbox_x_max = bbox_world[2],
    world_bbox_y_min = bbox_world[3],
    world_bbox_y_max = bbox_world[4],
    world_bbox_z_min = bbox_world[5],
    world_bbox_z_max = bbox_world[6],
    note = "Tractogram-derived density context is not anatomical brain volume.",
    stringsAsFactors = FALSE
  )
  metadata <- list(
    source = "tractogram-derived density context",
    warning = "This is not an anatomical image and is not brain, tissue, or lesion volume.",
    density_unit = "streamline occupancy per voxel",
    duplicate_point_policy = "Each streamline contributes at most once per voxel.",
    tractogram = tractogram,
    template_image = template_image,
    selected_track_space = selected_track_space,
    track_voxel_size_mm = as.numeric(track_voxel_size),
    axis_flip = as.logical(axis_flip),
    density_discretization = density_discretization,
    streamline_voxelization_mode = streamline_voxelization_mode,
    track_chunk_size = track_chunk_size,
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    summary = summary
  )

  RNifti::writeNifti(density, paths$density_image, template = template, datatype = "float")
  if (isTRUE(write_sidecars)) {
    utils::write.csv(summary, paths$summary_csv, row.names = FALSE)
    jsonlite::write_json(metadata, paths$metadata_json, auto_unbox = TRUE, pretty = TRUE, null = "null")
  }

  c(paths, list(summary = summary, metadata = metadata))
}

ftq_voxel_volume_mm3 <- function(image) {
  pts <- rbind(
    c(1, 1, 1),
    c(2, 1, 1),
    c(1, 2, 1),
    c(1, 1, 2)
  )
  world <- tryCatch(ftq_voxel_to_world_mm(pts, image), error = function(e) NULL)
  if (is.null(world) || !identical(dim(world), c(4L, 3L))) return(NA_real_)
  abs(det(rbind(
    world[2, ] - world[1, ],
    world[3, ] - world[1, ],
    world[4, ] - world[1, ]
  )))
}

ftq_bbox_world <- function(image, voxel_min, voxel_max) {
  corners <- expand.grid(
    x = c(voxel_min[1], voxel_max[1]),
    y = c(voxel_min[2], voxel_max[2]),
    z = c(voxel_min[3], voxel_max[3])
  )
  world <- tryCatch(RNifti::voxelToWorld(as.matrix(corners), image, simple = FALSE), error = function(e) NULL)
  if (is.null(world)) return(rep(NA_real_, 6L))
  c(
    min(world[, 1]), max(world[, 1]),
    min(world[, 2]), max(world[, 2]),
    min(world[, 3]), max(world[, 3])
  )
}

ftq_context_density_image <- function(context) {
  if (is.null(context)) return(NULL)
  if (is.character(context) && length(context) == 1L && nzchar(context)) return(context)
  if (is.list(context)) {
    path <- ftq_value(context$density_image, context$path, context$image, NULL)
    if (!is.null(path) && length(path) == 1L && !is.na(path) && nzchar(path)) return(path)
  }
  NULL
}
