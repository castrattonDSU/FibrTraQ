#' Create a connectivity table from tractography data (FibrTraQ)
#'
#' Orchestrates the workflow:
#' 1) Extract focal, local, and cNode ROI coordinates (threshold-aware).
#' 2) Propose a track-to-voxel mapping using a small sample of tracks when
#'    `track_space="auto"` and/or `axis_flip` is not explicitly provided.
#'    A reflected candidate requires explicit, independently verified axis flips.
#' 3) Identify all streamlines intersecting the focal ROI.
#' 4) Count focal-voxel/streamline associations and target-ROI connections
#'    using the selected target-assignment method.
#'
#' Coordinate handling notes:
#' - Native MRtrix TCK objects use RAS+ world millimetres and the declared NIfTI
#'   transform. No overlap-based space or reflection search is performed. The
#'   reference must declare spatial units and a valid, unambiguous qform/sform.
#' - Some tractography encodes voxel indices at a 0.1 scale ("vox10"): multiply
#'   coordinates by 10 before discretization.
#' - Some datasets require axis flips (e.g., flip_z) to align coordinate conventions
#'   between tractography and ROI arrays. When `axis_flip` is not provided, an
#'   automatic search over flips is performed for the sampled tracks, but a
#'   reflected candidate stops the analysis. Sample overlap cannot establish
#'   anatomical alignment; verify the mapping before setting flips explicitly.
#'
#' @param fNode 3-D array for the focal ROI and authoritative voxel grid (for
#'   example, a NIfTI image read by RNifti). A declared valid qform/sform
#'   enables automatic affine as well as dimension validation.
#' @param localROI A 3-D mask on the `fNode` grid or an n x 3 matrix/data frame
#'   of finite, integer-compatible, in-bounds voxel indices. If \code{NULL}, the
#'   focal ROI is used as the local mask. When both masks declare a valid
#'   qform/sform, their unit-normalized affines must match. Coordinate matrices
#'   carry no affine provenance. The \code{local} output is a non-exclusive
#'   flag count.
#' @param cNode_list List of target ROIs, each supplied as a 3-D mask on the
#'   `fNode` grid or an n x 3 coordinate matrix/data frame satisfying the same
#'   coordinate checks as `localROI`.
#' @param cNode_labels Character vector of names for cNodes (same length/order as list).
#' @param all_tracks Tractography object as read by freesurferformats::read.dti.trk()
#'   or freesurferformats::read.dti.tck(), an internal streaming `.trk` or `.tck`
#'   object, or a list of tracks (each track either
#'   list(num_points, coords) or an n x 3 matrix).
#' @param focal_threshold "nonzero" (default) or numeric threshold for ROI voxels.
#' @param discretization Discretization method(s) used by iso_tracks() when mapping
#'   continuous coordinates to voxel indices. May be "any" (default) or one/more of
#'   c("round","floor","ceil").
#' @param focal_assignment_mode How streamlines are assigned to focal voxels.
#'   `"all_focal_voxels"` counts each distinct focal voxel occupied by one or
#'   more discretized stored streamline sample points once and is invariant to
#'   streamline point-order reversal.
#'   `"first_intersection"` preserves the legacy behavior of assigning each
#'   streamline only to the first focal voxel encountered in stored point order.
#' @param streamline_voxelization_mode How streamline coordinates are converted
#'   to voxel hits. `"sample_points"` uses only stored streamline points.
#'   `"segment_interpolation"` linearly inserts points so no interpolated step
#'   exceeds 0.5 voxel along any axis before discretization, reducing missed
#'   voxel intersections when tractograms store sparse points.
#' @param reverse_streamlines Logical; if TRUE, reverse each streamline's point
#'   order in memory before counting. This is intended for orientation/QC
#'   sensitivity checks and does not write a reversed tractogram.
#' @param verbose Logical; diagnostics.
#' @param track_space Track coordinate space: one of c("auto","voxelmm",
#'   "vox10","voxel","world","trackvis"). Default "auto". `"voxelmm"`
#'   converts TrackVis voxel-mm coordinates to 1-based voxel-center indices
#'   using the header voxel sizes. `"trackvis"` applies a TrackVis header
#'   affine before converting world coordinates to ROI voxel indices.
#'   `"world"` and TrackVis affine outputs are interpreted as millimetres and
#'   converted to the reference NIfTI's declared native spatial unit.
#'   Recognized TCK objects accept only `"auto"` or `"world"`, both of which use
#'   the format's fixed RAS+ world-millimetre convention.
#' @param ref_image Reference image used when `track_space` is `"world"` or
#'   `"trackvis"`, considered by `"auto"`, and used for `"local_search"`
#'   geometry. If NULL, falls back to `fNode`. An explicitly used reference
#'   must have the focal dimensions and, when both objects declare a valid
#'   qform/sform, the same unit-normalized affine/grid.
#'   TCK requires a declared qform or sform and known spatial units; conflicting
#'   coded transforms are rejected rather than selected by overlap.
#' @param auto_n Number of tracks to sample when selecting mapping/flips. Default 200.
#' @param axis_flip Logical length-3 vector c(flip_x, flip_y, flip_z). If not supplied,
#'   flips are searched automatically, but analysis stops if the best candidate
#'   would reflect any axis. Verify anatomical alignment and explicitly set
#'   `track_space` and `axis_flip` before continuing. If supplied, the provided
#'   flips are used; they remain the caller's alignment decision.
#'   TCK objects accept only no flips; apply intended registration upstream.
#' @param track_affine Optional 4x4 affine that maps raw track coordinates to
#'   world coordinates. If omitted and `all_tracks` has a TrackVis
#'   `header$vox2ras`, that header affine is used.
#'   An explicit value is not accepted for recognized TCK objects.
#' @param track_voxel_size Optional numeric length-3 TrackVis header voxel sizes
#'   in mm. Required for `track_space = "voxelmm"`; inferred from
#'   `all_tracks$header` when available.
#'   An explicit value is not accepted for recognized TCK objects.
#' @param backend Counting backend. `"r"` uses the reference R implementation.
#'   `"cpp"` uses the compiled backend when available. `"auto"` uses the
#'   compiled backend when installed and falls back to the reference R backend.
#' @param progress Logical; emit chunk-level progress messages during the main
#'   tractogram pass. Defaults to the value of `verbose`.
#' @param progress_every Number of streamlines between progress messages.
#' @param profile Logical; attach stage timing information to the returned table.
#' @param endpoint_analysis Additive endpoint-versus-intersection analysis.
#'   `"exact_and_dilated1"` reports exact ROI endpoint membership and a fixed
#'   one-voxel 26-neighbor sensitivity analysis. `"none"` disables the
#'   additional analysis without changing existing FibrTraQ counts.
#' @param endpoint_context_list Optional list of finite binary 3-D anatomical
#'   context masks on the focal-image grid. Contexts audit endpoint containment
#'   and never alter streamline counting.
#' @param endpoint_context_labels Labels for `endpoint_context_list`. Defaults
#'   to the list names and must be non-empty and unique when contexts are used.
#' @param assignment_mode How a focal-positive streamline is connected
#'   to target ROIs. `"full_streamline"` is the backward-compatible
#'   all-voxel method: any selected path voxel in a target records a
#'   connection. `"endpoint_voxels"` requires the two endpoint voxels to
#'   resolve unambiguously to the focal ROI and one target ROI under the first
#'   focal-positive method in the ordered `discretization` candidates.
#'   `"local_search"` independently assigns each endpoint to the nearest
#'   labelled focal/target voxel centre within `search_radius_mm`, then
#'   requires an unambiguous focal-target endpoint pair. Equal-nearest distinct
#'   labels and voxels carrying overlapping node labels are ambiguous and left
#'   unassigned; equal-nearest voxels carrying the same label remain
#'   assignable. The aliases `"all_voxel"`, `"all_voxels"`, `"end_voxel"`,
#'   `"end_voxels"`, and `"endpoint_voxel"` are accepted. Target assignment
#'   does not change focal-path qualification or focal-voxel allocation.
#' @param search_radius_mm Positive endpoint search radius in physical
#'   millimetres for `assignment_mode = "local_search"`. The default
#'   2 mm matches the radius reported by Yeh et al. (2019); it is not asserted
#'   to be universally optimal or an exact reproduction of unspecified RC3
#'   implementation details. Distances are measured from the precise endpoint
#'   to labelled voxel centres in scanner space. NIfTI metre and micrometre
#'   spatial units are converted to millimetres. Missing spatial units are
#'   interpreted as millimetres with a warning. With a plain-array `fNode`,
#'   dimensions are checked but the caller must establish the reference
#'   affine/grid because the focal array carries no spatial metadata.
#'
#' @return Data frame: one row per focal voxel with focal-voxel/streamline
#'   association counts (\code{tracks}), the non-exclusive local flag count,
#'   non-local associations with no target-ROI connection under the selected
#'   target-assignment method
#'   (\code{leave_no_target}), one column per target ROI containing
#'   focal-voxel/streamline/target associations, and a \code{none} flag for
#'   focal voxels with zero associations. ROI-level unique-streamline and
#'   mutually exclusive fate summaries are attached as
#'   \code{fibrtraq_roi_summary}.
#'
#' @examples
#' # res <- FibrTraQ(fnode, NULL, cnodes, names(cnodes), trk)
#' @export
FibrTraQ <- function(fNode,
                     localROI = NULL,
                     cNode_list,
                     cNode_labels,
                     all_tracks,
                     focal_threshold = "nonzero",
	                     discretization = "any",
	                     focal_assignment_mode = c("all_focal_voxels", "first_intersection"),
	                     streamline_voxelization_mode = c("sample_points", "segment_interpolation"),
	                     reverse_streamlines = FALSE,
                     verbose = FALSE,
                     track_space = c("auto","voxelmm","vox10","voxel","world","trackvis"),
                     ref_image = NULL,
                     auto_n = 200L,
                     axis_flip = c(FALSE,FALSE,FALSE),
                     track_affine = NULL,
                     track_voxel_size = NULL,
                     backend = c("r", "cpp", "auto"),
                     progress = verbose,
                     progress_every = 100000L,
                     profile = FALSE,
                     endpoint_analysis = c("exact_and_dilated1", "none"),
                     endpoint_context_list = NULL,
                     endpoint_context_labels = names(endpoint_context_list),
                     assignment_mode = c(
                       "full_streamline",
                       "endpoint_voxels",
                       "local_search"
                     ),
                     search_radius_mm = 2) {

  track_space <- match.arg(track_space)
	  backend <- match.arg(backend)
	  focal_assignment_mode <- match.arg(focal_assignment_mode)
	  streamline_voxelization_mode <- match.arg(streamline_voxelization_mode)
	  endpoint_analysis <- match.arg(endpoint_analysis)
	  endpoint_enabled <- identical(endpoint_analysis, "exact_and_dilated1")
	  target_intersection_mode <- ftq_normalize_target_intersection_mode(
	    assignment_mode
	  )
	  local_search_radius_mm <- ftq_validate_local_search_radius(
	    search_radius_mm
	  )
  reverse_streamlines <- isTRUE(reverse_streamlines)
  auto_n <- as.integer(auto_n)
  progress_every <- as.integer(progress_every)[1]
  if (is.na(progress_every) || progress_every <= 0L) {
    stop("progress_every must be a positive integer.")
  }
  progress <- isTRUE(progress)
  profile <- isTRUE(profile)
  backend_used <- ftq_resolve_backend(backend)
  if (identical(backend_used, "cpp") &&
      !ftq_compiled_backend_available()) {
    stop("backend = 'cpp' requires the package compiled-code bindings to be available.")
  }

  ftq_elapsed <- function() as.numeric(proc.time()[["elapsed"]])
  profile_rows <- list()
  profile_stage_start <- ftq_elapsed()
  add_profile_stage <- function(stage) {
    if (!isTRUE(profile)) return(invisible(NULL))
    now <- ftq_elapsed()
    profile_rows[[length(profile_rows) + 1L]] <<- data.frame(
      stage = stage,
      seconds = now - profile_stage_start,
      backend = backend_used,
      stringsAsFactors = FALSE
    )
    profile_stage_start <<- now
    invisible(NULL)
  }

  axis_flip_user <- !missing(axis_flip)
  if (!is.logical(axis_flip) || length(axis_flip) != 3L || anyNA(axis_flip)) {
    stop("axis_flip must be a logical vector of length 3 with no missing values: c(flip_x, flip_y, flip_z).")
  }

  input_format <- ftq_track_format(all_tracks)
  tck_geometry <- NULL
  if (identical(input_format, "tck")) {
    if (!is.null(track_affine) || !is.null(track_voxel_size)) {
      stop("TCK coordinates are already in RAS+ world millimetres; ",
           "track_affine and track_voxel_size must be NULL.", call. = FALSE)
    }
    settings <- ftq_tck_coordinate_settings(
      track_space, if (axis_flip_user) axis_flip else NULL,
      if (is.null(ref_image)) fNode else ref_image
    )
    track_space <- settings$track_space
    axis_flip <- settings$axis_flip
    axis_flip_user <- TRUE
    ref_image <- settings$reference
    tck_geometry <- settings$metadata
    # Validate every supplied image's geometry before any coordinate mapping.
    # Plain coordinate arrays retain the explicit reference-grid contract.
    images <- c(list(fNode = fNode), cNode_list,
                if (!is.null(localROI)) list(localROI) else list(),
                endpoint_context_list)
    for (image in images) {
      if (inherits(image, "niftiImage")) {
        ftq_tck_validate_reference(image, "ROI image")
        if (!identical(as.integer(dim(image)), as.integer(dim(ref_image))) ||
            !ftq_local_search_affine_equal(image, ref_image)) {
          stop("TCK ROI images must use the same dimensions and unit-normalized ",
               "affine/grid as the reference image.", call. = FALSE)
        }
      }
    }
  }

  local_roi_defaulted_to_focal <- is.null(localROI)
  if (is.null(localROI)) localROI <- fNode
  if (length(cNode_list) != length(cNode_labels)) {
    stop("cNode_list and cNode_labels length mismatch")
  }
  if (anyDuplicated(cNode_labels)) {
    stop("cNode_labels must be unique.")
  }
  reserved_cols <- c("fNode_Coords", "tracks", "local", "leave_no_target",
                     "none", "cNode_Coords")
  if (any(cNode_labels %in% reserved_cols)) {
    stop("cNode_labels cannot use reserved output column names: ",
         paste(intersect(cNode_labels, reserved_cols), collapse = ", "))
  }

  # fNode must be a 3-D array so we know dims for key encoding and flips
  vol_dim <- dim(fNode)
  if (is.null(vol_dim) || length(vol_dim) != 3L) {
    stop("fNode must be a 3-D array (NIfTI volume).")
  }
  vol_dim <- as.integer(vol_dim[1:3])
  if (anyNA(vol_dim) || any(vol_dim <= 0L)) {
    stop("fNode must have three positive voxel-grid dimensions.", call. = FALSE)
  }

  # Every ROI used by the centralized counting engine must describe the same
  # voxel grid. Coordinate matrices are supported for local/target ROIs, but
  # they carry no affine provenance; validate their coordinates directly and
  # reserve affine comparisons for image objects with real spatial metadata.
  localROI <- ftq_validate_direct_roi_grid(
    localROI,
    fNode = fNode,
    vol_dim = vol_dim,
    input_name = "localROI"
  )
  if (length(cNode_list)) {
    cNode_list <- lapply(seq_along(cNode_list), function(i) {
      ftq_validate_direct_roi_grid(
        cNode_list[[i]],
        fNode = fNode,
        vol_dim = vol_dim,
        input_name = paste0("cNode_list[[", i, "]] ('", cNode_labels[[i]], "')")
      )
    })
  }

  ref_image_is_used <- !is.null(ref_image) && (
    track_space %in% c("world", "trackvis", "auto") ||
      identical(target_intersection_mode, "local_search")
  )
  if (isTRUE(ref_image_is_used)) {
    ftq_validate_direct_reference_grid(
      fNode = fNode,
      ref_image = ref_image,
      vol_dim = vol_dim
    )
  }

  if (!endpoint_enabled && length(endpoint_context_list)) {
    stop("endpoint_context_list requires endpoint_analysis = 'exact_and_dilated1'.")
  }
  endpoint_context_info <- ftq_endpoint_prepare_contexts(
    endpoint_context_list,
    endpoint_context_labels,
    vol_dim,
    reference = fNode
  )

  # 1) ROI coords (threshold-aware)
  fNode_coords   <- iso_coords(fNode,   threshold = focal_threshold)
  local_coords   <- iso_coords(localROI, threshold = focal_threshold)
  cNode_coords_list <- lapply(cNode_list, iso_coords, threshold = focal_threshold)

  # Track normalization: accept a supported in-memory object/list or stream.
  track_list <- NULL
  track_stream <- NULL
  if (ftq_is_track_stream(all_tracks)) {
    track_stream <- all_tracks
    if (!identical(input_format, "tck") && is.null(track_affine)) {
      track_affine <- ftq_trackvis_affine(track_stream$header)
    }
    if (!identical(input_format, "tck") && is.null(track_voxel_size)) {
      track_voxel_size <- ftq_trackvis_voxel_size(track_stream$header)
    }
  } else if (is.list(all_tracks) && !is.null(all_tracks$tracks) && is.list(all_tracks$tracks)) {
    if (!identical(input_format, "tck") && is.null(track_affine)) {
      track_affine <- ftq_trackvis_affine(all_tracks$header)
    }
    if (!identical(input_format, "tck") && is.null(track_voxel_size)) {
      track_voxel_size <- ftq_trackvis_voxel_size(all_tracks$header)
    }
    track_list <- all_tracks$tracks
  } else if (is.list(all_tracks)) {
    track_list <- all_tracks
  } else {
    stop("all_tracks must be a TRK/TCK object with $tracks, a supported stream object, or a list of tracks.")
  }
  sample_tracks <- if (!is.null(track_stream)) {
    ftq_track_stream_sample(track_stream, auto_n)
  } else {
    if (!length(track_list)) stop("all_tracks contains no tracks.")
    track_list[seq_len(min(auto_n, length(track_list)))]
  }
  if (!length(sample_tracks)) stop("all_tracks contains no tracks.")

  n_tracks <- if (!is.null(track_stream)) {
    ftq_track_count(track_stream)
  } else {
    length(track_list)
  }
  iterate_tracks <- function(callback) {
    if (!is.null(track_stream)) {
      ftq_track_stream_iterate(track_stream, callback = callback, chunk_size = track_stream$chunk_size)
    } else {
      callback(track_list, 1L, length(track_list))
      invisible(length(track_list))
    }
  }

  # Discretization methods for detection
  methods <- if (identical(discretization, "any") || any(discretization == "any")) {
    c("round","floor","ceil")
  } else {
    unique(match.arg(discretization, c("round","floor","ceil"), several.ok = TRUE))
  }

  # Reference image for world conversion (legacy fallback to fNode)
  ref_use <- ref_image
  if (is.null(ref_use)) ref_use <- fNode
  if (!is.null(track_affine)) {
    track_affine <- as.matrix(track_affine)
    if (!identical(dim(track_affine), c(4L, 4L)) || anyNA(track_affine)) {
      stop("track_affine must be a numeric 4x4 matrix.")
    }
    storage.mode(track_affine) <- "numeric"
  }
  if (!is.null(track_voxel_size)) {
    track_voxel_size <- ftq_validate_track_voxel_size(track_voxel_size)
  }
  if (identical(track_space, "voxelmm") && is.null(track_voxel_size)) {
    stop("track_space = 'voxelmm' requires TrackVis header voxel sizes or track_voxel_size.")
  }
  local_search_geometry <- if (identical(
    target_intersection_mode,
    "local_search"
  )) {
    ftq_local_search_geometry(fNode, ref_image = ref_image)
  } else {
    list(
      voxel_to_world = diag(3L),
      voxel_radius_per_mm = rep(1, 3L),
      voxel_size_mm = rep(NA_real_, 3L)
    )
  }
  if (identical(target_intersection_mode, "local_search") &&
      !isTRUE(local_search_geometry$spatial_units_known)) {
    warning(
      "The local-search reference image does not declare NIfTI spatial ",
      "units; affine coordinates are being interpreted as millimetres. ",
      "Set a valid spatial unit or verify that this assumption is correct.",
      call. = FALSE
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

  # Encode integer xyz -> linear index (NA if out-of-bounds)
  ftq_encode <- function(m) {
    mm <- as.matrix(m)[,1:3, drop=FALSE]
    storage.mode(mm) <- "integer"
    x <- mm[,1]; y <- mm[,2]; z <- mm[,3]
    inb <- x>=1L & x<=vol_dim[1] & y>=1L & y<=vol_dim[2] & z>=1L & z<=vol_dim[3]
    k <- x + (y-1L)*vol_dim[1] + (z-1L)*(vol_dim[1]*vol_dim[2])
    k[!inb] <- NA_integer_
    k
  }

  # Apply flips to integer voxel coords
  ftq_apply_flip <- function(m, flip) {
    mm <- m
    if (flip[1]) mm[,1] <- vol_dim[1] - mm[,1] + 1L
    if (flip[2]) mm[,2] <- vol_dim[2] - mm[,2] + 1L
    if (flip[3]) mm[,3] <- vol_dim[3] - mm[,3] + 1L
    storage.mode(mm) <- "integer"
    mm
  }

  # Build LUTs: focal-only for mapping selection and all-ROI overlap for
  # diagnostics. Target/local hits are reported but not used to choose mappings.
  lut_f <- rep(FALSE, prod(vol_dim))
  kf <- ftq_encode(fNode_coords)
  lut_f[kf[!is.na(kf)]] <- TRUE

  lut_u <- rep(FALSE, prod(vol_dim))
  set_lut <- function(coords) {
    if (is.null(coords) || !nrow(coords)) return(invisible(NULL))
    kk <- ftq_encode(coords)
    lut_u[kk[!is.na(kk)]] <<- TRUE
    invisible(NULL)
  }
  set_lut(fNode_coords)
  set_lut(local_coords)
  if (length(cNode_coords_list)) {
    for (cc in cNode_coords_list) set_lut(cc)
  }

  # Candidate track spaces for detection
  spaces <- if (track_space == "auto") {
    c(
      if (!is.null(track_voxel_size)) "voxelmm",
      "vox10", "voxel", "world",
      if (!is.null(track_affine)) "trackvis"
    )
  } else {
    track_space
  }

  # Candidate flips: if user supplied axis_flip, use only that; otherwise search all 8
  flip_list <- if (axis_flip_user) {
    list(as.logical(axis_flip))
  } else {
    list(
      c(FALSE,FALSE,FALSE),
      c(TRUE, FALSE,FALSE),
      c(FALSE,TRUE, FALSE),
      c(TRUE, TRUE, FALSE),
      c(FALSE,FALSE,TRUE),
      c(TRUE, FALSE,TRUE),
      c(FALSE,TRUE, TRUE),
      c(TRUE, TRUE, TRUE)
    )
  }

  # Convert a track coords matrix by space (returns numeric n x 3)
  ftq_map_coords <- function(tk, sp) {
    trk <- as_track(tk)
    c0 <- as.matrix(trk$coords)[,1:3, drop=FALSE]
    if (sp == "voxelmm") {
      return(ftq_trackvis_voxelmm_to_voxel(c0, track_voxel_size))
    }
    if (sp == "vox10") return(c0 * 10)
    if (sp == "voxel") return(c0)
    if (sp == "world") {
      return(ftq_world_mm_to_voxel(c0, ref_use))
    }
    if (sp == "trackvis") {
      if (is.null(track_affine)) {
        stop("track_space='trackvis' requires a TrackVis header or track_affine.")
      }
      pts_mm <- cbind(c0, 1) %*% t(track_affine)
      return(ftq_world_mm_to_voxel(
        pts_mm[, 1:3, drop = FALSE],
        ref_use
      ))
    }
    stop("Unknown track_space: ", sp)
  }

  ftq_reverse_track <- function(tk) {
    if (is.matrix(tk)) {
      return(tk[nrow(tk):1L, , drop = FALSE])
    }
    if (is.list(tk) && !is.null(tk$coords)) {
      tk$coords <- tk$coords[nrow(tk$coords):1L, , drop = FALSE]
      return(tk)
    }
    tk
  }

  # Build a linear raw-coordinate to voxel-coordinate transform for the C++
  # chunk backend. The transform is derived through the same mapper used by the
  # R backend, including unit-normalized world-to-voxel conversion, so TrackVis
  # and world-space millimetres keep the same 1-based voxel convention.
  ftq_raw_to_voxel_affine <- function(sp) {
    basis <- rbind(
      c(0, 0, 0),
      c(1, 0, 0),
      c(0, 1, 0),
      c(0, 0, 1)
    )
    mapped <- ftq_map_coords(basis, sp)
    origin <- mapped[1, ]
    aff <- matrix(0, nrow = 4L, ncol = 4L)
    aff[1:3, 1] <- mapped[2, ] - origin
    aff[1:3, 2] <- mapped[3, ] - origin
    aff[1:3, 3] <- mapped[4, ] - origin
    aff[1:3, 4] <- origin
    aff[4, 4] <- 1
    aff
  }

  # Score (space, flip) combos on a sample of tracks
  ftq_detect_mapping <- function() {
    n <- length(sample_tracks)
    idx <- seq_len(n)

    scores <- data.frame(
      track_space = character(0),
      flip_x = logical(0),
      flip_y = logical(0),
      flip_z = logical(0),
      hits_fnode = integer(0),
      hits_union = integer(0),
      in_bounds = numeric(0),
      stringsAsFactors = FALSE
    )

    for (sp in spaces) {
      for (flip in flip_list) {
        hf <- 0L
        hu <- 0L
        inb_num <- 0L
        inb_den <- 0L

        for (i in idx) {
          tk <- sample_tracks[[i]]
          c1 <- tryCatch(ftq_map_coords(tk, sp), error = function(e) NULL)
          if (is.null(c1)) next

          # in-bounds rate from round() only (stable metric)
          cr <- round(c1)
          storage.mode(cr) <- "integer"
          cr <- ftq_apply_flip(cr, flip)
          kr <- ftq_encode(cr)
          inb_num <- inb_num + sum(!is.na(kr))
          inb_den <- inb_den + length(kr)

          # hit tests (any discretization)
          hit_f <- FALSE
          hit_u <- FALSE
          for (mth in methods) {
            ci <- switch(mth,
              round = round(c1),
              floor = floor(c1),
              ceil  = ceiling(c1)
            )
            storage.mode(ci) <- "integer"
            ci <- ftq_apply_flip(ci, flip)
            kk <- ftq_encode(ci)

            if (!hit_f && any(!is.na(kk) & lut_f[kk])) hit_f <- TRUE
            if (!hit_u && any(!is.na(kk) & lut_u[kk])) hit_u <- TRUE
            if (hit_f && hit_u) break
          }

          if (hit_f) hf <- hf + 1L
          if (hit_u) hu <- hu + 1L
        }

        scores <- rbind(
          scores,
          data.frame(
            track_space = sp,
            flip_x = flip[1],
            flip_y = flip[2],
            flip_z = flip[3],
            hits_fnode = hf,
            hits_union = hu,
            in_bounds = if (inb_den > 0L) inb_num / inb_den else NA_real_,
            stringsAsFactors = FALSE
          )
        )
      }
    }

    # Choose best: maximize focal hits, then in-bounds rate, then prefer fewer
    # flips, then prefer the header-defined TrackVis voxel-mm convention before
    # legacy or less specific mappings. hits_union is kept only for audit output;
    # target/local ROI hits must not select mappings.
    pref_space <- c("voxelmm", "vox10", "voxel", "world", "trackvis")
    flip_n <- (scores$flip_x + scores$flip_y + scores$flip_z)
    space_rank <- match(scores$track_space, pref_space)

    o <- order(
      -scores$hits_fnode,
      -ifelse(is.na(scores$in_bounds), -Inf, scores$in_bounds),
      flip_n,
      space_rank
    )

    best <- scores[o[1], , drop = FALSE]

    if (isTRUE(verbose)) {
      message("Auto mapping candidate (n = ", n, "):")
      print(best, row.names = FALSE)
    }

    list(
      track_space = best$track_space,
      axis_flip = c(best$flip_x, best$flip_y, best$flip_z),
      scores = scores,
      selection_rule = "max_hits_fnode_then_in_bounds_then_fewest_flips_then_space_preference"
    )
  }

  # Decide final mapping + flips
  space_used <- track_space
  flip_used <- axis_flip
  mapping_scores <- NULL
  mapping_selection_rule <- if (identical(input_format, "tck")) {
    "tck_world_ras_mm_from_nifti_affine"
  } else "pinned"

  if (length(spaces) == 1L && axis_flip_user) {
    # Fully pinned: no detection
    space_used <- spaces
    flip_used <- axis_flip
  } else {
    det <- ftq_detect_mapping()
    space_used <- det$track_space
    flip_used <- det$axis_flip
    mapping_scores <- det$scores
    mapping_selection_rule <- det$selection_rule
  }

  # Focal overlap is a numerical score, not evidence of anatomical orientation.
  # Never turn a reflected candidate into a counting transform without an
  # explicit caller decision, and do not silently substitute another candidate.
  if (!axis_flip_user && any(flip_used)) {
    stop(
      "Automatic mapping would reflect one or more axes ",
      "(track_space = '", space_used, "', axis_flip = c(",
      paste(as.logical(flip_used), collapse = ", "), ")). ",
      "Analysis stopped before streamline counting. ",
      "Sampled focal-mask overlap cannot establish anatomical alignment. ",
      "Verify the tractogram and mask orientation against trusted geometry ",
      "or anatomy, then rerun FibrTraQ() or ftq_run_files() with track_space ",
      "and axis_flip explicitly set to the validated mapping. ",
      "For file workflows, use a new output directory.",
      call. = FALSE
    )
  }

  if (isTRUE(verbose)) {
    message("Using track_space = ", space_used,
            " ; axis_flip = c(", paste(as.logical(flip_used), collapse=","), ")")
  }
  add_profile_stage("mapping_detection")

	  # 3) Count exactly while streaming over tracks using the selected stored-point
	  # or segment-interpolation voxelization contract.
  n_fnode <- nrow(fNode_coords)
  f_labels <- if (n_fnode) {
    apply(fNode_coords, 1L, function(v) paste(v[1], v[2], v[3], sep = "_"))
  } else {
    character(0L)
  }
  tracks_counts <- integer(n_fnode)
  local_counts <- integer(n_fnode)
  leave_no_target_counts <- integer(n_fnode)
  target_counts <- matrix(0L, nrow = n_fnode, ncol = length(cNode_labels))
  colnames(target_counts) <- cNode_labels

  f_index <- rep.int(NA_integer_, prod(vol_dim))
  kf_valid <- kf[!is.na(kf)]
  if (length(kf_valid)) {
    row_idx <- which(!is.na(kf))
    keep_first <- !duplicated(kf_valid)
    f_index[kf_valid[keep_first]] <- row_idx[keep_first]
  }

  lut_local <- rep(FALSE, prod(vol_dim))
  kl <- ftq_encode(local_coords)
  lut_local[kl[!is.na(kl)]] <- TRUE

  target_lut <- matrix(FALSE, nrow = prod(vol_dim), ncol = length(cNode_labels))
  colnames(target_lut) <- cNode_labels
  if (length(cNode_coords_list)) {
    for (j in seq_along(cNode_coords_list)) {
      kj <- ftq_encode(cNode_coords_list[[j]])
      target_lut[kj[!is.na(kj)], j] <- TRUE
    }
  }
  lut_f_dilated <- if (endpoint_enabled) {
    ftq_endpoint_dilate_lut(lut_f, vol_dim)
  } else {
    rep(FALSE, length(lut_f))
  }
  target_lut_dilated <- matrix(
    FALSE,
    nrow = nrow(target_lut),
    ncol = ncol(target_lut),
    dimnames = dimnames(target_lut)
  )
  if (endpoint_enabled && ncol(target_lut)) {
    for (target in seq_len(ncol(target_lut))) {
      target_lut_dilated[, target] <-
        ftq_endpoint_dilate_lut(target_lut[, target], vol_dim)
    }
  }
  endpoint_accumulator <- ftq_endpoint_accumulator(
    length(cNode_labels),
    ncol(endpoint_context_info$lut)
  )
  add_profile_stage("prepare_rois_and_luts")

  local_hit_from_keys <- function(keys) {
    if (!nrow(local_coords)) return(0L)
    if (!length(keys) || anyNA(keys)) return(0L)
    as.integer(all(lut_local[keys]))
  }

  processed_tracks <- 0L
	  unique_focal_streamlines <- 0L
	  unique_local_streamlines <- 0L
	  unique_leave_no_target_streamlines <- 0L
	  unique_with_target_streamlines <- 0L
	  unique_local_target_positive_streamlines <- 0L
	  unique_local_target_negative_streamlines <- 0L
	  unique_nonlocal_target_positive_streamlines <- 0L
	  unique_nonlocal_target_negative_streamlines <- 0L
	  unique_any_target_positive_streamlines <- 0L
	  unique_any_target_negative_streamlines <- 0L
	  target_unique_counts <- integer(length(cNode_labels))
	  full_path_target_unique_counts <- integer(length(cNode_labels))
  local_search_assignment_rows <- list()
  total_focal_associations <- 0L
  total_focal_target_associations <- 0L

  count_track <- function(track, streamline_index = NA_real_) {
    mapped <- tryCatch(ftq_map_coords(track, space_used), error = function(e) NULL)
    if (is.null(mapped)) return(invisible(NULL))
    stored_track_points <- nrow(mapped)
    if (reverse_streamlines && nrow(mapped) > 1L) {
      mapped <- mapped[nrow(mapped):1L, , drop = FALSE]
    }
    if (endpoint_enabled) {
      ftq_endpoint_record_context_r(
        endpoint_accumulator,
        mapped,
        methods[[1L]],
        flip_used,
        vol_dim,
        endpoint_context_info$lut
      )
    }
    if (identical(streamline_voxelization_mode, "segment_interpolation")) {
      mapped <- ftq_densify_voxel_path(mapped)
    }

    chosen_k <- NULL
    focal_idxs <- integer(0L)
    for (mth in methods) {
      ci <- switch(mth,
        round = round(mapped),
        floor = floor(mapped),
        ceil = ceiling(mapped)
      )
      storage.mode(ci) <- "integer"
      ci <- ftq_apply_flip(ci, flip_used)
      kk <- ftq_encode(ci)
      hit_pos <- which(!is.na(kk) & lut_f[kk])
      if (length(hit_pos)) {
        chosen_k <- kk
        if (identical(focal_assignment_mode, "first_intersection")) {
          first_idx <- f_index[kk[hit_pos[1L]]]
          if (!is.na(first_idx)) focal_idxs <- first_idx
        } else {
          focal_hit_keys <- unique(kk[hit_pos])
          focal_idxs <- unique(f_index[focal_hit_keys])
          focal_idxs <- focal_idxs[!is.na(focal_idxs)]
        }
        break
      }
    }
    if (is.null(chosen_k) || !length(focal_idxs)) return(invisible(NULL))

    local_flag <- local_hit_from_keys(chosen_k)
    valid_k <- unique(chosen_k[!is.na(chosen_k)])
    full_path_target_hit <- if (length(valid_k) && length(cNode_labels)) {
      colSums(target_lut[valid_k, , drop = FALSE]) > 0L
    } else {
      rep(FALSE, length(cNode_labels))
    }
    target_assignment <- ftq_primary_target_hits_r(
      mode = target_intersection_mode,
      full_path_target_hit = full_path_target_hit,
      chosen_keys = chosen_k,
      mapped = mapped,
      two_endpoints = stored_track_points > 1L,
      lut_f = lut_f,
      target_lut = target_lut,
      axis_flip = flip_used,
      vol_dim = vol_dim,
      local_search_radius_mm = local_search_radius_mm,
      local_search_geometry = local_search_geometry
    )
    target_hit <- target_assignment$target_hit
    if (identical(target_intersection_mode, "local_search")) {
      endpoint_code <- function(endpoint) {
        if (!length(endpoint)) return(-3L)
        if (isTRUE(endpoint$ambiguous)) return(-2L)
        if (is.null(endpoint$label) || is.na(endpoint$label)) return(-1L)
        as.integer(endpoint$label) - 1L
      }
      endpoint_distance <- function(endpoint) {
        if (!length(endpoint) || is.null(endpoint$distance_mm)) return(NA_real_)
        as.numeric(endpoint$distance_mm)
      }
      assignments <- target_assignment$endpoints
      start_assignment <- if (length(assignments) >= 1L) assignments[[1L]] else list()
      end_assignment <- if (length(assignments) >= 2L) assignments[[2L]] else list()
      start_code <- endpoint_code(start_assignment)
      end_code <- endpoint_code(end_assignment)
      start_distance <- endpoint_distance(start_assignment)
      end_distance <- endpoint_distance(end_assignment)
      paired_target <- which(target_hit)
      paired_target <- if (length(paired_target)) paired_target[[1L]] else NA_integer_
      focal_distance <- NA_real_
      target_distance <- NA_real_
      max_distance <- NA_real_
      paired_full_path <- NA
      if (!is.na(paired_target)) {
        if (start_code == 0L && end_code == paired_target) {
          focal_distance <- start_distance
          target_distance <- end_distance
        } else if (end_code == 0L && start_code == paired_target) {
          focal_distance <- end_distance
          target_distance <- start_distance
        }
        if (is.finite(focal_distance) && is.finite(target_distance)) {
          max_distance <- max(focal_distance, target_distance)
        }
        paired_full_path <- isTRUE(full_path_target_hit[[paired_target]])
      }
      local_search_assignment_rows[[length(local_search_assignment_rows) + 1L]] <<-
        data.frame(
          streamline_index = as.numeric(streamline_index),
          start_status_code = start_code,
          end_status_code = end_code,
          start_distance_mm = start_distance,
          end_distance_mm = end_distance,
          start_ambiguous = isTRUE(start_assignment$ambiguous),
          end_ambiguous = isTRUE(end_assignment$ambiguous),
          paired_target_index = paired_target,
          paired_focal_distance_mm = focal_distance,
          paired_target_distance_mm = target_distance,
          paired_max_distance_mm = max_distance,
          paired_target_intersects_full_path = paired_full_path,
          stringsAsFactors = FALSE
        )
    }
    if (endpoint_enabled) {
      ftq_endpoint_record_focal_r(
        endpoint_accumulator,
        chosen_k,
        two_endpoints = stored_track_points > 1L,
        target_hit = full_path_target_hit,
        lut_f = lut_f,
        lut_f_dilated = lut_f_dilated,
        target_lut = target_lut,
        target_lut_dilated = target_lut_dilated
      )
    }

    unique_focal_streamlines <<- unique_focal_streamlines + 1L
	    unique_local_streamlines <<- unique_local_streamlines + local_flag
	    target_unique_counts <<- target_unique_counts + as.integer(target_hit)
	    full_path_target_unique_counts <<- full_path_target_unique_counts +
	      as.integer(full_path_target_hit)
	    any_target <- any(target_hit)
	    if (any_target) {
	      unique_with_target_streamlines <<- unique_with_target_streamlines + 1L
	    }
	    if (!any_target && local_flag == 0L) {
	      unique_leave_no_target_streamlines <<- unique_leave_no_target_streamlines + 1L
	    }
	    if (local_flag == 1L && any_target) {
	      unique_local_target_positive_streamlines <<- unique_local_target_positive_streamlines + 1L
	    } else if (local_flag == 1L && !any_target) {
	      unique_local_target_negative_streamlines <<- unique_local_target_negative_streamlines + 1L
	    } else if (local_flag == 0L && any_target) {
	      unique_nonlocal_target_positive_streamlines <<- unique_nonlocal_target_positive_streamlines + 1L
	    } else {
	      unique_nonlocal_target_negative_streamlines <<- unique_nonlocal_target_negative_streamlines + 1L
	    }
	    unique_any_target_positive_streamlines <<- unique_any_target_positive_streamlines + as.integer(any_target)
	    unique_any_target_negative_streamlines <<- unique_any_target_negative_streamlines + as.integer(!any_target)
	    total_focal_associations <<- total_focal_associations + length(focal_idxs)
    total_focal_target_associations <<- total_focal_target_associations +
      length(focal_idxs) * sum(target_hit)

    for (origin_idx in focal_idxs) {
      tracks_counts[origin_idx] <<- tracks_counts[origin_idx] + 1L
      local_counts[origin_idx] <<- local_counts[origin_idx] + local_flag
      if (length(cNode_labels)) {
        target_counts[origin_idx, ] <<- target_counts[origin_idx, ] + as.integer(target_hit)
      }
      if (!any_target && local_flag == 0L) {
        leave_no_target_counts[origin_idx] <<- leave_no_target_counts[origin_idx] + 1L
      }
    }
    invisible(NULL)
  }

  cpp_raw_to_voxel <- NULL
  count_chunk_cpp <- NULL
  if (identical(backend_used, "cpp")) {
    cpp_raw_to_voxel <- ftq_raw_to_voxel_affine(space_used)
    count_chunk_cpp <- function(chunk, start_index) {
      if (reverse_streamlines) {
        chunk <- lapply(chunk, ftq_reverse_track)
      }
      res <- ftq_count_chunk_cpp(
        tracks = chunk,
        raw_to_voxel = cpp_raw_to_voxel,
        methods = methods,
        axis_flip = as.logical(flip_used),
        vol_dim = vol_dim,
        f_index = f_index,
        lut_f = lut_f,
        lut_local = lut_local,
        has_local = nrow(local_coords) > 0L,
	        target_lut = target_lut,
	        n_fnode = n_fnode,
	        focal_assignment_mode = focal_assignment_mode,
	        streamline_voxelization_mode = streamline_voxelization_mode,
	        target_intersection_mode = target_intersection_mode,
	        local_search_radius_mm = local_search_radius_mm,
	        local_search_voxel_to_world =
	          local_search_geometry$voxel_to_world,
	        local_search_voxel_radius = local_search_voxel_radius,
	        endpoint_analysis = endpoint_enabled,
	        lut_f_dilated = lut_f_dilated,
	        target_lut_dilated = target_lut_dilated,
	        context_lut = endpoint_context_info$lut
	      )
      tracks_counts <<- tracks_counts + as.integer(res$tracks)
      local_counts <<- local_counts + as.integer(res$local)
      leave_no_target_counts <<- leave_no_target_counts + as.integer(res$leave_no_target)
      if (length(cNode_labels)) {
        target_counts <<- target_counts + res$targets
      }
      summary <- res$summary
      unique_focal_streamlines <<- unique_focal_streamlines +
        summary[["unique_focal_intersecting_streamlines"]]
      unique_local_streamlines <<- unique_local_streamlines +
        summary[["unique_local_streamlines"]]
      unique_leave_no_target_streamlines <<- unique_leave_no_target_streamlines +
        summary[["unique_leave_no_target_streamlines"]]
	      unique_with_target_streamlines <<- unique_with_target_streamlines +
	        summary[["unique_focal_intersecting_with_target_streamlines"]]
	      unique_local_target_positive_streamlines <<- unique_local_target_positive_streamlines +
	        summary[["unique_local_target_positive_streamlines"]]
	      unique_local_target_negative_streamlines <<- unique_local_target_negative_streamlines +
	        summary[["unique_local_target_negative_streamlines"]]
	      unique_nonlocal_target_positive_streamlines <<- unique_nonlocal_target_positive_streamlines +
	        summary[["unique_nonlocal_target_positive_streamlines"]]
	      unique_nonlocal_target_negative_streamlines <<- unique_nonlocal_target_negative_streamlines +
	        summary[["unique_nonlocal_target_negative_streamlines"]]
	      unique_any_target_positive_streamlines <<- unique_any_target_positive_streamlines +
	        summary[["unique_any_target_positive_streamlines"]]
	      unique_any_target_negative_streamlines <<- unique_any_target_negative_streamlines +
	        summary[["unique_any_target_negative_streamlines"]]
	      total_focal_associations <<- total_focal_associations +
        summary[["total_focal_voxel_streamline_associations"]]
	      total_focal_target_associations <<- total_focal_target_associations +
	        summary[["total_focal_voxel_streamline_target_associations"]]
	      target_unique_counts <<- target_unique_counts + as.integer(res$target_unique)
	      full_path_target_unique_counts <<-
	        full_path_target_unique_counts +
	        as.integer(res$full_path_target_unique)
	      if (identical(target_intersection_mode, "local_search") &&
	          is.data.frame(res$local_search_assignments) &&
	          nrow(res$local_search_assignments)) {
	        assignments <- res$local_search_assignments
	        assignments$streamline_index <-
	          as.numeric(start_index) +
	          as.numeric(assignments$chunk_streamline_index) - 1
	        assignments$chunk_streamline_index <- NULL
	        local_search_assignment_rows[[
	          length(local_search_assignment_rows) + 1L
	        ]] <<- assignments
	      }
	      if (endpoint_enabled) {
	        ftq_endpoint_merge_cpp(endpoint_accumulator, res$endpoint)
      }
      as.integer(res$processed)
    }
  }

  iterate_tracks(function(chunk, start_index, end_index) {
    if (identical(backend_used, "cpp")) {
      previous_tracks <- processed_tracks
      processed_tracks <<-
        processed_tracks + count_chunk_cpp(chunk, start_index)
      if (progress && (processed_tracks %/% progress_every > previous_tracks %/% progress_every ||
                       (!is.na(n_tracks) && processed_tracks >= n_tracks))) {
        total <- if (is.na(n_tracks)) "" else paste0(" / ", n_tracks)
        message("Processed ", processed_tracks, total, " tracks")
      }
    } else {
      for (j in seq_along(chunk)) {
        count_track(
          chunk[[j]],
          streamline_index = as.numeric(start_index) + j - 1
        )
        processed_tracks <<- processed_tracks + 1L
        if (progress && (processed_tracks %% progress_every == 0L ||
                         (!is.na(n_tracks) && processed_tracks == n_tracks))) {
          total <- if (is.na(n_tracks)) "" else paste0(" / ", n_tracks)
          message("Processed ", processed_tracks, total, " tracks")
        }
      }
    }
  })
  add_profile_stage("count_tracks")

  endpoint_tables <- if (endpoint_enabled) {
	    ftq_endpoint_finalize(
	      endpoint_accumulator,
	      target_labels = cNode_labels,
	      target_unique_counts = full_path_target_unique_counts,
      unique_focal_streamlines = unique_focal_streamlines,
      lut_f = lut_f,
      lut_f_dilated = lut_f_dilated,
      target_lut = target_lut,
      target_lut_dilated = target_lut_dilated,
      context_info = endpoint_context_info,
      vol_dim = vol_dim,
      voxel_size = ftq_endpoint_voxel_size(fNode, ref_use),
      methods = methods,
      streamline_voxelization_mode = streamline_voxelization_mode,
      track_space = space_used,
      axis_flip = flip_used
    )
  } else {
    NULL
  }
  add_profile_stage("assemble_endpoint_tables")

  local_search_assignments <- NULL
  if (identical(target_intersection_mode, "local_search")) {
    empty_assignments <- data.frame(
      streamline_index = numeric(),
      start_status_code = integer(),
      end_status_code = integer(),
      start_distance_mm = numeric(),
      end_distance_mm = numeric(),
      start_ambiguous = logical(),
      end_ambiguous = logical(),
      paired_target_index = integer(),
      paired_focal_distance_mm = numeric(),
      paired_target_distance_mm = numeric(),
      paired_max_distance_mm = numeric(),
      paired_target_intersects_full_path = logical(),
      stringsAsFactors = FALSE
    )
    local_search_assignments <- if (length(local_search_assignment_rows)) {
      out <- do.call(rbind, local_search_assignment_rows)
      rownames(out) <- NULL
      out
    } else {
      empty_assignments
    }
    status_label <- function(code) {
      code <- as.integer(code)
      out <- rep("invalid_status_code", length(code))
      out[code == -3L] <- "endpoint_unavailable"
      out[code == -2L] <- "ambiguous_nearest_node"
      out[code == -1L] <- "no_node_within_radius"
      out[code == 0L] <- "focal"
      target <- !is.na(code) & code >= 1L & code <= length(cNode_labels)
      out[target] <- cNode_labels[code[target]]
      out[is.na(code)] <- "endpoint_unavailable"
      out
    }
    local_search_assignments$start_status <- status_label(
      local_search_assignments$start_status_code
    )
    local_search_assignments$end_status <- status_label(
      local_search_assignments$end_status_code
    )
    paired_index <- as.integer(local_search_assignments$paired_target_index)
    local_search_assignments$paired_target <- NA_character_
    paired <- !is.na(paired_index) &
      paired_index >= 1L & paired_index <= length(cNode_labels)
    local_search_assignments$paired_target[paired] <-
      cNode_labels[paired_index[paired]]
    local_search_assignments$search_radius_mm <- local_search_radius_mm
    local_search_assignments <- local_search_assignments[c(
      "streamline_index",
      "start_status", "start_distance_mm", "start_ambiguous",
      "end_status", "end_distance_mm", "end_ambiguous",
      "paired_target", "paired_focal_distance_mm",
      "paired_target_distance_mm", "paired_max_distance_mm",
      "paired_target_intersects_full_path", "search_radius_mm"
    )]
  }
  add_profile_stage("assemble_local_search_assignments")

  out_df <- data.frame(
    fNode_Coords = f_labels,
    tracks = tracks_counts,
    local = local_counts,
    leave_no_target = leave_no_target_counts,
    stringsAsFactors = FALSE
  )
  if (length(cNode_labels)) {
    target_df <- as.data.frame(target_counts, stringsAsFactors = FALSE)
    names(target_df) <- cNode_labels
    out_df <- cbind(out_df, target_df)
  }
  out_df$none <- ifelse(out_df$tracks == 0L, 1L, 0L)
  rownames(out_df) <- NULL
  add_profile_stage("assemble_output")

	  unique_exclusive_fate_sum <- unique_local_target_positive_streamlines +
	    unique_local_target_negative_streamlines +
	    unique_nonlocal_target_positive_streamlines +
	    unique_nonlocal_target_negative_streamlines
	  base_roi_summary <- data.frame(
	    total_input_streamlines = if (is.na(n_tracks)) NA_real_ else as.numeric(n_tracks),
	    processed_streamlines = as.numeric(processed_tracks),
	    unique_focal_intersecting_streamlines = as.numeric(unique_focal_streamlines),
	    unique_local_streamlines = as.numeric(unique_local_streamlines),
	    unique_leave_no_target_streamlines = as.numeric(unique_leave_no_target_streamlines),
	    unique_leave_with_target_streamlines = as.numeric(unique_nonlocal_target_positive_streamlines),
	    unique_focal_intersecting_with_target_streamlines = as.numeric(unique_with_target_streamlines),
	    unique_local_target_positive_streamlines = as.numeric(unique_local_target_positive_streamlines),
	    unique_local_target_negative_streamlines = as.numeric(unique_local_target_negative_streamlines),
	    unique_nonlocal_target_positive_streamlines = as.numeric(unique_nonlocal_target_positive_streamlines),
	    unique_nonlocal_target_negative_streamlines = as.numeric(unique_nonlocal_target_negative_streamlines),
	    unique_any_target_positive_streamlines = as.numeric(unique_any_target_positive_streamlines),
	    unique_any_target_negative_streamlines = as.numeric(unique_any_target_negative_streamlines),
	    unique_exclusive_fate_sum = as.numeric(unique_exclusive_fate_sum),
	    unique_exclusive_fate_sum_matches_focal = identical(as.numeric(unique_exclusive_fate_sum), as.numeric(unique_focal_streamlines)),
	    total_focal_voxel_streamline_associations = as.numeric(total_focal_associations),
		    total_focal_voxel_streamline_target_associations = as.numeric(total_focal_target_associations),
		    focal_assignment_mode = focal_assignment_mode,
		    streamline_voxelization_mode = streamline_voxelization_mode,
		    target_intersection_mode = target_intersection_mode,
		    local_roi_defaulted_to_focal = local_roi_defaulted_to_focal,
	    selected_track_space = space_used,
	    track_voxel_size_mm = if (is.null(track_voxel_size)) NA_character_ else paste(signif(track_voxel_size, 8), collapse = ","),
    selected_axis_flip = paste(as.logical(flip_used), collapse = ","),
    selected_discretization = paste(methods, collapse = ","),
    reverse_streamlines = reverse_streamlines,
    backend = backend_used,
	    mapping_selection_rule = mapping_selection_rule,
		    stringsAsFactors = FALSE
		  )
		  if (identical(target_intersection_mode, "local_search")) {
		    base_roi_summary$local_search_radius_mm <- local_search_radius_mm
		    base_roi_summary$local_search_spatial_units <-
		      local_search_geometry$spatial_units
		    base_roi_summary$local_search_spatial_units_scale_to_mm <-
		      local_search_geometry$spatial_units_scale_to_mm
		    base_roi_summary$local_search_spatial_units_interpretation <-
		      local_search_geometry$spatial_units_interpretation
		  }
	  if (length(cNode_labels)) {
    target_summary <- as.data.frame(
      stats::setNames(as.list(as.numeric(target_unique_counts)),
                      paste0("unique_target_streamlines:", cNode_labels)),
      check.names = FALSE
    )
    if (!identical(target_intersection_mode, "full_streamline")) {
      full_path_target_summary <- as.data.frame(
        stats::setNames(
          as.list(as.numeric(full_path_target_unique_counts)),
          paste0("unique_full_path_target_streamlines:", cNode_labels)
        ),
        check.names = FALSE
      )
      target_summary <- cbind(target_summary, full_path_target_summary)
    }
    roi_summary <- cbind(base_roi_summary, target_summary)
  } else {
    roi_summary <- base_roi_summary
  }

	  attr(out_df, "fibrtraq_backend") <- backend_used
		  attr(out_df, "fibrtraq_focal_assignment_mode") <- focal_assignment_mode
		  attr(out_df, "fibrtraq_streamline_voxelization_mode") <- streamline_voxelization_mode
			  attr(out_df, "fibrtraq_target_intersection_mode") <-
			    target_intersection_mode
			  if (identical(target_intersection_mode, "local_search")) {
			    attr(out_df, "fibrtraq_local_search_radius_mm") <-
			      local_search_radius_mm
			    attr(out_df, "fibrtraq_local_search_spatial_units") <-
			      local_search_geometry$spatial_units
			    attr(
			      out_df,
			      "fibrtraq_local_search_spatial_units_scale_to_mm"
			    ) <- local_search_geometry$spatial_units_scale_to_mm
			    attr(
			      out_df,
			      "fibrtraq_local_search_spatial_units_interpretation"
			    ) <- local_search_geometry$spatial_units_interpretation
			  }
		  attr(out_df, "fibrtraq_local_roi_defaulted_to_focal") <- local_roi_defaulted_to_focal
  attr(out_df, "fibrtraq_track_space") <- space_used
  if (identical(input_format, "tck")) {
    attr(out_df, "fibrtraq_tractogram_format") <- input_format
    attr(out_df, "fibrtraq_coordinate_geometry") <- tck_geometry
  }
  attr(out_df, "fibrtraq_track_voxel_size") <- track_voxel_size
  attr(out_df, "fibrtraq_axis_flip") <- as.logical(flip_used)
  attr(out_df, "fibrtraq_discretization") <- methods
  attr(out_df, "fibrtraq_reverse_streamlines") <- reverse_streamlines
  attr(out_df, "fibrtraq_roi_summary") <- roi_summary
  attr(out_df, "fibrtraq_mapping_scores") <- mapping_scores
  attr(out_df, "fibrtraq_mapping_selection_rule") <- mapping_selection_rule
  attr(out_df, "fibrtraq_endpoint_analysis") <- endpoint_analysis
  if (!is.null(endpoint_tables)) {
    attr(out_df, "fibrtraq_endpoint_summary") <- endpoint_tables$summary
    attr(out_df, "fibrtraq_endpoint_penetration") <- endpoint_tables$penetration
    attr(out_df, "fibrtraq_endpoint_context") <- endpoint_tables$context
    attr(out_df, "fibrtraq_endpoint_metadata") <- endpoint_tables$metadata
  }
  if (!is.null(local_search_assignments)) {
    attr(out_df, "fibrtraq_local_search_assignments") <-
      local_search_assignments
  }
  if (isTRUE(profile)) {
    attr(out_df, "fibrtraq_profile") <- if (length(profile_rows)) {
      do.call(rbind, profile_rows)
    } else {
      data.frame(stage = character(), seconds = numeric(), backend = character())
    }
  }

  out_df
}

ftq_has_spatial_grid_metadata <- function(x) {
  if (!(
    inherits(x, "niftiImage") ||
      inherits(x, "nifti") ||
      inherits(x, "anlz")
  )) {
    return(FALSE)
  }
  transform <- tryCatch(unclass(RNifti::xform(x)), error = function(e) NULL)
  if (is.null(transform) || !is.matrix(transform) ||
      !identical(dim(transform), c(4L, 4L)) ||
      any(!is.finite(transform))) {
    return(FALSE)
  }
  transform_code <- attr(transform, "code", exact = TRUE)
  length(transform_code) == 1L && is.finite(transform_code) &&
    transform_code > 0
}

ftq_validate_direct_roi_coordinates <- function(x, vol_dim, input_name) {
  if (ncol(x) < 3L) {
    stop(input_name, " must have at least three coordinate columns.",
         call. = FALSE)
  }
  coords <- suppressWarnings(
    matrix(as.numeric(as.matrix(x)[, 1:3, drop = FALSE]), ncol = 3L)
  )
  if (any(!is.finite(coords))) {
    stop(input_name, " coordinates must all be finite numeric values.",
         call. = FALSE)
  }
  rounded <- round(coords)
  integer_tolerance <- sqrt(.Machine$double.eps) * pmax(1, abs(coords))
  if (any(abs(coords - rounded) > integer_tolerance)) {
    stop(input_name, " coordinates must be integer-compatible voxel indices.",
         call. = FALSE)
  }
  if (length(rounded) && any(
    rounded[, 1L] < 1L | rounded[, 1L] > vol_dim[[1L]] |
      rounded[, 2L] < 1L | rounded[, 2L] > vol_dim[[2L]] |
      rounded[, 3L] < 1L | rounded[, 3L] > vol_dim[[3L]]
  )) {
    stop(
      input_name,
      " coordinates must be within the fNode voxel grid (1..",
      paste(vol_dim, collapse = ", 1.."),
      ").",
      call. = FALSE
    )
  }
  storage.mode(rounded) <- "integer"
  colnames(rounded) <- c("x", "y", "z")
  rounded
}

ftq_validate_direct_roi_grid <- function(x, fNode, vol_dim, input_name) {
  if (is.matrix(x) || is.data.frame(x)) {
    return(ftq_validate_direct_roi_coordinates(x, vol_dim, input_name))
  }

  input_dims <- dim(x)
  if (is.null(input_dims) || length(input_dims) != 3L) {
    stop(
      input_name,
      " must be a 3-D array or an n x 3 coordinate matrix.",
      call. = FALSE
    )
  }
  input_dims <- as.integer(input_dims[1:3])
  if (!identical(input_dims, as.integer(vol_dim))) {
    stop(
      input_name,
      " dimensions must match the fNode grid (",
      input_name,
      " = ",
      paste(input_dims, collapse = "x"),
      "; fNode = ",
      paste(vol_dim, collapse = "x"),
      ").",
      call. = FALSE
    )
  }
  if (ftq_has_spatial_grid_metadata(fNode) &&
      ftq_has_spatial_grid_metadata(x) &&
      !ftq_local_search_affine_equal(fNode, x)) {
    stop(
      input_name,
      " must use the same unit-normalized affine/grid as fNode.",
      call. = FALSE
    )
  }
  x
}

ftq_reference_image_for_grid_validation <- function(ref_image) {
  if (is.character(ref_image)) {
    if (length(ref_image) != 1L || is.na(ref_image) ||
        !nzchar(trimws(ref_image))) {
      stop("ref_image must be one non-empty path or a NIfTI image.",
           call. = FALSE)
    }
    return(tryCatch(
      RNifti::readNifti(path.expand(ref_image), internal = FALSE),
      error = function(e) {
        stop("ref_image could not be read: ", conditionMessage(e),
             call. = FALSE)
      }
    ))
  }
  ref_image
}

ftq_validate_direct_reference_grid <- function(fNode, ref_image, vol_dim) {
  reference <- ftq_reference_image_for_grid_validation(ref_image)
  reference_dims <- dim(reference)
  if (is.null(reference_dims) || length(reference_dims) != 3L) {
    stop("ref_image must provide a three-dimensional voxel grid.",
         call. = FALSE)
  }
  reference_dims <- as.integer(reference_dims)
  if (!identical(reference_dims, as.integer(vol_dim))) {
    stop(
      "ref_image dimensions must match the fNode grid (ref_image = ",
      paste(reference_dims, collapse = "x"),
      "; fNode = ",
      paste(vol_dim, collapse = "x"),
      ").",
      call. = FALSE
    )
  }
  if (ftq_has_spatial_grid_metadata(fNode) &&
      ftq_has_spatial_grid_metadata(reference) &&
      !ftq_local_search_affine_equal(fNode, reference)) {
    stop(
      "ref_image must use the same unit-normalized affine/grid as fNode.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

ftq_compiled_backend_available <- function() {
  is.function(get0("ftq_count_chunk_cpp", mode = "function"))
}

ftq_resolve_backend <- function(backend) {
  backend <- match.arg(backend, c("r", "cpp", "auto"))
  if (identical(backend, "auto")) {
    if (ftq_compiled_backend_available()) "cpp" else "r"
  } else {
    backend
  }
}
