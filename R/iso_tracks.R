#' Identify tracks intersecting an ROI
#'
#' Determine whether a single streamline (track) intersects a focal ROI
#' and, if so, return (i) the matched focal-voxel coordinate, (ii) the full
#' track in voxel-index coordinates, and (iii) a local flag indicating whether
#' all track points lie within a provided local ROI.
#'
#' Track coordinates may be provided in different spaces:
#' - "voxel": coordinates are already voxel indices.
#' - "voxelmm": coordinates are TrackVis voxel-mm coordinates measured from
#'   the voxel corner. They are divided by `track_voxel_size` and shifted by
#'   0.5 to obtain 1-based voxel-center indices.
#' - "vox10": coordinates must be multiplied by 10 before discretization
#'   (legacy behavior used in earlier FibrTraQ development).
#' - "world": coordinates are in world millimetres and are converted to voxel
#'   indices using the reference NIfTI image. Metre and micrometre NIfTI
#'   spatial units are normalized automatically.
#' - "trackvis": coordinates are raw TrackVis points and are first transformed
#'   to world millimetres with `track_affine`, then converted to voxel indices.
#' - "auto": try mappings in the order vox10 -> voxel -> world -> trackvis
#'   (when required references are provided) and accept the first mapping that
#'   yields an fNode intersection.
#'
#' Discretization from floating-point coordinates to voxel indices is controlled
#' by `discretization`. You may pass a single method ("round", "floor", "ceil"),
#' a character vector of methods (tried in order), or "any" (meaning try all three).
#' A hit is accepted if any method yields an intersection. The `local` flag is set
#' to 1 if, under any method, all points fall within the `localROI`.
#'
#' Some tractography/ROI pairs require a flip in one or more voxel axes to align
#' coordinate conventions. Use `axis_flip` together with `vol_dim` to apply flips
#' after discretization. For example, if `axis_flip` is `c(FALSE, FALSE, TRUE)`
#' (flip_z), then z is transformed as \code{z <- nz - z + 1} where \code{nz} is
#' the z dimension of the volume.
#'
#' @param track A single track as `list(num_points, coords[n,3])` or a 3-col
#'   matrix/data.frame of coordinates (coerced via [as_track()]).
#' @param fNode ROI voxel coordinates (matrix with 3 cols) or a 3-D array mask.
#'   Arrays are converted via [iso_coords()].
#' @param localROI Optional local ROI coordinates (matrix) or array mask.
#'   If omitted or empty, the `local` flag will be 0.
#' @param verbose Logical; emit diagnostic messages. Default `FALSE`.
#' @param track_space One of `"voxel"` (default), `"voxelmm"`, `"vox10"`,
#'   `"world"`, `"trackvis"`, or `"auto"`.
#' @param ref_image A NIfTI reference (path or RNifti image) required if
#'   `track_space` is `"world"` or `"trackvis"` and used if `track_space = "auto"`
#'   considers either mapping.
#' @param discretization One or more of `"round"`, `"floor"`, `"ceil"`, or `"any"`.
#'   If `"any"` is present, all three methods are tried.
#' @param axis_flip Logical length-3 vector \code{c(flip_x, flip_y, flip_z)}.
#'   Default no flips.
#' @param vol_dim Integer length-3 vector giving the voxel dimensions
#'   \code{c(nx, ny, nz)}. Required if any element of `axis_flip` is TRUE and
#'   dims cannot be inferred.
#' @param track_affine Optional 4x4 affine mapping raw TrackVis coordinates to
#'   world coordinates. Required when `track_space = "trackvis"`.
#' @param track_voxel_size Numeric length-3 TrackVis header voxel sizes in mm.
#'   Required when `track_space = "voxelmm"`.
#'
#' @return `NULL` if no intersection; otherwise a list with:
#' \describe{
#'   \item{coords}{1x3 integer matrix with the matched focal voxel.}
#'   \item{tracks}{data.frame of integer voxel indices for the track.}
#'   \item{local}{integer 0/1 non-exclusive flag indicating that stored sample
#'     points are local under the selected voxelization.}
#' }
#' @seealso [iso_coords()], [as_track()]
#' @importFrom RNifti worldToVoxel
#' @export
iso_tracks <- function(track,
                       fNode,
                       localROI = NULL,
                       verbose = FALSE,
                       track_space = c("voxel", "voxelmm", "vox10", "world", "trackvis", "auto"),
                       ref_image = NULL,
                       discretization = c("round","floor","ceil"),
                       axis_flip = c(FALSE,FALSE,FALSE),
                       vol_dim = NULL,
                       track_affine = NULL,
                       track_voxel_size = NULL) {

  track_space <- match.arg(track_space)

  # Discretization methods
  if (length(discretization) == 1L && identical(discretization, "any")) {
    methods <- c("round","floor","ceil")
  } else if (any(discretization == "any")) {
    methods <- c("round","floor","ceil")
  } else {
    methods <- unique(match.arg(discretization,
                                c("round","floor","ceil"),
                                several.ok = TRUE))
  }

  # Capture dim guesses before we coerce ROI arrays to coordinate matrices
  vol_guess <- NULL
  if (is.null(vol_dim)) {
    if (!is.null(dim(fNode)) && length(dim(fNode)) == 3L) {
      vol_guess <- dim(fNode)[1:3]
    } else if (!is.null(dim(localROI)) && length(dim(localROI)) == 3L) {
      vol_guess <- dim(localROI)[1:3]
    } else if (!is.null(ref_image) && !is.character(ref_image) &&
               !is.null(dim(ref_image)) && length(dim(ref_image)) >= 3L) {
      vol_guess <- dim(ref_image)[1:3]
    }
  }

  # axis_flip validation
  if (is.null(axis_flip)) axis_flip <- c(FALSE,FALSE,FALSE)
  if (!is.logical(axis_flip) || length(axis_flip) != 3L) {
    stop("iso_tracks: axis_flip must be a logical vector of length 3: c(flip_x, flip_y, flip_z).")
  }
  need_flip <- any(axis_flip)
  if (!is.null(track_affine)) {
    track_affine <- as.matrix(track_affine)
    if (!identical(dim(track_affine), c(4L, 4L)) || anyNA(track_affine)) {
      stop("iso_tracks: track_affine must be a numeric 4x4 matrix.")
    }
    storage.mode(track_affine) <- "numeric"
  }
  if (!is.null(track_voxel_size)) {
    track_voxel_size <- ftq_validate_track_voxel_size(
      track_voxel_size,
      caller = "iso_tracks"
    )
  }
  if (identical(track_space, "voxelmm") && is.null(track_voxel_size)) {
    stop("iso_tracks: track_space = 'voxelmm' requires track_voxel_size from the TrackVis header.")
  }

  # vol_dim validation (only if needed)
  if (need_flip) {
    if (is.null(vol_dim) && !is.null(vol_guess)) vol_dim <- vol_guess
    if (is.null(vol_dim) || length(vol_dim) < 3L) {
      stop("iso_tracks: vol_dim (length 3) is required when axis_flip contains TRUE and dims cannot be inferred.")
    }
    vol_dim <- as.integer(vol_dim[1:3])
    if (anyNA(vol_dim) || any(vol_dim <= 0L)) stop("iso_tracks: vol_dim must be positive integers.")
  }

  trk <- tryCatch(as_track(track), error = function(e) NULL)
  if (is.null(trk) || is.null(trk$coords) || ncol(as.matrix(trk$coords)) < 3L) {
    if (isTRUE(verbose)) message("iso_tracks: invalid track input")
    return(NULL)
  }

  # Normalize ROI inputs to integer coordinate matrices
  to_coords <- function(x) {
    if (is.null(x)) {
      return(matrix(integer(0), ncol = 3L,
                    dimnames = list(NULL, c("x","y","z"))))
    }
    if (length(dim(x)) == 3L) x <- iso_coords(x)
    xm <- as.matrix(x)[, 1:3, drop = FALSE]
    storage.mode(xm) <- "integer"
    colnames(xm) <- c("x","y","z")
    xm
  }
  fNode    <- to_coords(fNode)
  localROI <- to_coords(localROI)

  if (!nrow(fNode)) {
    if (isTRUE(verbose)) message("iso_tracks: empty fNode")
    return(NULL)
  }

  # Raw coordinates from the track (first 3 columns)
  coords0 <- as.matrix(trk$coords)[, 1:3, drop = FALSE]

  # Helper: string key for fast set operations
  key <- function(m) paste(m[,1], m[,2], m[,3], sep = ",")

  f_keys <- key(fNode)
  l_keys <- if (nrow(localROI)) key(localROI) else character(0)

  # Apply axis flips to integer voxel coordinates
  apply_flip <- function(m) {
    if (!need_flip) return(m)
    mm <- m
    if (axis_flip[1]) mm[,1] <- vol_dim[1] - mm[,1] + 1L
    if (axis_flip[2]) mm[,2] <- vol_dim[2] - mm[,2] + 1L
    if (axis_flip[3]) mm[,3] <- vol_dim[3] - mm[,3] + 1L
    storage.mode(mm) <- "integer"
    colnames(mm) <- c("x","y","z")
    mm
  }

  # Discretization helper
  discretize <- function(m, mth) {
    out <- switch(mth,
                  round = round(m),
                  floor = floor(m),
                  ceil  = ceiling(m))
    storage.mode(out) <- "integer"
    colnames(out) <- c("x","y","z")
    out <- apply_flip(out)
    out
  }

  # Candidate spaces (for auto)
  spaces <- if (track_space == "auto") {
    sp <- c(if (!is.null(track_voxel_size)) "voxelmm", "vox10", "voxel")
    if (!is.null(ref_image)) sp <- c(sp, "world")
    if (!is.null(ref_image) && !is.null(track_affine)) sp <- c(sp, "trackvis")
    sp
  } else {
    track_space
  }

  chosen <- NULL
  chosen_coords <- NULL

  # Try each mapping, then each discretization, until a hit is found
  for (sp in spaces) {
    coords <- coords0

    if (sp == "voxelmm") {
      coords <- ftq_trackvis_voxelmm_to_voxel(
        coords0,
        track_voxel_size,
        caller = "iso_tracks"
      )
    } else if (sp == "vox10") {
      coords <- coords0 * 10
    } else if (sp == "world") {
      if (is.null(ref_image)) {
        if (isTRUE(verbose)) message("iso_tracks: ref_image missing; skipping 'world'")
        next
      }
      coords <- ftq_world_mm_to_voxel(coords0, ref_image)
    } else if (sp == "trackvis") {
      if (is.null(ref_image) || is.null(track_affine)) {
        if (isTRUE(verbose)) message("iso_tracks: ref_image or track_affine missing; skipping 'trackvis'")
        next
      }
      pts_mm <- cbind(coords0, 1) %*% t(track_affine)
      coords <- ftq_world_mm_to_voxel(
        pts_mm[, 1:3, drop = FALSE],
        ref_image
      )
    } else if (sp == "voxel") {
      coords <- coords0
    } else {
      stop("iso_tracks: unknown track_space: ", sp)
    }

    for (mth in methods) {
      tmat  <- discretize(coords, mth)
      tkeys <- key(tmat)
      if (any(tkeys %in% f_keys)) {
        chosen <- list(method = mth, tmat = tmat, tkeys = tkeys)
        chosen_coords <- coords
        break
      }
    }

    if (!is.null(chosen)) break
  }

  if (is.null(chosen)) {
    if (isTRUE(verbose)) {
      message("iso_tracks: no fNode intersection. track_space tried: ",
              paste(spaces, collapse = ","),
              "; methods: ", paste(methods, collapse = ","),
              if (need_flip) paste0("; axis_flip=", paste(axis_flip, collapse=",")) else "")
    }
    return(NULL)
  }

  # Matched focal voxel
  first_hit <- which(chosen$tkeys %in% f_keys)[1]
  f_match   <- match(chosen$tkeys[first_hit], f_keys, nomatch = 0L)
  origin_vox <- fNode[f_match, , drop = FALSE]

  # Local flag: 1 if ANY method yields all points within localROI, using the chosen mapping
  local_flag <- 0L
  if (nrow(localROI)) {
    for (mth in methods) {
      tmat_m <- discretize(chosen_coords, mth)
      if (all(key(tmat_m) %in% l_keys)) {
        local_flag <- 1L
        break
      }
    }
  }

  tracks_df <- data.frame(x = chosen$tmat[,1],
                          y = chosen$tmat[,2],
                          z = chosen$tmat[,3],
                          stringsAsFactors = FALSE)

  list(coords = origin_vox, tracks = tracks_df, local = local_flag)
}
