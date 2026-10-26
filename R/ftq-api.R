#' FibrTraQ (branded entry point)
#'
#' @description Branded wrapper around \code{\link{FibrTraQ}}.
#' @param ... Additional arguments forwarded to \code{\link{FibrTraQ}}.
#'   This is primarily for backward compatibility. Unknown arguments will error
#'   in \code{FibrTraQ()} because it does not accept arbitrary \code{...}.
#' @return Per-voxel table of focal-voxel target associations.
#' @rdname FibrTraQ
#' @export
fibrtraq <- function(...) FibrTraQ(...)

#' ROI voxel coordinates (branded)
#'
#' @description Wrapper around \code{\link{iso_coords}}.
#' @param ... Additional arguments forwarded to \code{\link{iso_coords}}.
#'   Unknown arguments will error in \code{iso_coords()} because it does not
#'   accept arbitrary \code{...}.
#' @return Integer matrix with columns \code{x}, \code{y}, \code{z}.
#' @export
ftq_roi_coords <- function(...) iso_coords(...)

#' Tracks intersecting an ROI (branded)
#'
#' @description Wrapper around \code{\link{iso_tracks}}.
#' @inheritParams iso_tracks
#' @return If the track intersects \code{fNode}, a list with elements
#'   \code{coords}, \code{tracks}, and \code{local}. Otherwise \code{NULL}.
#' @export
ftq_tracks_through <- function(track,
                               fNode,
                               localROI = NULL,
                               verbose = FALSE,
                               track_space = c("voxel","voxelmm","vox10","world","trackvis","auto"),
                               ref_image = NULL,
                               discretization = c("round","floor","ceil"),
                               axis_flip = c(FALSE,FALSE,FALSE),
                               vol_dim = NULL,
                               track_affine = NULL,
                               track_voxel_size = NULL) {
  iso_tracks(track = track,
             fNode = fNode,
             localROI = localROI,
             verbose = verbose,
             track_space = track_space,
             ref_image = ref_image,
             discretization = discretization,
             axis_flip = axis_flip,
             vol_dim = vol_dim,
             track_affine = track_affine,
             track_voxel_size = track_voxel_size)
}

#' Group tracks by focal voxel (branded)
#'
#' @description Wrapper around \code{\link{track_relister}}.
#' @param ... Additional arguments forwarded to \code{\link{track_relister}}.
#'   Unknown arguments will error if not supported by \code{track_relister()}.
#' @return See \code{\link{track_relister}}.
#' @export
ftq_group_by_voxel <- function(...) track_relister(...)

#' Count target-ROI intersections per focal voxel (branded)
#'
#' @description Wrapper around \code{\link{cNode_tracks}}.
#' @param ... Additional arguments forwarded to \code{\link{cNode_tracks}}.
#'   Unknown arguments will error if not supported by \code{cNode_tracks()}.
#' @return See \code{\link{cNode_tracks}}.
#' @export
ftq_target_hits <- function(...) cNode_tracks(...)

#' Summarize voxel-resolved target associations and RSC
#'
#' @description Wrapper around \code{\link{summary_connectivity}}.
#' @param ... Additional arguments forwarded to \code{\link{summary_connectivity}}.
#'   Unknown arguments will error in \code{summary_connectivity()} because it does not
#'   accept arbitrary \code{...}.
#' @return See \code{\link{summary_connectivity}}.
#' @export
rsc_summary <- function(...) summary_connectivity(...)
