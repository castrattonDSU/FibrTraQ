#' Summarize voxel-resolved target associations and RSC
#'
#' Compute, for each target ROI column, the number of
#' focal-voxel/streamline/target associations and RSC. RSC is defined here as
#' the fraction of focal voxels with at least one focal-voxel/streamline/target
#' association for that target.
#'
#' Under the primary \code{all_focal_voxels} mode, a single unique streamline can
#' occupy multiple focal voxels. Therefore, summing a target column across focal
#' voxels yields an association total, not a unique target-positive streamline
#' count. Use \code{\link{ftq_target_summary_tables}} with the ROI-level summary
#' attached to a fresh \code{\link{FibrTraQ}} result when unique target-positive
#' focal-intersecting streamline counts are required.
#'
#' For a fresh \code{\link{FibrTraQ}} result, the target columns reflect its
#' selected \code{target_intersection_mode}: full-streamline path voxels, exact
#' endpoint voxels, or local endpoint search. This summarizer does not reassign
#' targets. Legacy \code{\link{cNode_tracks}} inputs retain their historical
#' full-streamline intersection meaning.
#'
#' The input fibr_table is the per-voxel table produced by
#' \code{\link{FibrTraQ}} and \code{\link{cNode_tracks}} with columns:
#' \itemize{
#'   \item fNode_Coords: character key "x_y_z" for the focal voxel.
#'   \item tracks: focal-voxel/streamline association count for the voxel.
#'   \item local: count of streamlines with the non-exclusive local flag.
#'   \item leave_no_target: focal-voxel/streamline associations that are
#'         non-local and have no target-ROI connection under the selected
#'         primary assignment mode.
#'   \item one column per target label: number of focal-voxel/streamline/target
#'         associations under the selected primary assignment mode.
#'   \item none: 1 if tracks == 0, else 0.
#' }
#'
#' @param fibr_table Data frame as returned by \code{\link{FibrTraQ}} or
#'   \code{\link{cNode_tracks}}.
#' @param percent Logical; if TRUE, also return rsc_percent (0 to 100).
#'   If FALSE, return rsc as a proportion on 0 to 1.
#' @param include_legacy_tracks Logical; include deprecated \code{tracks} as an
#'   alias of \code{target_associations}. This alias is retained for backward
#'   compatibility only.
#' @param warn_legacy Logical; warn when the deprecated \code{tracks} alias is
#'   returned.
#'
#' @return A data frame with one row per target (column order follows the input),
#'   with columns:
#'   \describe{
#'     \item{node}{Target label (column name).}
#'     \item{target_associations}{Number of focal-voxel/streamline/target
#'       associations for that target.}
#'     \item{focal_voxel_streamline_target_associations}{Same value as
#'       \code{target_associations}, with an explicit denominator name for
#'       publication tables.}
#'     \item{tracks}{Deprecated alias of \code{target_associations}, present only
#'       when \code{include_legacy_tracks = TRUE}.}
#'     \item{voxels_connected}{Number of focal voxels with at least one
#'       association for that target.}
#'     \item{rsc or rsc_percent}{RSC as a proportion (0 to 1) or percent (0 to 100).}
#'     \item{total_focal_voxels}{Total number of focal voxels.}
#'     \item{total_origin_voxels}{Deprecated alias of \code{total_focal_voxels},
#'       retained for backward compatibility.}
#'   }
#'
#' @examples
#' # res <- FibrTraQ(...)
#' # summary_connectivity(res, percent = TRUE, include_legacy_tracks = FALSE)
#'
#' @seealso \code{\link{FibrTraQ}}, \code{\link{cNode_tracks}},
#'   \code{\link{rsc_summary}}
#' @export
summary_connectivity <- function(fibr_table,
                                 percent = FALSE,
                                 include_legacy_tracks = TRUE,
                                 warn_legacy = include_legacy_tracks) {
  if (!is.data.frame(fibr_table)) {
    stop("fibr_table must be a data.frame.")
  }
  include_legacy_tracks <- isTRUE(include_legacy_tracks)
  warn_legacy <- isTRUE(warn_legacy)

  reserved <- c("fNode_Coords", "tracks", "local", "leave_no_target",
                "none", "cNode_Coords")
  node_cols <- setdiff(colnames(fibr_table), reserved)
  metric_col <- if (isTRUE(percent)) "rsc_percent" else "rsc"
  base_cols <- c(
    "node",
    "target_associations",
    "focal_voxel_streamline_target_associations",
    if (include_legacy_tracks) "tracks",
    "voxels_connected",
    metric_col,
    "total_focal_voxels",
    "total_origin_voxels"
  )
  maybe_warn_legacy <- function() {
    if (include_legacy_tracks && warn_legacy) {
      warning(
        "`tracks` in summary_connectivity() is deprecated and is an alias for ",
        "`target_associations`; it is not a unique target-positive streamline count. ",
        "Use `target_associations` or `focal_voxel_streamline_target_associations`, ",
        "and use ftq_target_summary_tables() for unique ROI-level streamline counts.",
        call. = FALSE
      )
    }
  }

  if (length(node_cols) == 0L) {
    out <- data.frame(
      node = character(0L),
      target_associations = numeric(0L),
      focal_voxel_streamline_target_associations = numeric(0L),
      voxels_connected = integer(0L),
      total_focal_voxels = integer(0L),
      total_origin_voxels = integer(0L),
      stringsAsFactors = FALSE
    )
    if (include_legacy_tracks) out$tracks <- numeric(0L)
    if (isTRUE(percent)) out$rsc_percent <- numeric(0L) else out$rsc <- numeric(0L)
    maybe_warn_legacy()
    return(out[, base_cols, drop = FALSE])
  }

  n_vox <- nrow(fibr_table)
  if (n_vox == 0L) {
    out <- data.frame(
      node = node_cols,
      target_associations = numeric(length(node_cols)),
      focal_voxel_streamline_target_associations = numeric(length(node_cols)),
      voxels_connected = integer(length(node_cols)),
      total_focal_voxels = integer(length(node_cols)),
      total_origin_voxels = integer(length(node_cols)),
      stringsAsFactors = FALSE
    )
    if (include_legacy_tracks) out$tracks <- numeric(length(node_cols))
    if (isTRUE(percent)) out$rsc_percent <- numeric(length(node_cols)) else out$rsc <- numeric(length(node_cols))
    maybe_warn_legacy()
    return(out[, base_cols, drop = FALSE])
  }

  mat <- as.matrix(fibr_table[, node_cols, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")

  association_totals <- colSums(mat, na.rm = TRUE)
  vox_connected <- colSums(mat > 0, na.rm = TRUE)
  rsc_prop      <- vox_connected / n_vox

  out <- data.frame(
    node = node_cols,
    target_associations = as.numeric(association_totals),
    focal_voxel_streamline_target_associations = as.numeric(association_totals),
    voxels_connected = as.integer(vox_connected),
    total_focal_voxels = rep.int(n_vox, length(node_cols)),
    total_origin_voxels = rep.int(n_vox, length(node_cols)),
    stringsAsFactors = FALSE
  )
  if (include_legacy_tracks) out$tracks <- out$target_associations
  if (isTRUE(percent)) {
    out$rsc_percent <- 100 * rsc_prop
  } else {
    out$rsc <- rsc_prop
  }
  out <- out[, base_cols, drop = FALSE]

  rownames(out) <- NULL
  maybe_warn_legacy()
  out
}
