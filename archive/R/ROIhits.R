#' Legacy voxel-wise target-intersection table (compatibility shim)
#'
#' `ROIhits()` is retained for backward compatibility and now delegates to
#' `FibrTraQ()` (the modern pipeline). It computes, for each focal voxel, which
#' targets have at least one streamline intersection (binary per-target flags) and whether
#' the voxel's streamlines remain entirely within the *local* ROI.
#'
#' For legacy behavior, the *local* ROI is taken to be the focal ROI
#' (`localROI = fNode`), which matches historical usage of `ROIhits()`.
#'
#' @param fNode 3-D binary array (e.g., NIfTI volume read to an array) defining
#'   the focal ROI.
#' @param tracks Tractography object or list as returned by
#'   \code{freesurferformats::read.dti.trk()}, i.e., a named list containing
#'   \code{$tracks} (each a matrix of xyz points). A plain list of tracks is
#'   also accepted.
#' @param cNode_list List of 3-D binary arrays for each target ROI.
#'
#' @return A data.frame with one row per focal voxel and columns:
#'   \itemize{
#'     \item \strong{Target-label columns}: 1 if any streamline associated with the
#'       voxel intersects that target, else 0.
#'     \item \strong{"Voxel Coords"}: character key "x_y_z" for the focal voxel.
#'     \item \strong{"Local"}: 1 if all track points for that voxel fall within
#'       \code{fNode} (legacy behavior), else 0.
#'   }
#'
#' @examples
#' # Example (structure only):
#' # roi_df <- ROIhits(fNode = fnode_array,
#' #                   tracks = freesurferformats::read.dti.trk("subj.trk"),
#' #                   cNode_list = list(T1 = t1_array, T2 = t2_array))
#'
#' @seealso \code{\link{FibrTraQ}}, \code{\link{fibrtraq}},
#'   \code{\link{ftq_target_hits}}, \code{freesurferformats::read.dti.trk()}
#' @export
ROIhits <- function(fNode, tracks, cNode_list) {
  # Derive target labels (fall back to target_1, target_2, ...)
  labs <- names(cNode_list)
  if (is.null(labs) || any(!nzchar(labs))) {
    labs <- paste0("target_", seq_along(cNode_list))
  }
  # Ensure syntactically valid, unique column names
  labs <- make.names(labs, unique = TRUE)

  # Run the modern pipeline (legacy local behavior: localROI == fNode)
  out <- FibrTraQ(
    fNode         = fNode,
    localROI      = fNode,
    cNode_list    = cNode_list,
    cNode_labels  = labs,
    all_tracks    = tracks
  )

  # If the modern pipeline returned nothing, keep legacy NULL
  if (is.null(out) || !is.data.frame(out) || nrow(out) == 0L) {
    return(out)
  }

  # Convert per-target counts to binary flags (>=1 -> 1), keep legacy column names
  cn <- labs
  cn <- intersect(cn, colnames(out))
  hit_bin <- as.data.frame((out[, cn, drop = FALSE] > 0) + 0L, stringsAsFactors = FALSE)

  # Assemble legacy-shaped output: <targets>, "Voxel Coords", "Local"
  res <- cbind(hit_bin,
               "Voxel Coords" = out$fNode_Coords,
               Local = as.integer(out$local))
  rownames(res) <- NULL
  res
}
