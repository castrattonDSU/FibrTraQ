#' Extract voxel coordinates from an ROI
#'
#' Given a binary or label ROI volume read with \code{freesurferformats} (NIfTI array),
#' return the coordinates for voxels that belong to the ROI. By default any nonzero
#' voxel is considered part of the ROI (robust for masks that are not strictly 1/0).
#'
#' If a matrix/data.frame of coordinates is supplied, the first three columns are
#' returned unchanged (coerced to a numeric matrix).
#'
#' @param input_ROI Either a 3-D array (NIfTI data) or an \code{n x 3} coordinate matrix.
#' @param threshold Either "nonzero" (default) or a numeric value; when numeric,
#'   voxels with values >= \code{threshold} are included.
#' @return An \code{n x 3} integer matrix with columns \code{x,y,z}.
#' @examples
#' # roi_coords <- iso_coords(my_roi)  # any nonzero voxels
#' # roi_coords <- iso_coords(my_roi, threshold = 0.5)  # numeric threshold
#' @export
iso_coords <- function(input_ROI, threshold = "nonzero") {
  # If coords already
  if (is.matrix(input_ROI) || is.data.frame(input_ROI)) {
    coords <- as.matrix(input_ROI)[, 1:3, drop = FALSE]
    storage.mode(coords) <- "integer"
    colnames(coords) <- c("x","y","z")
    return(coords)
  }

  # Must be a 3D array
  if (length(dim(input_ROI)) != 3L) stop("input_ROI must be a 3-D array or an n x 3 coordinate matrix")

  sel <- if (identical(threshold, "nonzero")) {
    which(input_ROI != 0, arr.ind = TRUE)
  } else if (is.numeric(threshold) && length(threshold) == 1L) {
    which(input_ROI >= threshold, arr.ind = TRUE)
  } else {
    stop("threshold must be \"nonzero\" or a single numeric value")
  }

  if (nrow(sel) == 0L) {
    # return empty 3-col matrix with correct type/colnames
    out <- matrix(integer(0L), nrow = 0L, ncol = 3L)
    colnames(out) <- c("x","y","z")
    return(out)
  }

  storage.mode(sel) <- "integer"
  colnames(sel) <- c("x","y","z")
  sel
}
