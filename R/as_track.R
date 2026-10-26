#' Normalize a streamline (track) to a standard structure
#'
#' Convert various single‑track representations into a standard list with
#' fields \code{num_points} (integer) and \code{coords} (n x 3 numeric matrix).
#' This helper makes downstream code robust to different tractography sources.
#'
#' Accepted inputs (single track):
#' \itemize{
#'   \item A list with \code{num_points} and \code{coords}.
#'   \item A matrix or data.frame with at least 3 columns (first 3 taken as
#'         x, y, z in voxel index space).
#' }
#' (When reading \code{.trk} files via \code{freesurferformats::read.dti.trk()},
#' you typically iterate over the returned \code{$tracks} list, whose elements
#' are n x 3 coordinate matrices suitable for this function.)
#'
#' The function does not round to integers; rounding to voxel indices is handled
#' by \code{\link{iso_tracks}}.
#'
#' @param track Single streamline. See “Accepted inputs” above.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{num_points}}{Integer number of points along the track.}
#'   \item{\code{coords}}{n x 3 numeric matrix of coordinates (columns x, y, z).}
#' }
#'
#' @examples
#' as_track(matrix(c(1,2,3, 1,2,4, 1,2,5), ncol = 3, byrow = TRUE))
#'
#' @seealso \code{\link{iso_tracks}},
#'   \code{freesurferformats::read.dti.trk()}
#' @export
as_track <- function(track) {
  # Case 1: already a list with the expected fields ---------------------------
  if (is.list(track) && !is.null(track$num_points) && !is.null(track$coords)) {
    coords <- as.matrix(track$coords)
    if (ncol(coords) < 3L || nrow(coords) < 1L) {
      stop("Invalid 'coords': must be an n x 3 matrix with n >= 1.")
    }
    coords <- coords[, 1:3, drop = FALSE]
    storage.mode(coords) <- "numeric"
    colnames(coords) <- c("x","y","z")
    np <- as.integer(track$num_points)
    if (is.na(np) || np < 1L) np <- nrow(coords)
    return(list(num_points = np, coords = coords))
  }

  # Case 2: matrix/data.frame of coordinates ----------------------------------
  if (is.matrix(track) || is.data.frame(track)) {
    coords <- as.matrix(track)
    if (ncol(coords) < 3L || nrow(coords) < 1L) {
      stop("Coordinate input must be an n x 3 matrix/data.frame with n >= 1.")
    }
    coords <- coords[, 1:3, drop = FALSE]
    storage.mode(coords) <- "numeric"
    colnames(coords) <- c("x","y","z")
    return(list(num_points = nrow(coords), coords = coords))
  }

  stop("Unsupported track input: provide a list(num_points, coords) or an n x 3 matrix/data.frame.")
}
