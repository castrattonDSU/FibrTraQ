#' Group tracks by focal voxel
#'
#' Re-index a list of streamline objects (from \code{iso_tracks()}) by the focal
#' voxel they intersect. Returns (1) a list of per-voxel track groups ordered to
#' match the rows of \code{fNode}, and (2) a per-voxel summary data frame with
#' focal-voxel/streamline association counts and streamlines with the
#' non-exclusive local flag.
#'
#' The focal-voxel key is encoded as an underscore-delimited string "x_y_z" to
#' ensure stable joins with downstream tables.
#'
#' @param iso_track_list List of streamline objects produced by
#'   \code{\link{iso_tracks}()} (NULLs are allowed and ignored).
#' @param fNode Integer matrix of focal ROI voxel coordinates with columns
#'   \code{x}, \code{y}, \code{z}, typically from \code{\link{iso_coords}()}.
#'
#' @return A list of length 2:
#' \enumerate{
#'   \item \code{per_voxel_groups}: list where element \code{i} contains a list
#'         of streamline objects whose \code{$coords} match \code{fNode[i, ]}.
#'   \item \code{summary}: data.frame with columns
#'         \code{fNode_Coords} (character key "x_y_z"),
#'         \code{tracks} (integer count of focal-voxel/streamline associations),
#'         \code{local} (integer count of streamlines with the non-exclusive
#'         local flag).
#' }
#'
#' @examples
#' # Suppose iso_track_list contains results from iso_tracks(), and fNode is
#' # an integer matrix of ROI voxel coordinates:
#' # out <- track_relister(iso_track_list, fNode)
#' # str(out[[1]])  # per-voxel groups
#' # head(out[[2]]) # per-voxel summary table
#'
#' @seealso \code{\link{iso_coords}}, \code{\link{iso_tracks}}, \code{\link{cNode_tracks}}
#' @export
track_relister <- function(iso_track_list, fNode) {
  # Validate fNode
  if (is.null(fNode) || !is.matrix(fNode) || ncol(fNode) < 3L) {
    stop("fNode must be an integer matrix with at least three columns: x, y, z.")
  }
  fn <- as.matrix(fNode)[, 1:3, drop = FALSE]
  storage.mode(fn) <- "integer"
  colnames(fn) <- c("x", "y", "z")

  n <- nrow(fn)
  if (n == 0L) {
    per_voxel_groups <- list()
    summary <- data.frame(
      fNode_Coords = character(0L),
      tracks = integer(0L),
      local  = integer(0L),
      stringsAsFactors = FALSE
    )
    return(list(per_voxel_groups, summary))
  }

  # Keys for fNode rows and an index map for fast assignment
  f_labels <- apply(fn, 1L, function(v) paste(v[1], v[2], v[3], sep = "_"))
  index_map <- structure(seq_len(n), names = f_labels)

  # Initialize empty groups aligned to fNode row order
  per_voxel_groups <- vector(mode = "list", length = n)

  # Helper to coerce '$coords' from an iso_tracks() result to a 1x3 integer vector
  coords_to_key <- function(coords) {
    if (is.null(coords)) return(NA_character_)
    if (is.matrix(coords)) {
      if (ncol(coords) < 3L || nrow(coords) < 1L) return(NA_character_)
      vv <- coords[1, 1:3, drop = TRUE]
    } else {
      vv <- coords
    }
    vv <- as.integer(vv[1:3])
    paste(vv[1], vv[2], vv[3], sep = "_")
  }

  # Assign each non-NULL track to its focal-voxel group
  if (is.list(iso_track_list) && length(iso_track_list)) {
    for (obj in iso_track_list) {
      if (is.null(obj)) next
      k <- coords_to_key(obj$coords)
      if (!is.na(k)) {
        i <- index_map[[k]]
        if (!is.null(i)) {
          per_voxel_groups[[i]] <- c(per_voxel_groups[[i]], list(obj))
        }
      }
    }
  }

  # Build per-voxel summary
  tracks_counts <- integer(n)
  local_counts  <- integer(n)
  for (i in seq_len(n)) {
    li <- per_voxel_groups[[i]]
    tracks_counts[i] <- length(li)
    if (tracks_counts[i] > 0L) {
      local_counts[i] <- sum(vapply(li, function(x) {
        v <- x$local
        if (is.null(v)) 0L else as.integer(v[1])
      }, integer(1L)))
    }
  }

  summary <- data.frame(
    fNode_Coords = f_labels,
    tracks = tracks_counts,
    local  = local_counts,
    stringsAsFactors = FALSE
  )

  list(per_voxel_groups, summary)
}
