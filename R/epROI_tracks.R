#' Count target-ROI intersections per focal voxel
#'
#' For each focal-voxel's set of streamlines (from \code{track_relister()}),
#' count how many focal-voxel/streamline associations intersect each target ROI.
#' Target ROIs are
#' supplied as coordinate matrices (x, y, z) in the same voxel index space as
#' the streamline coordinates stored with each track.
#'
#' The function updates a per-voxel summary table by adding one column per
#' target label and filling in the number of associations for streamlines that
#' intersect that target. \code{leave_no_target} counts associations that are not
#' local and do not intersect any target. A "none" column is also provided (1 if
#' the voxel has zero tracks,
#' else 0).
#'
#' Assumptions:
#' \itemize{
#'   \item \code{track_relist} is the output of \code{track_relister()}, i.e.,
#'         a list where \code{[[1]]} is the list of per-voxel track groups and
#'         \code{[[2]]} is a data.frame with at least \code{fNode_Coords},
#'         \code{tracks}, and \code{local}.
#'   \item Track coordinates are integer voxel indices (or have already been
#'         rounded to nearest integers upstream).
#' }
#'
#' @param cNode_coords_list List of integer matrices (columns x, y, z), one per
#'   target ROI; typically created by \code{\link{iso_coords}}.
#' @param cNode_labels Character vector of target ROI labels; must match the
#'   length and order of \code{cNode_coords_list}. Labels should be unique and
#'   valid column names.
#' @param track_relist Output from \code{\link{track_relister}}, i.e., a list
#'   where \code{[[1]]} is the per-voxel track grouping and \code{[[2]]} is the
#'   voxel summary data.frame.
#'
#' @return A data.frame with one row per focal voxel including columns:
#'   \code{fNode_Coords} (character key "x_y_z"), \code{tracks}
#'   (focal-voxel/streamline associations), \code{local} (count of streamlines
#'   with the non-exclusive local flag), \code{leave_no_target} (count of
#'   non-local associations with no target-ROI intersection), one column per
#'   target label with the number of focal-voxel/streamline/target associations,
#'   and \code{none} (1 if \code{tracks == 0}, else 0).
#'
#' @examples
#' # cn_df <- cNode_tracks(cnodes, labels, relist)
#'
#' @seealso \code{\link{iso_coords}}, \code{\link{track_relister}},
#'   \code{\link{rsc_summary}}
#' @export
cNode_tracks <- function(cNode_coords_list, cNode_labels, track_relist) {
  # Basic input checks ---------------------------------------------------------
  if (is.null(cNode_coords_list) || length(cNode_coords_list) == 0L) return(NULL)
  if (is.null(cNode_labels) || length(cNode_labels) != length(cNode_coords_list)) {
    stop("cNode_labels must have the same length as cNode_coords_list.")
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
  if (is.null(track_relist) || length(track_relist) < 2L) {
    stop("track_relist must be a list with [[1]] (per-voxel groups) and [[2]] (summary data.frame).")
  }

  per_voxel_groups <- track_relist[[1]]
  coord_match_df   <- as.data.frame(track_relist[[2]], stringsAsFactors = FALSE)

  # Ensure summary table has one numeric column per target label ---------------
  for (lab in cNode_labels) {
    if (!lab %in% names(coord_match_df)) {
      coord_match_df[[lab]] <- 0L
    } else {
      if (!is.numeric(coord_match_df[[lab]])) {
        coord_match_df[[lab]] <- as.integer(coord_match_df[[lab]])
      }
      coord_match_df[[lab]][is.na(coord_match_df[[lab]])] <- 0L
    }
  }
  coord_match_df$leave_no_target <- 0L

  # Helpers --------------------------------------------------------------------
  norm_xyz <- function(m) {
    if (is.null(m) || (is.matrix(m) && nrow(m) == 0L)) {
      return(matrix(integer(0L), nrow = 0L, ncol = 3L,
                    dimnames = list(NULL, c("x","y","z"))))
    }
    mm <- as.matrix(m)[, 1:3, drop = FALSE]
    storage.mode(mm) <- "integer"
    colnames(mm) <- c("x","y","z")
    mm
  }
  key_of <- function(m) paste(m[,1], m[,2], m[,3], sep = ",")

  # Precompute key-sets for all targets ---------------------------------------
  cnode_keys <- lapply(cNode_coords_list, function(mat) {
    mm <- norm_xyz(mat)
    if (nrow(mm) == 0L) character(0L) else unique(key_of(mm))
  })
  names(cnode_keys) <- cNode_labels
  order_output <- function(df) {
    cols <- c("fNode_Coords", "tracks", "local", "leave_no_target",
              cNode_labels, "none")
    cols <- cols[cols %in% names(df)]
    extra <- setdiff(names(df), cols)
    df[, c(cols, extra), drop = FALSE]
  }

  # Iterate over focal-voxel groups -------------------------------------------
  if (!is.list(per_voxel_groups) || length(per_voxel_groups) == 0L) {
    # Still finalize "none" column and return
    coord_match_df$local[is.na(coord_match_df$local)] <- 0L
    coord_match_df$leave_no_target[is.na(coord_match_df$leave_no_target)] <- 0L
    coord_match_df$none <- ifelse(is.na(coord_match_df$tracks) | coord_match_df$tracks == 0L, 1L, 0L)
    return(order_output(coord_match_df))
  }

  for (m in seq_along(per_voxel_groups)) {
    voxel_tracks <- per_voxel_groups[[m]]
    if (length(voxel_tracks) == 0L) next

    # Focal-voxel row key "x_y_z" (from the first track's stored coords)
    focal_coords <- voxel_tracks[[1]]$coords
    if (is.null(focal_coords) || length(focal_coords) == 0L) next
    if (is.matrix(focal_coords)) focal_coords <- focal_coords[1, , drop = TRUE]
    row_key <- paste(focal_coords, collapse = "_")

    # Accumulators per target for this voxel
    hit_counts <- integer(length(cNode_labels))
    leave_no_target_count <- 0L

    # For each track in this voxel group, test intersection against each target
    for (n in seq_along(voxel_tracks)) {
      tr <- voxel_tracks[[n]]
      pts <- tr$tracks
      if (is.null(pts)) next
      pts <- as.matrix(pts)[, 1:3, drop = FALSE]
      storage.mode(pts) <- "integer"
      if (!nrow(pts)) next
      t_keys <- unique(key_of(pts))
      target_hit_any <- FALSE

      # Count once per track per target if any point hits that target ROI
      for (k in seq_along(cNode_labels)) {
        if (length(cnode_keys[[k]]) && any(t_keys %in% cnode_keys[[k]])) {
          hit_counts[k] <- hit_counts[k] + 1L
          target_hit_any <- TRUE
        }
      }

      local_flag <- tr$local
      if (is.null(local_flag) || is.na(local_flag[1])) local_flag <- 0L
      if (!target_hit_any && as.integer(local_flag[1]) == 0L) {
        leave_no_target_count <- leave_no_target_count + 1L
      }
    }

    # Update the corresponding row in coord_match_df
    idx <- which(coord_match_df$fNode_Coords == row_key)
    if (length(idx) == 1L) {
      coord_match_df[idx, cNode_labels] <- coord_match_df[idx, cNode_labels] + hit_counts
      coord_match_df$leave_no_target[idx] <- coord_match_df$leave_no_target[idx] + leave_no_target_count
    }
  }

  # Finalize 'local' and 'none' flags -----------------------------------------
  coord_match_df$local[is.na(coord_match_df$local)] <- 0L
  coord_match_df$leave_no_target[is.na(coord_match_df$leave_no_target)] <- 0L
  coord_match_df$none <- ifelse(is.na(coord_match_df$tracks) | coord_match_df$tracks == 0L, 1L, 0L)

  order_output(coord_match_df)
}
