# Convert all streamline points from world (mm) to voxel indices using a reference image.
# Requires RNifti; respects qform/sform.
# Returns a deep copy of 'trk' with coords replaced.
#' @importFrom RNifti worldToVoxel
ftq_trk_world_to_voxel <- function(trk, ref_image, add_one = FALSE) {
  if (is.null(trk) || is.null(trk$tracks)) stop("Invalid TRK object.")
  trk2 <- trk
  for (i in seq_along(trk2$tracks)) {
    pts_mm <- trk2$tracks[[i]]$coords
    if (!is.matrix(pts_mm) || ncol(pts_mm) < 3L) next
    vox <- ftq_world_mm_to_voxel(pts_mm, ref_image)
    if (add_one) vox <- vox + 1
    trk2$tracks[[i]]$coords <- vox
  }
  trk2
}

# Build a key vector "x,y,z" for integer coords
ftq_keys <- function(m) {
  if (!is.matrix(m) || ncol(m) < 3L) stop("coords matrix with 3 columns required")
  paste(m[,1], m[,2], m[,3], sep = ",")
}
