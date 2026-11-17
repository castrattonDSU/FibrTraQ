# Build the archived native routines only when a research audit first uses one.
ftq_native_audits <- local({
  compiled <- NULL
  function() {
    if (is.null(compiled)) {
      routines <- new.env(parent = baseenv())
      Rcpp::sourceCpp(file.path(.archive_root, "src", "native_audits.cpp"),
                      env = routines)
      compiled <<- routines
    }
    compiled
  }
})

ftq_spacing_chunk_cpp <- function(tracks, raw_to_voxel, voxel_size,
                                  bin_width_mm = 0.01, bin_width_vox = 0.01) {
  ftq_native_audits()$ftq_spacing_chunk_cpp(
    tracks, raw_to_voxel, voxel_size, bin_width_mm, bin_width_vox
  )
}

ftq_endpoint_audit_chunk_cpp <- function(tracks, raw_to_voxel, method, axis_flip,
                                         vol_dim, region_lut, region_lut_dilated,
                                         context_lut, streamline_voxelization_mode,
                                         penetration_hist_max = 500L) {
  ftq_native_audits()$ftq_endpoint_audit_chunk_cpp(
    tracks, raw_to_voxel, method, axis_flip, vol_dim, region_lut,
    region_lut_dilated, context_lut, streamline_voxelization_mode, penetration_hist_max
  )
}
