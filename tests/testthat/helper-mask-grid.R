ftq_test_nifti_image <- function(arr,
                                 shift = c(0, 0, 0),
                                 scale = c(1, 1, 1),
                                 spatial_units = NULL) {
  scale <- rep_len(as.numeric(scale), 3L)
  shift <- rep_len(as.numeric(shift), 3L)
  image <- RNifti::asNifti(
    arr,
    reference = list(
      sform_code = 2L,
      qform_code = 0L,
      srow_x = c(scale[[1L]], 0, 0, shift[[1L]]),
      srow_y = c(0, scale[[2L]], 0, shift[[2L]]),
      srow_z = c(0, 0, scale[[3L]], shift[[3L]])
    )
  )
  if (!is.null(spatial_units)) RNifti::pixunits(image) <- spatial_units
  image
}

ftq_test_write_nifti <- function(path, arr, shift = c(0, 0, 0), datatype = "uint8") {
  RNifti::writeNifti(
    ftq_test_nifti_image(arr, shift = shift),
    path,
    datatype = datatype
  )
  invisible(path)
}
