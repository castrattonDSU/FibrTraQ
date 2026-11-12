# TCK points are RAS+ world millimetres. Never reinterpret them by ROI overlap.
ftq_tck_validate_reference <- function(reference, context = "reference_image") {
  if (is.character(reference) && length(reference) == 1L &&
      !is.na(reference) && nzchar(reference)) {
    reference <- RNifti::readNifti(reference, internal = FALSE)
  }
  if (!inherits(reference, "niftiImage")) {
    stop("TCK coordinates require a NIfTI ", context,
         " with declared spatial geometry; a plain array is insufficient.",
         call. = FALSE)
  }
  dimensions <- dim(reference)
  if (length(dimensions) != 3L || any(!is.finite(dimensions)) ||
      any(dimensions <= 0)) {
    stop("TCK ", context, " must define a three-dimensional voxel grid.",
         call. = FALSE)
  }
  header <- RNifti::niftiHeader(reference)
  codes <- c(qform = header$qform_code, sform = header$sform_code)
  if (length(codes) != 2L || any(!is.finite(codes)) || any(codes < 0) ||
      !any(codes > 0)) {
    stop("TCK ", context, " must declare a valid NIfTI qform or sform. ",
         "Verify the source geometry; no affine will be guessed from ",
         "voxel size or streamline overlap.", call. = FALSE)
  }
  units <- ftq_local_search_unit_info(reference)
  if (!isTRUE(units$spatial_units_known)) {
    stop("TCK ", context, " must declare NIfTI spatial units (mm, m, or um). ",
         "TCK points are in millimetres; unknown image units are ambiguous.",
         call. = FALSE)
  }
  read_transform <- function(prefer_qform) {
    transform <- unclass(RNifti::xform(reference,
                                     useQuaternionFirst = prefer_qform))
    if (!is.matrix(transform) || !identical(dim(transform), c(4L, 4L)) ||
        any(!is.finite(transform)) ||
        !isTRUE(all.equal(as.numeric(transform[4L, ]), c(0, 0, 0, 1),
                          tolerance = 1e-8))) {
      stop("TCK ", context, " has an invalid NIfTI affine.", call. = FALSE)
    }
    transform[1:3, ] <- transform[1:3, , drop = FALSE] * units$scale_to_mm
    inverse <- tryCatch(solve(transform), error = function(e) NULL)
    if (is.null(inverse) || any(!is.finite(inverse))) {
      stop("TCK ", context, " has a singular or non-invertible NIfTI affine.",
           call. = FALSE)
    }
    transform
  }
  affine <- read_transform(codes[["qform"]] > 0)
  if (all(codes > 0)) {
    alternate <- read_transform(FALSE)
    # Compare voxel-centre positions in millimetres, including far corners.
    corners <- as.matrix(expand.grid(lapply(dimensions, function(n) c(0, n - 1))))
    disagreement <- max(abs(cbind(corners, 1) %*% t(affine - alternate)))
    if (!is.finite(disagreement) || disagreement > 1e-3) {
      stop("TCK ", context, " has conflicting qform and sform geometry ",
           "(voxel-centre disagreement exceeds 0.001 mm). Verify which ",
           "coordinate system matches the tractogram and resolve the image ",
           "metadata explicitly before analysis.", call. = FALSE)
    }
  }
  list(reference = reference, metadata = list(
    input_coordinates = "RAS+ world millimetres",
    voxel_coordinates = "one-based voxel centres",
    reference_spatial_units = units$spatial_units,
    reference_units_to_mm = units$scale_to_mm,
    qform_code = unname(codes[["qform"]]),
    sform_code = unname(codes[["sform"]]),
    selected_nifti_form = if (codes[["qform"]] > 0) "qform" else "sform",
    reference_affine_mm = affine,
    mapping_source = "tck_world_ras_mm_from_nifti_affine"
  ))
}

ftq_tck_coordinate_settings <- function(track_space = "auto", axis_flip = NULL,
                                        reference) {
  if (length(track_space) != 1L || is.na(track_space) ||
      !track_space %in% c("auto", "world")) {
    stop("TCK files contain RAS+ world coordinates in millimetres; ",
         "track_space must be 'auto' or 'world'.", call. = FALSE)
  }
  if (!is.null(axis_flip) &&
      (!is.logical(axis_flip) || length(axis_flip) != 3L ||
       anyNA(axis_flip) || any(axis_flip))) {
    stop("TCK files use their declared RAS+ world coordinates; axis_flip ",
         "must be NULL or c(FALSE, FALSE, FALSE). Verify input alignment ",
         "and apply any intended registration outside this analysis.",
         call. = FALSE)
  }
  validated <- ftq_tck_validate_reference(reference)
  c(list(track_space = "world", axis_flip = rep(FALSE, 3L)), validated)
}
