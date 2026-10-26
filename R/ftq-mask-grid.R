#' Report mask/template voxel-grid compatibility
#'
#' Inspect a NIfTI ROI mask against a template image without modifying files.
#' The report distinguishes dimension mismatch from affine/grid mismatch. When
#' reference masks are supplied, overlap and Jaccard metrics are computed only
#' for references that already share the template grid.
#'
#' This function is diagnostic only. It does not resample, rewrite, or
#' otherwise modify masks. Use it before a FibrTraQ run when a mask appears to
#' have the wrong dimensions or when two masks have matching dimensions but may
#' not occupy the same spatial grid.
#'
#' @param mask Path to a NIfTI mask.
#' @param template_image Path to the NIfTI image defining the intended grid.
#' @param reference_masks Optional named vector/list of same-grid masks used for
#'   overlap QC.
#' @param output_csv Optional CSV path for the report.
#'
#' @return A data.frame with grid, bounding-box, and optional overlap metrics.
#' @export
ftq_mask_grid_report <- function(mask,
                                 template_image,
                                 reference_masks = NULL,
                                 output_csv = NULL) {
  mask <- normalizePath(path.expand(mask), mustWork = TRUE)
  template_image <- normalizePath(path.expand(template_image), mustWork = TRUE)
  mask_img <- RNifti::readNifti(mask, internal = FALSE)
  template <- RNifti::readNifti(template_image, internal = FALSE)
  base <- ftq_mask_grid_base_row(mask_img, mask, template, template_image)

  reference_paths <- ftq_named_paths(reference_masks, prefix = "reference")
  if (!length(reference_paths)) {
    out <- cbind(base, ftq_empty_reference_grid_row(), stringsAsFactors = FALSE)
  } else {
    rows <- lapply(names(reference_paths), function(label) {
      ref_path <- normalizePath(path.expand(reference_paths[[label]]), mustWork = TRUE)
      ref_img <- RNifti::readNifti(ref_path, internal = FALSE)
      cbind(
        base,
        ftq_reference_grid_row(mask_img, template, ref_img, label, ref_path),
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
  }

  if (!is.null(output_csv)) {
    output_csv <- normalizePath(path.expand(output_csv), mustWork = FALSE)
    dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(out, output_csv, row.names = FALSE)
  }
  out
}

#' Make a same-grid binary union mask
#'
#' Create a biologically interpretable union mask from ROI masks that already
#' share the template grid. This is preferred over resampling when aligned
#' regional masks are available.
#'
#' Every input mask must match `template_image` in both dimensions and affine
#' spatial grid. The output is binary and is written with the template NIfTI
#' header. This is the recommended correction when a dataset provides several
#' trusted aligned regional masks whose union defines a local or whole-brain
#' analysis context.
#'
#' @param masks Named vector/list of NIfTI mask paths.
#' @param template_image NIfTI image defining the output grid and header.
#' @param output_path Output NIfTI path.
#' @param label Label recorded in returned metadata.
#' @param overwrite Whether to replace an existing output.
#'
#' @return A list containing `path`, `summary`, `inputs`, and
#'   `pairwise_overlaps`.
#' @export
ftq_make_union_mask <- function(masks,
                                template_image,
                                output_path,
                                label = "union",
                                overwrite = FALSE) {
  mask_paths <- ftq_named_paths(masks, prefix = "mask")
  if (!length(mask_paths)) stop("masks must contain one or more NIfTI mask paths.")
  template_image <- normalizePath(path.expand(template_image), mustWork = TRUE)
  output_path <- normalizePath(path.expand(output_path), mustWork = FALSE)
  if (file.exists(output_path) && !isTRUE(overwrite)) {
    stop("Union mask already exists. Use overwrite = TRUE to replace: ", output_path)
  }

  template <- RNifti::readNifti(template_image, internal = FALSE)
  template_dims <- dim(template)[1:3]
  union <- array(FALSE, dim = template_dims)
  mask_arrays <- list()
  input_rows <- list()

  for (nm in names(mask_paths)) {
    path <- normalizePath(path.expand(mask_paths[[nm]]), mustWork = TRUE)
    img <- RNifti::readNifti(path, internal = FALSE)
    ftq_assert_same_grid(
      img,
      template,
      context = paste0("Mask '", nm, "'"),
      template_context = "template_image"
    )
    arr <- as.array(img) != 0
    mask_arrays[[nm]] <- arr
    union <- union | arr
    input_rows[[nm]] <- data.frame(
      label = nm,
      path = path,
      nonzero_voxels = sum(arr, na.rm = TRUE),
      md5 = unname(tools::md5sum(path)),
      stringsAsFactors = FALSE
    )
  }

  pairwise <- ftq_pairwise_overlaps(mask_arrays)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  RNifti::writeNifti(array(as.integer(union), dim = template_dims), output_path,
                     template = template, datatype = "uint8")

  summary <- data.frame(
    label = as.character(label)[1],
    path = output_path,
    template_image = template_image,
    input_masks = length(mask_paths),
    union_voxels = sum(union, na.rm = TRUE),
    output_md5 = unname(tools::md5sum(output_path)),
    stringsAsFactors = FALSE
  )
  list(
    path = output_path,
    summary = summary,
    inputs = do.call(rbind, input_rows),
    pairwise_overlaps = pairwise
  )
}

#' Regrid a mask onto a template using trusted affine metadata
#'
#' This conservative fallback samples the source mask in world coordinates with
#' nearest-neighbor interpolation. It is appropriate only when the source
#' mask's NIfTI affine/header correctly describes its physical space. It should
#' not be used to "fix" masks whose spatial metadata are known to be wrong.
#'
#' `ftq_regrid_mask()` keeps binary/label-mask semantics by thresholding the
#' source and writing a binary output on the template grid. It refuses to
#' overwrite existing files unless `overwrite = TRUE` and emits QC warnings
#' when the output/source voxel-count ratio is suspicious or supplied reference
#' masks show poor overlap.
#'
#' @param mask Path to a NIfTI binary/label mask.
#' @param template_image Path to the output-grid template image.
#' @param output_path Output NIfTI path.
#' @param method Regridding method. Currently only `"nearest_world"` is
#'   supported.
#' @param threshold Source voxels greater than this value are included.
#' @param overwrite Whether to replace an existing output.
#' @param reference_masks Optional same-template-grid masks used for overlap QC.
#'
#' @return A list containing `path`, `summary`, and optional `reference_overlap`.
#' @export
ftq_regrid_mask <- function(mask,
                            template_image,
                            output_path,
                            method = "nearest_world",
                            threshold = 0,
                            overwrite = FALSE,
                            reference_masks = NULL) {
  method <- match.arg(method, "nearest_world")
  mask <- normalizePath(path.expand(mask), mustWork = TRUE)
  template_image <- normalizePath(path.expand(template_image), mustWork = TRUE)
  output_path <- normalizePath(path.expand(output_path), mustWork = FALSE)
  if (file.exists(output_path) && !isTRUE(overwrite)) {
    stop("Regridded mask already exists. Use overwrite = TRUE to replace: ", output_path)
  }

  src <- RNifti::readNifti(mask, internal = FALSE)
  template <- RNifti::readNifti(template_image, internal = FALSE)
  src_mask <- is.finite(as.array(src)) & as.array(src) > threshold
  out <- ftq_regrid_mask_nearest_world(src_mask, src, template)

  source_voxels <- sum(src_mask, na.rm = TRUE)
  output_voxels <- sum(out, na.rm = TRUE)
  ratio <- if (source_voxels > 0) output_voxels / source_voxels else NA_real_
  warnings <- character()
  if (is.finite(ratio) && (ratio < 0.5 || ratio > 2)) {
    warnings <- c(warnings, paste0(
      "Suspicious voxel-count change after regridding: output/source ratio = ",
      signif(ratio, 4), ". Confirm that the source affine/header is trusted."
    ))
  }

  reference_overlap <- ftq_regrid_reference_overlap(out, template, reference_masks)
  if (nrow(reference_overlap) && all(reference_overlap$jaccard < 0.2, na.rm = TRUE)) {
    warnings <- c(warnings,
                  "Regridded mask has low overlap with all supplied reference masks.")
  }
  for (msg in warnings) warning(msg, call. = FALSE)

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  RNifti::writeNifti(array(as.integer(out), dim = dim(out)), output_path,
                     template = template, datatype = "uint8")
  summary <- data.frame(
    path = output_path,
    source_mask = mask,
    template_image = template_image,
    method = method,
    threshold = threshold,
    source_nonzero_voxels = source_voxels,
    output_nonzero_voxels = output_voxels,
    output_source_voxel_ratio = ratio,
    warnings = paste(warnings, collapse = " | "),
    output_md5 = unname(tools::md5sum(output_path)),
    stringsAsFactors = FALSE
  )
  list(path = output_path, summary = summary, reference_overlap = reference_overlap)
}

ftq_named_paths <- function(x, prefix) {
  if (is.null(x)) return(character())
  if (is.list(x)) x <- unlist(x, use.names = TRUE)
  x <- as.character(x)
  if (!length(x)) return(character())
  nms <- names(x)
  if (is.null(nms) || any(!nzchar(nms))) {
    nms <- paste0(prefix, "_", seq_along(x))
  }
  stats::setNames(x, nms)
}

ftq_mask_grid_base_row <- function(mask_img, mask_path, template, template_path) {
  mask_info <- ftq_mask_grid_info(mask_img)
  template_info <- ftq_mask_grid_info(template)
  dims_match <- identical(mask_info$dims, template_info$dims)
  affine_match <- ftq_affine_equal(mask_img, template)
  cbind(
    data.frame(
      mask = mask_path,
      template_image = template_path,
      mask_dims = ftq_dim_label(mask_info$dims),
      template_dims = ftq_dim_label(template_info$dims),
      dims_match = dims_match,
      mask_voxel_size = paste(signif(mask_info$pixdim, 6), collapse = "x"),
      template_voxel_size = paste(signif(template_info$pixdim, 6), collapse = "x"),
      voxel_size_match = isTRUE(all.equal(mask_info$pixdim, template_info$pixdim, tolerance = 1e-6)),
      affine_match = affine_match,
      grid_match = dims_match && affine_match,
      nonzero_voxels = mask_info$nonzero_voxels,
      stringsAsFactors = FALSE
    ),
    as.data.frame(as.list(mask_info$voxel_bbox)),
    as.data.frame(as.list(mask_info$world_bbox)),
    stringsAsFactors = FALSE
  )
}

ftq_empty_reference_grid_row <- function() {
  data.frame(
    reference_label = NA_character_,
    reference_mask = NA_character_,
    reference_grid_match_template = NA,
    overlap_voxels = NA_integer_,
    union_voxels = NA_integer_,
    jaccard = NA_real_,
    stringsAsFactors = FALSE
  )
}

ftq_reference_grid_row <- function(mask_img, template, ref_img, label, ref_path) {
  ref_grid_match <- ftq_same_grid(ref_img, template)
  overlap <- union <- jaccard <- NA_real_
  if (ref_grid_match && ftq_same_grid(mask_img, template)) {
    a <- as.array(mask_img) != 0
    b <- as.array(ref_img) != 0
    overlap <- sum(a & b, na.rm = TRUE)
    union <- sum(a | b, na.rm = TRUE)
    jaccard <- if (union > 0) overlap / union else NA_real_
  }
  data.frame(
    reference_label = label,
    reference_mask = ref_path,
    reference_grid_match_template = ref_grid_match,
    overlap_voxels = as.integer(overlap),
    union_voxels = as.integer(union),
    jaccard = as.numeric(jaccard),
    stringsAsFactors = FALSE
  )
}

ftq_mask_grid_info <- function(img) {
  arr <- as.array(img)
  coords <- which(is.finite(arr) & arr != 0, arr.ind = TRUE)
  list(
    dims = as.integer(dim(img)[1:3]),
    pixdim = as.numeric(RNifti::pixdim(img)[1:3]),
    nonzero_voxels = nrow(coords),
    voxel_bbox = ftq_voxel_bbox(coords),
    world_bbox = ftq_world_bbox(coords, img)
  )
}

ftq_voxel_bbox <- function(coords) {
  nms <- c("voxel_bbox_x_min", "voxel_bbox_x_max",
           "voxel_bbox_y_min", "voxel_bbox_y_max",
           "voxel_bbox_z_min", "voxel_bbox_z_max")
  if (!nrow(coords)) return(stats::setNames(as.list(rep(NA_integer_, 6L)), nms))
  rng <- apply(coords[, 1:3, drop = FALSE], 2L, range)
  stats::setNames(as.list(as.integer(c(rng[, 1], rng[, 2], rng[, 3]))), nms)
}

ftq_world_bbox <- function(coords, img) {
  nms <- c("world_bbox_x_min", "world_bbox_x_max",
           "world_bbox_y_min", "world_bbox_y_max",
           "world_bbox_z_min", "world_bbox_z_max")
  if (!nrow(coords)) return(stats::setNames(as.list(rep(NA_real_, 6L)), nms))
  world <- RNifti::voxelToWorld(coords[, 1:3, drop = FALSE], img, simple = FALSE)
  world <- matrix(world, ncol = 3L)
  rng <- apply(world, 2L, range)
  stats::setNames(as.list(as.numeric(c(rng[, 1], rng[, 2], rng[, 3]))), nms)
}

ftq_same_grid <- function(a, b, tolerance = 1e-6) {
  identical(as.integer(dim(a)[1:3]), as.integer(dim(b)[1:3])) &&
    ftq_affine_equal(a, b, tolerance = tolerance)
}

ftq_affine_equal <- function(a, b, tolerance = 1e-6) {
  ftq_local_search_affine_equal(a, b, tolerance = tolerance)
}

ftq_assert_same_grid <- function(img, template, context = "mask", template_context = "template") {
  dims_match <- identical(as.integer(dim(img)[1:3]), as.integer(dim(template)[1:3]))
  affine_match <- ftq_affine_equal(img, template)
  if (dims_match && affine_match) return(invisible(TRUE))
  reason <- if (!dims_match) {
    paste0("dimensions differ: ", context, " = ", ftq_dim_label(dim(img)[1:3]),
           ", ", template_context, " = ", ftq_dim_label(dim(template)[1:3]))
  } else {
    paste0("dimensions match but affine/grid differs between ", context,
           " and ", template_context)
  }
  stop(reason, ". Use ftq_make_union_mask() when multiple aligned regional masks ",
       "are available; use ftq_regrid_mask() only when the source affine/header ",
       "is trusted.", call. = FALSE)
}

ftq_pairwise_overlaps <- function(mask_arrays) {
  if (length(mask_arrays) < 2L) {
    return(data.frame(mask_a = character(), mask_b = character(),
                      overlap_voxels = integer(), stringsAsFactors = FALSE))
  }
  nms <- names(mask_arrays)
  rows <- list()
  k <- 1L
  for (i in seq_len(length(mask_arrays) - 1L)) {
    for (j in seq.int(i + 1L, length(mask_arrays))) {
      rows[[k]] <- data.frame(
        mask_a = nms[[i]],
        mask_b = nms[[j]],
        overlap_voxels = sum(mask_arrays[[i]] & mask_arrays[[j]], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  do.call(rbind, rows)
}

ftq_regrid_mask_nearest_world <- function(src_mask, src_img, template) {
  target_dims <- as.integer(dim(template)[1:3])
  grid <- which(array(TRUE, dim = target_dims), arr.ind = TRUE)
  world_mm <- ftq_voxel_to_world_mm(grid, template)
  src_vox <- round(ftq_world_mm_to_voxel(world_mm, src_img))
  src_dims <- as.integer(dim(src_img)[1:3])
  in_bounds <- src_vox[, 1] >= 1L & src_vox[, 1] <= src_dims[1] &
    src_vox[, 2] >= 1L & src_vox[, 2] <= src_dims[2] &
    src_vox[, 3] >= 1L & src_vox[, 3] <= src_dims[3]
  out <- logical(nrow(grid))
  idx <- src_vox[in_bounds, , drop = FALSE]
  out[in_bounds] <- src_mask[cbind(idx[, 1], idx[, 2], idx[, 3])]
  array(out, dim = target_dims)
}

ftq_regrid_reference_overlap <- function(out, template, reference_masks) {
  refs <- ftq_named_paths(reference_masks, prefix = "reference")
  if (!length(refs)) {
    return(data.frame(reference_label = character(), reference_mask = character(),
                      overlap_voxels = integer(), union_voxels = integer(),
                      jaccard = numeric(), stringsAsFactors = FALSE))
  }
  rows <- lapply(names(refs), function(label) {
    path <- normalizePath(path.expand(refs[[label]]), mustWork = TRUE)
    ref <- RNifti::readNifti(path, internal = FALSE)
    ftq_assert_same_grid(ref, template, context = paste0("reference mask '", label, "'"),
                         template_context = "template_image")
    ref_arr <- as.array(ref) != 0
    overlap <- sum(out & ref_arr, na.rm = TRUE)
    union <- sum(out | ref_arr, na.rm = TRUE)
    data.frame(
      reference_label = label,
      reference_mask = path,
      overlap_voxels = overlap,
      union_voxels = union,
      jaccard = if (union > 0) overlap / union else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
