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
#' @noRd
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
  stop(reason, ". Supply images aligned to the same voxel grid and verify ",
       "their alignment against the subject's anatomy.", call. = FALSE)
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
