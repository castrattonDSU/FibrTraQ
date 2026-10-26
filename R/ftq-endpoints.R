#' Retrieve endpoint-analysis tables from a FibrTraQ result
#'
#' Endpoint connectivity is an additive analysis that distinguishes full-path
#' target intersections from streamlines whose opposite endpoints lie in the
#' focal and target ROIs. The returned tables also contain target-penetration
#' and endpoint anatomical-context summaries.
#'
#' @param result A data frame returned directly by [FibrTraQ()]. Endpoint
#'   attributes are not retained by a CSV round trip; use the dedicated
#'   endpoint CSV files for persisted results.
#'
#' @return A list with `summary`, `penetration`, `context`, and `metadata`.
#' @export
ftq_endpoint_tables <- function(result) {
  if (!is.data.frame(result)) {
    stop("result must be a data frame returned directly by FibrTraQ().")
  }
  tables <- list(
    summary = attr(result, "fibrtraq_endpoint_summary"),
    penetration = attr(result, "fibrtraq_endpoint_penetration"),
    context = attr(result, "fibrtraq_endpoint_context"),
    metadata = attr(result, "fibrtraq_endpoint_metadata")
  )
  if (all(vapply(tables, is.null, logical(1)))) {
    stop("Endpoint analysis is unavailable for this result.")
  }
  tables
}

ftq_endpoint_prepare_contexts <- function(context_list,
                                          context_labels,
                                          vol_dim,
                                          reference = NULL) {
  if (is.null(context_list)) context_list <- list()
  if (!is.list(context_list)) {
    stop("endpoint_context_list must be NULL or a list of binary 3-D masks.")
  }
  if (!length(context_list)) {
    return(list(
      lut = matrix(FALSE, nrow = prod(vol_dim), ncol = 0L),
      labels = character(),
      metadata = data.frame(
        context = character(), context_type = character(), source_path = character(),
        checksum_md5 = character(), mask_voxels = numeric(), dimensions = character(),
        affine_match_reference = logical(), stringsAsFactors = FALSE
      )
    ))
  }

  if (is.null(context_labels) || length(context_labels) != length(context_list)) {
    stop("endpoint_context_labels must contain one label per endpoint context mask.")
  }
  context_labels <- as.character(context_labels)
  if (anyNA(context_labels) || any(!nzchar(context_labels)) || anyDuplicated(context_labels)) {
    stop("endpoint_context_labels must be non-empty, non-missing, and unique.")
  }
  if ("reference_grid" %in% context_labels) {
    stop("endpoint_context_labels cannot use the reserved label 'reference_grid'.")
  }

  lut <- matrix(FALSE, nrow = prod(vol_dim), ncol = length(context_list))
  metadata <- vector("list", length(context_list))
  for (i in seq_along(context_list)) {
    mask <- context_list[[i]]
    mask_dim <- dim(mask)
    if (is.null(mask_dim) || length(mask_dim) != 3L ||
        !identical(as.integer(mask_dim), as.integer(vol_dim))) {
      stop(
        "Endpoint context '", context_labels[[i]],
        "' must be a 3-D mask with dimensions ", paste(vol_dim, collapse = "x"), "."
      )
    }
    values <- as.numeric(mask)
    if (any(!is.finite(values)) || any(!values %in% c(0, 1))) {
      stop("Endpoint context '", context_labels[[i]], "' must contain only finite binary 0/1 values.")
    }
    lut[, i] <- values != 0

    affine_match <- NA
    if (!is.null(reference)) {
      affine_match <- tryCatch(ftq_affine_equal(mask, reference), error = function(e) NA)
      if (identical(affine_match, FALSE)) {
        stop("Endpoint context '", context_labels[[i]], "' affine does not match the focal reference.")
      }
    }
    source_path <- attr(mask, "fibrtraq_source_path", exact = TRUE)
    checksum <- attr(mask, "fibrtraq_checksum_md5", exact = TRUE)
    metadata[[i]] <- data.frame(
      context = context_labels[[i]],
      context_type = "binary_mask",
      source_path = if (is.null(source_path)) NA_character_ else as.character(source_path)[1],
      checksum_md5 = if (is.null(checksum)) NA_character_ else as.character(checksum)[1],
      mask_voxels = sum(lut[, i]),
      dimensions = paste(vol_dim, collapse = "x"),
      affine_match_reference = if (length(affine_match)) as.logical(affine_match)[1] else NA,
      stringsAsFactors = FALSE
    )
  }
  colnames(lut) <- context_labels
  list(lut = lut, labels = context_labels, metadata = do.call(rbind, metadata))
}

ftq_endpoint_dilate_lut <- function(lut, vol_dim) {
  if (length(lut) != prod(vol_dim)) stop("Endpoint LUT dimensions are inconsistent.")
  lut <- as.logical(lut)
  if (is.function(get0("ftq_dilate_lut_cpp", mode = "function"))) {
    return(as.logical(ftq_dilate_lut_cpp(lut, as.integer(vol_dim))))
  }
  active <- which(lut)
  if (!length(active)) return(lut)
  coordinates <- arrayInd(active, .dim = as.integer(vol_dim))
  out <- lut
  for (dx in -1:1) {
    for (dy in -1:1) {
      for (dz in -1:1) {
        shifted <- sweep(coordinates, 2L, c(dx, dy, dz), `+`)
        inside <- shifted[, 1] >= 1L & shifted[, 1] <= vol_dim[[1]] &
          shifted[, 2] >= 1L & shifted[, 2] <= vol_dim[[2]] &
          shifted[, 3] >= 1L & shifted[, 3] <= vol_dim[[3]]
        if (!any(inside)) next
        shifted <- shifted[inside, , drop = FALSE]
        keys <- shifted[, 1] +
          (shifted[, 2] - 1L) * vol_dim[[1]] +
          (shifted[, 3] - 1L) * vol_dim[[1]] * vol_dim[[2]]
        out[keys] <- TRUE
      }
    }
  }
  out
}

ftq_endpoint_accumulator <- function(n_targets, n_contexts) {
  classes <- c(
    "shared_mask_ambiguous",
    "paired_focal_target_endpoints",
    "target_endpoint_only",
    "focal_endpoint_only",
    "neither_relevant_endpoint"
  )
  acc <- new.env(parent = emptyenv())
  acc$valid_streamlines <- 0
  acc$endpoint_count <- 0
  acc$endpoint_out_of_bounds <- 0
  acc$endpoint_nonfinite <- 0
  acc$endpoint_outside_one_voxel <- 0
  acc$endpoint_outside_more_than_one_voxel <- 0
  acc$endpoint_out_of_bounds_low <- numeric(3L)
  acc$endpoint_out_of_bounds_high <- numeric(3L)
  acc$endpoint_out_of_bounds_low_distance <- matrix(0, nrow = 3L, ncol = 2L)
  acc$endpoint_out_of_bounds_high_distance <- matrix(0, nrow = 3L, ncol = 2L)
  acc$grid_any_endpoint_inside <- 0
  acc$grid_both_endpoints_inside <- 0
  acc$context_endpoint_inside <- numeric(n_contexts)
  acc$context_any_endpoint_inside <- numeric(n_contexts)
  acc$context_both_endpoints_inside <- numeric(n_contexts)
  acc$focal_endpoint_streamlines_exact <- 0
  acc$focal_endpoint_streamlines_dilated1 <- 0
  acc$paired_endpoint_streamlines_exact <- numeric(n_targets)
  acc$paired_endpoint_streamlines_dilated1 <- numeric(n_targets)
  acc$paired_intersection_streamlines_exact <- numeric(n_targets)
  acc$paired_intersection_streamlines_dilated1 <- numeric(n_targets)
  acc$focal_endpoint_target_intersection_exact <- numeric(n_targets)
  acc$focal_endpoint_target_intersection_dilated1 <- numeric(n_targets)
  acc$class_exact <- matrix(0, nrow = n_targets, ncol = length(classes),
                            dimnames = list(NULL, classes))
  acc$class_dilated1 <- matrix(0, nrow = n_targets, ncol = length(classes),
                               dimnames = list(NULL, classes))
  acc$penetration <- replicate(n_targets, numeric(), simplify = FALSE)
  acc
}

ftq_endpoint_discretize <- function(coords, method, axis_flip, vol_dim) {
  discrete <- switch(
    method,
    round = round(coords),
    floor = floor(coords),
    ceil = ceiling(coords),
    stop("Unsupported endpoint discretization method: ", method)
  )
  for (axis in seq_len(3L)) {
    if (isTRUE(axis_flip[[axis]])) {
      discrete[, axis] <- vol_dim[[axis]] - discrete[, axis] + 1
    }
  }
  discrete
}

ftq_endpoint_record_context_r <- function(acc,
                                          mapped,
                                          method,
                                          axis_flip,
                                          vol_dim,
                                          context_lut) {
  if (!nrow(mapped)) return(invisible(NULL))
  endpoint_rows <- if (nrow(mapped) == 1L) 1L else c(1L, nrow(mapped))
  endpoints <- mapped[endpoint_rows, 1:3, drop = FALSE]
  discrete <- ftq_endpoint_discretize(endpoints, method, axis_flip, vol_dim)
  finite <- apply(discrete, 1L, function(x) all(is.finite(x)))
  in_bounds <- finite &
    discrete[, 1] >= 1 & discrete[, 1] <= vol_dim[1] &
    discrete[, 2] >= 1 & discrete[, 2] <= vol_dim[2] &
    discrete[, 3] >= 1 & discrete[, 3] <= vol_dim[3]

  acc$valid_streamlines <- acc$valid_streamlines + 1
  acc$endpoint_count <- acc$endpoint_count + nrow(discrete)
  acc$endpoint_out_of_bounds <- acc$endpoint_out_of_bounds + sum(!in_bounds)
  acc$endpoint_nonfinite <- acc$endpoint_nonfinite + sum(!finite)
  acc$grid_any_endpoint_inside <- acc$grid_any_endpoint_inside + as.numeric(any(in_bounds))
  acc$grid_both_endpoints_inside <- acc$grid_both_endpoints_inside + as.numeric(all(in_bounds))

  for (row in seq_len(nrow(discrete))) {
    if (!finite[[row]]) next
    max_distance <- 0
    for (axis in seq_len(3L)) {
      coordinate <- discrete[row, axis]
      if (coordinate < 1) {
        distance <- 1 - coordinate
        acc$endpoint_out_of_bounds_low[axis] <- acc$endpoint_out_of_bounds_low[axis] + 1
        bin <- if (distance == 1) 1L else 2L
        acc$endpoint_out_of_bounds_low_distance[axis, bin] <-
          acc$endpoint_out_of_bounds_low_distance[axis, bin] + 1
        max_distance <- max(max_distance, distance)
      } else if (coordinate > vol_dim[[axis]]) {
        distance <- coordinate - vol_dim[[axis]]
        acc$endpoint_out_of_bounds_high[axis] <- acc$endpoint_out_of_bounds_high[axis] + 1
        bin <- if (distance == 1) 1L else 2L
        acc$endpoint_out_of_bounds_high_distance[axis, bin] <-
          acc$endpoint_out_of_bounds_high_distance[axis, bin] + 1
        max_distance <- max(max_distance, distance)
      }
    }
    if (max_distance == 1) {
      acc$endpoint_outside_one_voxel <- acc$endpoint_outside_one_voxel + 1
    } else if (max_distance > 1) {
      acc$endpoint_outside_more_than_one_voxel <-
        acc$endpoint_outside_more_than_one_voxel + 1
    }
  }

  keys <- rep(NA_integer_, nrow(discrete))
  if (any(in_bounds)) {
    coords_in <- discrete[in_bounds, , drop = FALSE]
    keys[in_bounds] <- as.integer(
      coords_in[, 1] + (coords_in[, 2] - 1) * vol_dim[1] +
        (coords_in[, 3] - 1) * vol_dim[1] * vol_dim[2]
    )
  }
  if (ncol(context_lut)) {
    for (context in seq_len(ncol(context_lut))) {
      inside <- rep(FALSE, length(keys))
      valid <- !is.na(keys)
      inside[valid] <- context_lut[keys[valid], context]
      acc$context_endpoint_inside[context] <-
        acc$context_endpoint_inside[context] + sum(inside)
      acc$context_any_endpoint_inside[context] <-
        acc$context_any_endpoint_inside[context] + as.numeric(any(inside))
      acc$context_both_endpoints_inside[context] <-
        acc$context_both_endpoints_inside[context] + as.numeric(all(inside))
    }
  }
  invisible(NULL)
}

ftq_endpoint_membership <- function(key, lut) {
  !is.na(key) && key >= 1L && key <= length(lut) && isTRUE(lut[[key]])
}

ftq_endpoint_classification <- function(start_focal,
                                        end_focal,
                                        start_target,
                                        end_target,
                                        two_endpoints) {
  focal_any <- start_focal || end_focal
  target_any <- start_target || end_target
  shared <- (start_focal && start_target) || (end_focal && end_target)
  paired <- isTRUE(two_endpoints) &&
    ((start_focal && end_target) || (end_focal && start_target))
  class <- if (shared) {
    1L
  } else if (paired) {
    2L
  } else if (target_any) {
    3L
  } else if (focal_any) {
    4L
  } else {
    5L
  }
  list(class = class, paired = paired, focal_any = focal_any, target_any = target_any)
}

ftq_endpoint_record_focal_r <- function(acc,
                                        chosen_keys,
                                        two_endpoints,
                                        target_hit,
                                        lut_f,
                                        lut_f_dilated,
                                        target_lut,
                                        target_lut_dilated) {
  if (!length(chosen_keys)) return(invisible(NULL))
  start_key <- chosen_keys[[1L]]
  end_key <- if (isTRUE(two_endpoints)) chosen_keys[[length(chosen_keys)]] else start_key
  start_focal_exact <- ftq_endpoint_membership(start_key, lut_f)
  end_focal_exact <- ftq_endpoint_membership(end_key, lut_f)
  start_focal_d1 <- ftq_endpoint_membership(start_key, lut_f_dilated)
  end_focal_d1 <- ftq_endpoint_membership(end_key, lut_f_dilated)
  focal_any_exact <- start_focal_exact || end_focal_exact
  focal_any_d1 <- start_focal_d1 || end_focal_d1
  acc$focal_endpoint_streamlines_exact <-
    acc$focal_endpoint_streamlines_exact + as.numeric(focal_any_exact)
  acc$focal_endpoint_streamlines_dilated1 <-
    acc$focal_endpoint_streamlines_dilated1 + as.numeric(focal_any_d1)

  valid_keys <- unique(chosen_keys[!is.na(chosen_keys)])
  for (target in seq_len(ncol(target_lut))) {
    exact <- ftq_endpoint_classification(
      start_focal_exact,
      end_focal_exact,
      ftq_endpoint_membership(start_key, target_lut[, target]),
      ftq_endpoint_membership(end_key, target_lut[, target]),
      two_endpoints
    )
    dilated1 <- ftq_endpoint_classification(
      start_focal_d1,
      end_focal_d1,
      ftq_endpoint_membership(start_key, target_lut_dilated[, target]),
      ftq_endpoint_membership(end_key, target_lut_dilated[, target]),
      two_endpoints
    )
    acc$paired_endpoint_streamlines_exact[target] <-
      acc$paired_endpoint_streamlines_exact[target] + as.numeric(exact$paired)
    acc$paired_endpoint_streamlines_dilated1[target] <-
      acc$paired_endpoint_streamlines_dilated1[target] + as.numeric(dilated1$paired)
    if (!isTRUE(target_hit[[target]])) next

    acc$class_exact[target, exact$class] <- acc$class_exact[target, exact$class] + 1
    acc$class_dilated1[target, dilated1$class] <-
      acc$class_dilated1[target, dilated1$class] + 1
    acc$paired_intersection_streamlines_exact[target] <-
      acc$paired_intersection_streamlines_exact[target] + as.numeric(exact$paired)
    acc$paired_intersection_streamlines_dilated1[target] <-
      acc$paired_intersection_streamlines_dilated1[target] + as.numeric(dilated1$paired)
    acc$focal_endpoint_target_intersection_exact[target] <-
      acc$focal_endpoint_target_intersection_exact[target] + as.numeric(exact$focal_any)
    acc$focal_endpoint_target_intersection_dilated1[target] <-
      acc$focal_endpoint_target_intersection_dilated1[target] + as.numeric(dilated1$focal_any)

    penetration <- if (length(valid_keys)) sum(target_lut[valid_keys, target]) else 0
    key <- as.character(penetration)
    hist <- acc$penetration[[target]]
    if (!length(hist)) {
      hist <- stats::setNames(1, key)
    } else if (key %in% names(hist)) {
      hist[[key]] <- hist[[key]] + 1
    } else {
      hist[[key]] <- 1
    }
    acc$penetration[[target]] <- hist
  }
  invisible(NULL)
}

ftq_endpoint_merge_cpp <- function(acc, endpoint) {
  if (is.null(endpoint)) return(invisible(NULL))
  scalar_fields <- c(
    "valid_streamlines", "endpoint_count", "endpoint_out_of_bounds",
    "endpoint_nonfinite", "endpoint_outside_one_voxel",
    "endpoint_outside_more_than_one_voxel", "grid_any_endpoint_inside",
    "grid_both_endpoints_inside", "focal_endpoint_streamlines_exact",
    "focal_endpoint_streamlines_dilated1"
  )
  vector_fields <- c(
    "endpoint_out_of_bounds_low", "endpoint_out_of_bounds_high",
    "context_endpoint_inside", "context_any_endpoint_inside",
    "context_both_endpoints_inside", "paired_endpoint_streamlines_exact",
    "paired_endpoint_streamlines_dilated1",
    "paired_intersection_streamlines_exact",
    "paired_intersection_streamlines_dilated1",
    "focal_endpoint_target_intersection_exact",
    "focal_endpoint_target_intersection_dilated1"
  )
  matrix_fields <- c(
    "endpoint_out_of_bounds_low_distance",
    "endpoint_out_of_bounds_high_distance",
    "class_exact", "class_dilated1"
  )
  for (field in scalar_fields) acc[[field]] <- acc[[field]] + as.numeric(endpoint[[field]])
  for (field in vector_fields) acc[[field]] <- acc[[field]] + as.numeric(endpoint[[field]])
  for (field in matrix_fields) acc[[field]] <- acc[[field]] + as.matrix(endpoint[[field]])

  penetration <- endpoint$penetration
  if (!is.null(penetration) && nrow(penetration)) {
    for (row in seq_len(nrow(penetration))) {
      target <- as.integer(penetration$target_index[[row]])
      key <- as.character(penetration$unique_target_voxels[[row]])
      count <- as.numeric(penetration$streamline_count[[row]])
      hist <- acc$penetration[[target]]
      if (!length(hist)) {
        hist <- stats::setNames(count, key)
      } else if (key %in% names(hist)) {
        hist[[key]] <- hist[[key]] + count
      } else {
        hist[[key]] <- count
      }
      acc$penetration[[target]] <- hist
    }
  }
  invisible(NULL)
}

ftq_endpoint_percent <- function(numerator, denominator) {
  out <- 100 * numerator / denominator
  invalid <- !is.finite(denominator) | denominator <= 0
  if (length(invalid) == 1L) {
    if (invalid) out[] <- NA_real_
  } else {
    out[invalid] <- NA_real_
  }
  out
}

ftq_endpoint_voxel_size <- function(fNode, ref_image = NULL) {
  stored <- attr(fNode, "fibrtraq_voxel_size_mm", exact = TRUE)
  if (!is.null(stored) && length(stored) >= 3L) return(as.numeric(stored[1:3]))
  source <- if (!is.null(ref_image)) ref_image else fNode
  size <- tryCatch(as.numeric(RNifti::pixdim(source)[1:3]), error = function(e) rep(NA_real_, 3L))
  size <- size * ftq_local_search_unit_info(source)$scale_to_mm
  if (length(size) < 3L || any(!is.finite(size)) || any(size <= 0)) rep(NA_real_, 3L) else size
}

ftq_endpoint_finalize <- function(acc,
                                  target_labels,
                                  target_unique_counts,
                                  unique_focal_streamlines,
                                  lut_f,
                                  lut_f_dilated,
                                  target_lut,
                                  target_lut_dilated,
                                  context_info,
                                  vol_dim,
                                  voxel_size,
                                  methods,
                                  streamline_voxelization_mode,
                                  track_space,
                                  axis_flip) {
  classes <- colnames(acc$class_exact)
  tolerance_data <- list(
    exact = list(
      classes = acc$class_exact,
      focal_endpoints = acc$focal_endpoint_streamlines_exact,
      paired = acc$paired_endpoint_streamlines_exact,
      paired_intersection = acc$paired_intersection_streamlines_exact,
      focal_endpoint_target_intersection = acc$focal_endpoint_target_intersection_exact,
      focal_lut = lut_f,
      target_lut = target_lut,
      definition = "endpoint voxel lies inside the supplied ROI mask"
    ),
    dilated1 = list(
      classes = acc$class_dilated1,
      focal_endpoints = acc$focal_endpoint_streamlines_dilated1,
      paired = acc$paired_endpoint_streamlines_dilated1,
      paired_intersection = acc$paired_intersection_streamlines_dilated1,
      focal_endpoint_target_intersection = acc$focal_endpoint_target_intersection_dilated1,
      focal_lut = lut_f_dilated,
      target_lut = target_lut_dilated,
      definition = "voxel-grid sensitivity: endpoint lies within Chebyshev radius 1 (26-neighbor dilation)"
    )
  )

  summary_rows <- list()
  row_index <- 0L
  for (tolerance in names(tolerance_data)) {
    values <- tolerance_data[[tolerance]]
    for (target in seq_along(target_labels)) {
      row_index <- row_index + 1L
      class_counts <- stats::setNames(as.numeric(values$classes[target, ]), classes)
      intersections <- as.numeric(target_unique_counts[[target]])
      focal_endpoints <- as.numeric(values$focal_endpoints)
      paired <- as.numeric(values$paired[[target]])
      paired_intersection <- as.numeric(values$paired_intersection[[target]])
      focal_endpoint_intersection <-
        as.numeric(values$focal_endpoint_target_intersection[[target]])
      target_endpoint_intersection <- class_counts[["shared_mask_ambiguous"]] +
        class_counts[["paired_focal_target_endpoints"]] +
        class_counts[["target_endpoint_only"]]
      focal_endpoint_intersection_from_classes <-
        class_counts[["shared_mask_ambiguous"]] +
        class_counts[["paired_focal_target_endpoints"]] +
        class_counts[["focal_endpoint_only"]]
      overlap <- sum(values$focal_lut & values$target_lut[, target])
      focal_voxels <- sum(values$focal_lut)
      target_voxels <- sum(values$target_lut[, target])
      smaller <- min(focal_voxels, target_voxels)

      summary_rows[[row_index]] <- data.frame(
        target = target_labels[[target]],
        tolerance = tolerance,
        tolerance_definition = values$definition,
        unique_focal_intersecting_streamlines = as.numeric(unique_focal_streamlines),
        target_intersection_streamlines = intersections,
        target_intersection_percent_of_unique_focal_intersecting_streamlines =
          ftq_endpoint_percent(intersections, unique_focal_streamlines),
        focal_endpoint_streamlines = focal_endpoints,
        paired_endpoint_streamlines = paired,
        target_intersection_with_paired_endpoints = paired_intersection,
        paired_percent_of_target_intersections =
          ftq_endpoint_percent(paired_intersection, intersections),
        paired_percent_of_focal_endpoint_streamlines =
          ftq_endpoint_percent(paired, focal_endpoints),
        target_intersection_with_any_target_endpoint = target_endpoint_intersection,
        target_endpoint_percent_of_target_intersections =
          ftq_endpoint_percent(target_endpoint_intersection, intersections),
        target_intersection_with_any_focal_endpoint =
          focal_endpoint_intersection_from_classes,
        focal_endpoint_percent_of_target_intersections =
          ftq_endpoint_percent(focal_endpoint_intersection_from_classes, intersections),
        focal_endpoint_with_target_intersection_streamlines = focal_endpoint_intersection,
        target_intersection_percent_of_focal_endpoint_streamlines =
          ftq_endpoint_percent(focal_endpoint_intersection, focal_endpoints),
        shared_mask_ambiguous = class_counts[["shared_mask_ambiguous"]],
        paired_focal_target_endpoints = class_counts[["paired_focal_target_endpoints"]],
        target_endpoint_only = class_counts[["target_endpoint_only"]],
        focal_endpoint_only = class_counts[["focal_endpoint_only"]],
        neither_relevant_endpoint = class_counts[["neither_relevant_endpoint"]],
        shared_mask_ambiguous_percent_of_target_intersections =
          ftq_endpoint_percent(class_counts[["shared_mask_ambiguous"]], intersections),
        target_endpoint_only_percent_of_target_intersections =
          ftq_endpoint_percent(class_counts[["target_endpoint_only"]], intersections),
        focal_endpoint_only_percent_of_target_intersections =
          ftq_endpoint_percent(class_counts[["focal_endpoint_only"]], intersections),
        neither_relevant_endpoint_percent_of_target_intersections =
          ftq_endpoint_percent(class_counts[["neither_relevant_endpoint"]], intersections),
        focal_mask_voxels = focal_voxels,
        target_mask_voxels = target_voxels,
        focal_target_overlap_voxels = overlap,
        focal_target_overlap_percent_of_smaller_mask =
          ftq_endpoint_percent(overlap, smaller),
        voxel_size_x_mm = voxel_size[[1]],
        voxel_size_y_mm = voxel_size[[2]],
        voxel_size_z_mm = voxel_size[[3]],
        streamline_voxelization_mode = streamline_voxelization_mode,
        target_intersection_mode = "full_streamline",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  summary <- if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame()
  rownames(summary) <- NULL

  penetration_rows <- list()
  row_index <- 0L
  for (target in seq_along(target_labels)) {
    hist <- acc$penetration[[target]]
    if (!length(hist)) next
    bins <- as.integer(names(hist))
    order_index <- order(bins)
    bins <- bins[order_index]
    counts <- as.numeric(hist[order_index])
    denominator <- as.numeric(target_unique_counts[[target]])
    row_index <- row_index + 1L
    penetration_rows[[row_index]] <- data.frame(
      target = target_labels[[target]],
      unique_target_voxels = bins,
      streamline_count = counts,
      target_intersection_streamlines = denominator,
      percent_of_target_intersections = ftq_endpoint_percent(counts, denominator),
      streamline_voxelization_mode = streamline_voxelization_mode,
      penetration_interpretation = if (identical(streamline_voxelization_mode, "sample_points")) {
        "stored-point target penetration; sparse streamline points can undercount traversed voxels"
      } else {
        "segment-interpolated target penetration"
      },
      stringsAsFactors = FALSE
    )
  }
  penetration <- if (length(penetration_rows)) {
    do.call(rbind, penetration_rows)
  } else {
    data.frame(
      target = character(), unique_target_voxels = integer(), streamline_count = numeric(),
      target_intersection_streamlines = numeric(), percent_of_target_intersections = numeric(),
      streamline_voxelization_mode = character(), penetration_interpretation = character(),
      stringsAsFactors = FALSE
    )
  }
  rownames(penetration) <- NULL

  context_base <- data.frame(
    context = "reference_grid",
    context_type = "reference_grid",
    valid_streamlines = acc$valid_streamlines,
    endpoint_count = acc$endpoint_count,
    endpoints_inside = acc$endpoint_count - acc$endpoint_out_of_bounds,
    endpoints_outside = acc$endpoint_out_of_bounds,
    endpoints_inside_percent = ftq_endpoint_percent(
      acc$endpoint_count - acc$endpoint_out_of_bounds, acc$endpoint_count
    ),
    endpoints_nonfinite = acc$endpoint_nonfinite,
    endpoints_one_voxel_outside = acc$endpoint_outside_one_voxel,
    endpoints_more_than_one_voxel_outside = acc$endpoint_outside_more_than_one_voxel,
    streamlines_with_any_endpoint_inside = acc$grid_any_endpoint_inside,
    streamlines_with_both_endpoints_inside = acc$grid_both_endpoints_inside,
    streamlines_with_both_endpoints_inside_percent =
      ftq_endpoint_percent(acc$grid_both_endpoints_inside, acc$valid_streamlines),
    x_below_grid = acc$endpoint_out_of_bounds_low[[1]],
    x_above_grid = acc$endpoint_out_of_bounds_high[[1]],
    y_below_grid = acc$endpoint_out_of_bounds_low[[2]],
    y_above_grid = acc$endpoint_out_of_bounds_high[[2]],
    z_below_grid = acc$endpoint_out_of_bounds_low[[3]],
    z_above_grid = acc$endpoint_out_of_bounds_high[[3]],
    x_one_voxel_below = acc$endpoint_out_of_bounds_low_distance[1, 1],
    x_more_than_one_voxel_below = acc$endpoint_out_of_bounds_low_distance[1, 2],
    x_one_voxel_above = acc$endpoint_out_of_bounds_high_distance[1, 1],
    x_more_than_one_voxel_above = acc$endpoint_out_of_bounds_high_distance[1, 2],
    y_one_voxel_below = acc$endpoint_out_of_bounds_low_distance[2, 1],
    y_more_than_one_voxel_below = acc$endpoint_out_of_bounds_low_distance[2, 2],
    y_one_voxel_above = acc$endpoint_out_of_bounds_high_distance[2, 1],
    y_more_than_one_voxel_above = acc$endpoint_out_of_bounds_high_distance[2, 2],
    z_one_voxel_below = acc$endpoint_out_of_bounds_low_distance[3, 1],
    z_more_than_one_voxel_below = acc$endpoint_out_of_bounds_low_distance[3, 2],
    z_one_voxel_above = acc$endpoint_out_of_bounds_high_distance[3, 1],
    z_more_than_one_voxel_above = acc$endpoint_out_of_bounds_high_distance[3, 2],
    source_path = NA_character_,
    checksum_md5 = NA_character_,
    mask_voxels = prod(vol_dim),
    dimensions = paste(vol_dim, collapse = "x"),
    affine_match_reference = TRUE,
    stringsAsFactors = FALSE
  )
  context <- context_base
  if (length(context_info$labels)) {
    mask_rows <- lapply(seq_along(context_info$labels), function(i) {
      inside <- acc$context_endpoint_inside[[i]]
      metadata_row <- context_info$metadata[i, , drop = FALSE]
      row <- context_base
      row$context <- context_info$labels[[i]]
      row$context_type <- "binary_mask"
      row$endpoints_inside <- inside
      row$endpoints_outside <- acc$endpoint_count - inside
      row$endpoints_inside_percent <- ftq_endpoint_percent(inside, acc$endpoint_count)
      row$streamlines_with_any_endpoint_inside <- acc$context_any_endpoint_inside[[i]]
      row$streamlines_with_both_endpoints_inside <- acc$context_both_endpoints_inside[[i]]
      row$streamlines_with_both_endpoints_inside_percent <- ftq_endpoint_percent(
        acc$context_both_endpoints_inside[[i]], acc$valid_streamlines
      )
      boundary_columns <- c(
        "endpoints_nonfinite", "endpoints_one_voxel_outside",
        "endpoints_more_than_one_voxel_outside", "x_below_grid", "x_above_grid",
        "y_below_grid", "y_above_grid", "z_below_grid", "z_above_grid",
        grep("_voxel_(below|above)$", names(row), value = TRUE),
        grep("more_than_one_voxel_(below|above)$", names(row), value = TRUE)
      )
      row[unique(boundary_columns)] <- NA_real_
      row$source_path <- metadata_row$source_path
      row$checksum_md5 <- metadata_row$checksum_md5
      row$mask_voxels <- metadata_row$mask_voxels
      row$dimensions <- metadata_row$dimensions
      row$affine_match_reference <- metadata_row$affine_match_reference
      row
    })
    context <- rbind(context_base, do.call(rbind, mask_rows))
  }
  rownames(context) <- NULL

  metadata <- list(
    endpoint_analysis = "exact_and_dilated1",
    primary_endpoint_definition = "exact ROI membership",
    sensitivity_endpoint_definition =
      "Chebyshev radius 1 in voxel coordinates (26-neighbor dilation); not a uniform mm tolerance",
    focal_conditioned_discretization = paste0(
      "first focal-positive method in order: ", paste(methods, collapse = ",")
    ),
    unconditional_context_discretization = methods[[1]],
    streamline_voxelization_mode = streamline_voxelization_mode,
    penetration_definition = "unique exact target voxels on the selected full streamline path",
    target_intersection_mode = "full_streamline",
    selected_track_space = track_space,
    selected_axis_flip = as.logical(axis_flip),
    reference_dimensions = as.integer(vol_dim),
    reference_voxel_size_mm = as.numeric(voxel_size),
    context_labels = context_info$labels
  )

  list(summary = summary, penetration = penetration, context = context, metadata = metadata)
}

ftq_validate_endpoint_tables <- function(endpoint_summary,
                                         endpoint_penetration,
                                         endpoint_context,
                                         roi_summary = NULL,
                                         tolerance = 1e-8) {
  checks <- list()
  add <- function(check, passed, detail = "") {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = check,
      passed = isTRUE(passed),
      detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  }
  warnings <- character()
  summary <- ftq_read_table(endpoint_summary)
  penetration <- ftq_read_table(endpoint_penetration)
  context <- ftq_read_table(endpoint_context)
  roi <- if (is.null(roi_summary)) NULL else ftq_read_table(roi_summary)

  summary_required <- c(
    "target", "tolerance", "target_intersection_streamlines",
    "focal_endpoint_streamlines", "paired_endpoint_streamlines",
    "target_intersection_with_paired_endpoints", "shared_mask_ambiguous",
    "paired_focal_target_endpoints", "target_endpoint_only",
    "focal_endpoint_only", "neither_relevant_endpoint",
    "paired_percent_of_target_intersections",
    "paired_percent_of_focal_endpoint_streamlines",
    "focal_target_overlap_voxels"
  )
  add(
    "endpoint_summary_required_columns",
    all(summary_required %in% names(summary)),
    paste(setdiff(summary_required, names(summary)), collapse = ",")
  )
  if (!all(summary_required %in% names(summary))) {
    return(list(
      checks = do.call(rbind, checks),
      biological_qc_status = "unavailable",
      biological_qc_warnings = "Endpoint summary is incomplete."
    ))
  }

  key <- paste(summary$target, summary$tolerance, sep = "\r")
  add("endpoint_summary_unique_target_tolerance", !anyDuplicated(key))
  observed_tolerances <- split(as.character(summary$tolerance), summary$target)
  add(
    "endpoint_summary_exact_and_dilated1_present",
    length(observed_tolerances) > 0L && all(vapply(
      observed_tolerances,
      function(x) setequal(x, c("exact", "dilated1")),
      logical(1)
    ))
  )

  count_columns <- c(
    "target_intersection_streamlines", "focal_endpoint_streamlines",
    "paired_endpoint_streamlines", "target_intersection_with_paired_endpoints",
    "shared_mask_ambiguous", "paired_focal_target_endpoints",
    "target_endpoint_only", "focal_endpoint_only", "neither_relevant_endpoint",
    "focal_target_overlap_voxels"
  )
  counts <- as.matrix(summary[count_columns])
  storage.mode(counts) <- "numeric"
  add("endpoint_summary_counts_finite_nonnegative", all(is.finite(counts) & counts >= 0))

  class_sum <- rowSums(summary[c(
    "shared_mask_ambiguous", "paired_focal_target_endpoints",
    "target_endpoint_only", "focal_endpoint_only", "neither_relevant_endpoint"
  )])
  add(
    "endpoint_categories_partition_target_intersections",
    all(abs(class_sum - summary$target_intersection_streamlines) <= tolerance),
    paste("max_abs_difference =", max(abs(class_sum - summary$target_intersection_streamlines)))
  )
  paired_class_difference <-
    summary$target_intersection_with_paired_endpoints -
      summary$paired_focal_target_endpoints
  add(
    "unambiguous_paired_class_subset_of_all_paired_intersections",
    all(paired_class_difference >= -tolerance &
          paired_class_difference <= summary$shared_mask_ambiguous + tolerance)
  )
  add(
    "paired_intersection_subset_of_target_intersections",
    all(summary$target_intersection_with_paired_endpoints <=
          summary$target_intersection_streamlines + tolerance)
  )

  percent_equal <- function(observed, numerator, denominator) {
    expected <- ftq_endpoint_percent(numerator, denominator)
    same_na <- is.na(observed) & is.na(expected)
    same_value <- is.finite(observed) & is.finite(expected) &
      abs(observed - expected) <= tolerance
    all(same_na | same_value)
  }
  add(
    "paired_percent_of_target_intersections_recomputes",
    percent_equal(
      summary$paired_percent_of_target_intersections,
      summary$target_intersection_with_paired_endpoints,
      summary$target_intersection_streamlines
    )
  )
  add(
    "paired_percent_of_focal_endpoint_streamlines_recomputes",
    percent_equal(
      summary$paired_percent_of_focal_endpoint_streamlines,
      summary$paired_endpoint_streamlines,
      summary$focal_endpoint_streamlines
    )
  )
  percent_columns <- grep("percent", names(summary), value = TRUE)
  percent_matrix <- as.matrix(summary[percent_columns])
  storage.mode(percent_matrix) <- "numeric"
  add(
    "endpoint_summary_percentages_finite_or_missing",
    all(is.na(percent_matrix) | is.finite(percent_matrix))
  )

  exact <- summary[summary$tolerance == "exact", , drop = FALSE]
  d1 <- summary[summary$tolerance == "dilated1", , drop = FALSE]
  d1 <- d1[match(exact$target, d1$target), , drop = FALSE]
  matched <- nrow(exact) == nrow(d1) && !anyNA(d1$target) &&
    identical(as.character(exact$target), as.character(d1$target))
  add("endpoint_tolerances_target_matched", matched)
  if (matched) {
    add(
      "exact_focal_endpoint_count_not_above_dilated1",
      all(exact$focal_endpoint_streamlines <= d1$focal_endpoint_streamlines + tolerance)
    )
    add(
      "exact_paired_endpoint_count_not_above_dilated1",
      all(exact$paired_endpoint_streamlines <= d1$paired_endpoint_streamlines + tolerance)
    )
    add(
      "exact_paired_endpoints_equal_exact_paired_intersections",
      all(abs(
        exact$paired_endpoint_streamlines -
          exact$target_intersection_with_paired_endpoints
      ) <= tolerance)
    )
    intersection_consistent <- vapply(exact$target, function(target) {
      values <- summary$target_intersection_streamlines[summary$target == target]
      length(values) == 2L && max(values) - min(values) <= tolerance
    }, logical(1))
    add("target_intersection_denominator_tolerance_invariant", all(intersection_consistent))
  }

  if (!is.null(roi) && nrow(roi)) {
    roi_mode <- if ("target_intersection_mode" %in% names(roi)) {
      tryCatch(
        ftq_normalize_target_intersection_mode(roi$target_intersection_mode[[1L]]),
        error = function(e) NA_character_
      )
    } else {
      "full_streamline"
    }
    if (identical(roi_mode, "full_streamline")) {
      roi_ok <- vapply(unique(summary$target), function(target) {
        column <- paste0("unique_target_streamlines:", target)
        if (!column %in% names(roi)) return(FALSE)
        observed <- unique(
          summary$target_intersection_streamlines[summary$target == target]
        )
        length(observed) == 1L &&
          abs(as.numeric(observed) - as.numeric(roi[[column]][1])) <= tolerance
      }, logical(1))
      add("endpoint_intersections_match_roi_unique_counts", all(roi_ok))
    } else {
      roi_ok <- vapply(unique(summary$target), function(target) {
        column <- paste0("unique_full_path_target_streamlines:", target)
        if (!column %in% names(roi)) return(FALSE)
        observed <- unique(
          summary$target_intersection_streamlines[summary$target == target]
        )
        length(observed) == 1L &&
          abs(as.numeric(observed) - as.numeric(roi[[column]][1])) <=
            tolerance
      }, logical(1))
      add(
        "endpoint_denominator_preserved_as_full_path",
        !is.na(roi_mode) && all(roi_ok),
        paste0(
          "endpoint denominators were reconciled to independently persisted ",
          "full-streamline counts while primary ROI target counts use ",
          roi_mode
        )
      )
    }
  }

  penetration_required <- c(
    "target", "unique_target_voxels", "streamline_count",
    "target_intersection_streamlines", "percent_of_target_intersections"
  )
  add(
    "endpoint_penetration_required_columns",
    all(penetration_required %in% names(penetration)),
    paste(setdiff(penetration_required, names(penetration)), collapse = ",")
  )
  if (all(penetration_required %in% names(penetration))) {
    add(
      "endpoint_penetration_values_valid",
      !nrow(penetration) || all(
        is.finite(penetration$unique_target_voxels) &
          penetration$unique_target_voxels >= 1 &
          is.finite(penetration$streamline_count) &
          penetration$streamline_count >= 0
      )
    )
    penetration_targets <- unique(summary$target)
    penetration_sum_ok <- vapply(penetration_targets, function(target) {
      expected <- unique(summary$target_intersection_streamlines[
        summary$target == target & summary$tolerance == "exact"
      ])
      observed <- sum(penetration$streamline_count[penetration$target == target])
      length(expected) == 1L && abs(observed - expected) <= tolerance
    }, logical(1))
    add("endpoint_penetration_sums_to_target_intersections", all(penetration_sum_ok))
    penetration_percent_ok <- !nrow(penetration) || percent_equal(
      penetration$percent_of_target_intersections,
      penetration$streamline_count,
      penetration$target_intersection_streamlines
    )
    add("endpoint_penetration_percentages_recompute", penetration_percent_ok)
  }

  context_required <- c(
    "context", "context_type", "valid_streamlines", "endpoint_count",
    "endpoints_inside", "endpoints_outside",
    "streamlines_with_any_endpoint_inside",
    "streamlines_with_both_endpoints_inside",
    "endpoints_nonfinite", "endpoints_one_voxel_outside",
    "endpoints_more_than_one_voxel_outside"
  )
  add(
    "endpoint_context_required_columns",
    all(context_required %in% names(context)),
    paste(setdiff(context_required, names(context)), collapse = ",")
  )
  status <- "unavailable"
  if (all(context_required %in% names(context))) {
    grid <- context[context$context_type == "reference_grid", , drop = FALSE]
    add("endpoint_context_has_one_reference_grid_row", nrow(grid) == 1L)
    context_arithmetic <- all(abs(
      context$endpoints_inside + context$endpoints_outside - context$endpoint_count
    ) <= tolerance)
    add("endpoint_context_inside_outside_reconciles", context_arithmetic)
    bounds_ok <- all(
      context$streamlines_with_both_endpoints_inside <=
        context$streamlines_with_any_endpoint_inside + tolerance &
        context$streamlines_with_any_endpoint_inside <=
        context$valid_streamlines + tolerance
    )
    add("endpoint_context_streamline_counts_bounded", bounds_ok)

    if (nrow(grid) == 1L) {
      if (grid$endpoints_one_voxel_outside > 0) {
        warnings <- c(
          warnings,
          paste(grid$endpoints_one_voxel_outside, "endpoint(s) are exactly one voxel outside the reference grid.")
        )
      }
      if (grid$endpoints_more_than_one_voxel_outside > 0) {
        warnings <- c(
          warnings,
          paste(grid$endpoints_more_than_one_voxel_outside, "endpoint(s) are more than one voxel outside the reference grid.")
        )
      }
      if (grid$endpoints_nonfinite > 0) {
        warnings <- c(warnings, paste(grid$endpoints_nonfinite, "endpoint(s) have non-finite mapped coordinates."))
      }
      has_anatomical_context <- any(context$context_type == "binary_mask")
      has_shared_ambiguity <- any(
        summary$shared_mask_ambiguous[summary$tolerance == "exact"] > 0
      )
      if (!has_anatomical_context) {
        warnings <- c(warnings, "No explicit binary anatomical endpoint context mask was supplied.")
        status <- "unavailable"
      } else if (grid$endpoints_more_than_one_voxel_outside > 0 ||
                 grid$endpoints_nonfinite > 0 || has_shared_ambiguity) {
        status <- "review"
      } else {
        status <- "pass"
      }
      if (has_shared_ambiguity) {
        warnings <- c(warnings, "Focal-target endpoint membership is ambiguous for at least one intersection because masks overlap.")
      }
    }
  }

  checks_df <- if (length(checks)) do.call(rbind, checks) else data.frame(
    check = character(), passed = logical(), detail = character(), stringsAsFactors = FALSE
  )
  rownames(checks_df) <- NULL
  list(
    checks = checks_df,
    biological_qc_status = status,
    biological_qc_warnings = unique(warnings)
  )
}

ftq_endpoint_weighted_quantile <- function(values, weights, probability) {
  if (!length(values) || !sum(weights) || !is.finite(probability)) return(NA_real_)
  order_index <- order(values)
  values <- as.numeric(values[order_index])
  weights <- as.numeric(weights[order_index])
  values[[which(cumsum(weights) >= probability * sum(weights))[1L]]]
}

ftq_endpoint_penetration_summary <- function(penetration, target_labels = NULL) {
  if (is.null(target_labels)) target_labels <- unique(as.character(penetration$target))
  rows <- lapply(target_labels, function(target) {
    part <- penetration[penetration$target == target, , drop = FALSE]
    denominator <- if (nrow(part)) unique(as.numeric(part$target_intersection_streamlines)) else 0
    if (length(denominator) != 1L) denominator <- sum(as.numeric(part$streamline_count))
    values <- as.numeric(part$unique_target_voxels)
    weights <- as.numeric(part$streamline_count)
    count_where <- function(predicate) sum(weights[predicate], na.rm = TRUE)
    one <- count_where(values == 1)
    two <- count_where(values == 2)
    three_to_four <- count_where(values >= 3 & values <= 4)
    five_to_nine <- count_where(values >= 5 & values <= 9)
    ten_or_more <- count_where(values >= 10)
    at_most_two <- count_where(values <= 2)
    at_least_three <- count_where(values >= 3)
    at_least_five <- count_where(values >= 5)
    at_least_ten <- count_where(values >= 10)
    data.frame(
      target = target,
      target_intersection_streamlines = denominator,
      one_voxel_touch_streamlines = one,
      one_voxel_touch_percent_of_target_intersections =
        ftq_endpoint_percent(one, denominator),
      two_voxel_only_streamlines = two,
      two_voxel_only_percent_of_target_intersections =
        ftq_endpoint_percent(two, denominator),
      three_to_four_voxels_streamlines = three_to_four,
      three_to_four_voxels_percent_of_target_intersections =
        ftq_endpoint_percent(three_to_four, denominator),
      five_to_nine_voxels_streamlines = five_to_nine,
      five_to_nine_voxels_percent_of_target_intersections =
        ftq_endpoint_percent(five_to_nine, denominator),
      ten_or_more_voxels_streamlines = ten_or_more,
      ten_or_more_voxels_percent_of_target_intersections =
        ftq_endpoint_percent(ten_or_more, denominator),
      at_most_two_voxels_streamlines = at_most_two,
      at_most_two_voxels_percent_of_target_intersections =
        ftq_endpoint_percent(at_most_two, denominator),
      at_least_three_voxels_streamlines = at_least_three,
      at_least_three_voxels_percent_of_target_intersections =
        ftq_endpoint_percent(at_least_three, denominator),
      at_least_five_voxels_streamlines = at_least_five,
      at_least_five_voxels_percent_of_target_intersections =
        ftq_endpoint_percent(at_least_five, denominator),
      at_least_ten_voxels_streamlines = at_least_ten,
      at_least_ten_voxels_percent_of_target_intersections =
        ftq_endpoint_percent(at_least_ten, denominator),
      penetration_p25_unique_target_voxels =
        ftq_endpoint_weighted_quantile(values, weights, 0.25),
      penetration_median_unique_target_voxels =
        ftq_endpoint_weighted_quantile(values, weights, 0.50),
      penetration_p75_unique_target_voxels =
        ftq_endpoint_weighted_quantile(values, weights, 0.75),
      penetration_p95_unique_target_voxels =
        ftq_endpoint_weighted_quantile(values, weights, 0.95),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

ftq_plot_endpoint_classification <- function(summary, path) {
  exact <- summary[summary$tolerance == "exact", , drop = FALSE]
  if (!nrow(exact)) return(invisible(NULL))
  classes <- c(
    "paired_focal_target_endpoints",
    "shared_mask_ambiguous",
    "target_endpoint_only",
    "focal_endpoint_only",
    "neither_relevant_endpoint"
  )
  values <- t(as.matrix(exact[classes]))
  denominators <- as.numeric(exact$target_intersection_streamlines)
  values <- 100 * values / matrix(
    rep(denominators, each = nrow(values)),
    nrow = nrow(values)
  )
  undefined <- !is.finite(denominators) | denominators <= 0 |
    colSums(!is.finite(values)) > 0L
  values[, undefined] <- NA_real_
  colors <- c("#16827A", "#C24156", "#D8872D", "#4B78A8", "#7A7F87")
  grDevices::png(ftq_png_filename(path), width = 2200, height = 1500, res = 220)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(matrix(c(1, 2), nrow = 2L), heights = c(1, 8))
  graphics::par(oma = c(2, 0, 3, 0))
  ftq_endpoint_plot_legend(
    legend = c(
      "Unambiguous paired",
      "Shared-mask ambiguity",
      "cNode end voxel only",
      "fNode end voxel only",
      "Neither"
    ),
    colors = colors,
    cex = 0.82
  )
  graphics::par(mar = c(9, 5, 2, 2))
  bp <- graphics::barplot(
    values,
    names.arg = exact$target,
    las = 2,
    col = colors,
    border = NA,
    ylim = c(0, 100),
    ylab = "Percent of full-path target intersections"
  )
  if (any(undefined)) graphics::text(bp[undefined], 4, "NA", col = "#9f1239", cex = 0.9)
  ftq_endpoint_plot_availability_note()
  graphics::mtext(
    "Exact end-voxel classification of full-path target intersections",
    outer = TRUE,
    line = 0.5,
    font = 2
  )
  invisible(path)
}

ftq_plot_endpoint_comparison <- function(summary, path) {
  exact <- summary[summary$tolerance == "exact", , drop = FALSE]
  d1 <- summary[summary$tolerance == "dilated1", , drop = FALSE]
  d1 <- d1[match(exact$target, d1$target), , drop = FALSE]
  if (!nrow(exact) || nrow(exact) != nrow(d1)) return(invisible(NULL))
  intersection_values <- rbind(
    Exact = exact$paired_percent_of_target_intersections,
    `1-voxel sensitivity` = d1$paired_percent_of_target_intersections
  )
  focal_values <- rbind(
    Exact = exact$paired_percent_of_focal_endpoint_streamlines,
    `1-voxel sensitivity` = d1$paired_percent_of_focal_endpoint_streamlines
  )
  intersection_denominators <- rbind(exact$target_intersection_streamlines,
                                      d1$target_intersection_streamlines)
  focal_denominators <- rbind(exact$focal_endpoint_streamlines,
                             d1$focal_endpoint_streamlines)
  intersection_values[!is.finite(intersection_denominators) |
                        intersection_denominators <= 0 | !is.finite(intersection_values)] <- NA_real_
  focal_values[!is.finite(focal_denominators) |
                 focal_denominators <= 0 | !is.finite(focal_values)] <- NA_real_
  finite_values <- c(intersection_values, focal_values)
  finite_values <- finite_values[is.finite(finite_values)]
  ymax <- if (length(finite_values)) max(finite_values) else 100
  if (ymax <= 0) ymax <- 100
  upper <- min(100, max(5, ymax * 1.1))
  grDevices::png(ftq_png_filename(path), width = 2400, height = 1500, res = 220)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(
    matrix(c(1, 1, 2, 3), nrow = 2L, byrow = TRUE),
    heights = c(1, 8)
  )
  graphics::par(oma = c(2, 0, 4, 0))
  colors <- c("#16827A", "#D8872D")
  ftq_endpoint_plot_legend(
    legend = rownames(intersection_values),
    colors = colors
  )
  graphics::par(mar = c(9, 5, 3, 1))
  intersection_bp <- graphics::barplot(
    intersection_values,
    beside = TRUE,
    names.arg = exact$target,
    las = 2,
    col = colors,
    border = NA,
    ylim = c(0, upper),
    ylab = "Percent",
    main = "Among full-path intersections"
  )
  missing_intersections <- !is.finite(intersection_values)
  if (any(missing_intersections)) graphics::text(
    intersection_bp[missing_intersections], upper * 0.04, "NA", col = "#9f1239", cex = 0.75
  )
  graphics::par(mar = c(9, 5, 3, 1))
  focal_bp <- graphics::barplot(
    focal_values,
    beside = TRUE,
    names.arg = exact$target,
    las = 2,
    col = colors,
    border = NA,
    ylim = c(0, upper),
    ylab = "Percent",
    main = "Among fNode-end streamlines"
  )
  missing_focal <- !is.finite(focal_values)
  if (any(missing_focal)) graphics::text(
    focal_bp[missing_focal], upper * 0.04, "NA", col = "#9f1239", cex = 0.75
  )
  ftq_endpoint_plot_availability_note()
  graphics::mtext(
    "Geometric fNode-cNode end pairs: exact and 1-voxel dilation",
    outer = TRUE,
    line = 0.5,
    font = 2
  )
  graphics::mtext(
    "Geometric pairs include shared-mask membership; the classification's paired class excludes ambiguity.",
    outer = TRUE, line = 1.7, cex = 0.74
  )
  invisible(path)
}

ftq_plot_endpoint_penetration <- function(penetration_summary, path) {
  if (!nrow(penetration_summary)) return(invisible(NULL))
  values <- rbind(
    `1 voxel` = penetration_summary$one_voxel_touch_percent_of_target_intersections,
    `2 voxels` = penetration_summary$two_voxel_only_percent_of_target_intersections,
    `3-4 voxels` = penetration_summary$three_to_four_voxels_percent_of_target_intersections,
    `5-9 voxels` = penetration_summary$five_to_nine_voxels_percent_of_target_intersections,
    `10+ voxels` = penetration_summary$ten_or_more_voxels_percent_of_target_intersections
  )
  denominator <- as.numeric(penetration_summary$target_intersection_streamlines)
  undefined <- !is.finite(denominator) | denominator <= 0 | colSums(!is.finite(values)) > 0L
  values[, undefined] <- NA_real_
  grDevices::png(ftq_png_filename(path), width = 2200, height = 1500, res = 220)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(matrix(c(1, 2), nrow = 2L), heights = c(1, 8))
  graphics::par(oma = c(2, 0, 4, 0))
  colors <- c("#C24156", "#D8872D", "#16827A", "#4B78A8", "#6F5B95")
  ftq_endpoint_plot_legend(
    legend = rownames(values),
    colors = colors,
    cex = 0.86
  )
  graphics::par(mar = c(9, 5, 2, 2))
  bp <- graphics::barplot(
    values,
    beside = TRUE,
    col = colors,
    border = NA,
    ylim = c(0, 100),
    axes = FALSE,
    axisnames = FALSE
  )
  graphics::axis(2, at = seq(0, 100, by = 20), las = 1)
  graphics::axis(
    1,
    at = colMeans(bp),
    labels = penetration_summary$target,
    las = 2,
    tick = FALSE
  )
  if (any(undefined)) graphics::text(colMeans(bp)[undefined], 4, "NA", col = "#9f1239", cex = 0.9)
  ftq_endpoint_plot_availability_note()
  graphics::box(bty = "l")
  graphics::mtext(
    "Percent of full-path cNode intersections",
    side = 2,
    line = 3
  )
  graphics::mtext(
    "Mutually exclusive cNode-penetration bins for full-path intersections",
    outer = TRUE,
    line = 0.5,
    font = 2
  )
  graphics::mtext("Penetration counts distinct exact target voxels, not distance in millimetres.",
                  outer = TRUE, line = 1.7, cex = 0.76)
  invisible(path)
}

ftq_endpoint_plot_availability_note <- function() {
  graphics::mtext(
    "NA: percentage unavailable (zero/absent denominator or missing data). Unmarked zero-height bars mean 0%.",
    side = 1, outer = TRUE, line = 0.5, cex = 0.72
  )
  invisible(NULL)
}

ftq_endpoint_plot_legend <- function(legend, colors, cex = 0.9) {
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::legend(
    "center",
    horiz = TRUE,
    legend = legend,
    fill = colors,
    border = NA,
    bty = "n",
    cex = cex
  )
  invisible(NULL)
}
