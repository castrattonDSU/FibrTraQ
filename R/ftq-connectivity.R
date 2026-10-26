ftq_target_intersection_modes <- function() {
  c("full_streamline", "endpoint_voxels", "local_search")
}

ftq_normalize_target_intersection_mode <- function(mode) {
  if (missing(mode) || is.null(mode)) mode <- "full_streamline"
  mode <- as.character(mode)
  choices <- ftq_target_intersection_modes()
  if (length(mode) > 1L) {
    if (identical(mode, choices)) {
      mode <- mode[[1L]]
    } else {
      stop(
        "target_intersection_mode must be one non-empty value.",
        call. = FALSE
      )
    }
  }
  if (!length(mode) || is.na(mode[[1L]]) || !nzchar(trimws(mode[[1L]]))) {
    stop("target_intersection_mode must be one non-empty value.", call. = FALSE)
  }
  mode <- gsub("-", "_", tolower(trimws(mode[[1L]])), fixed = TRUE)
  aliases <- c(
    all_voxel = "full_streamline",
    all_voxels = "full_streamline",
    end_voxel = "endpoint_voxels",
    end_voxels = "endpoint_voxels",
    endpoint_voxel = "endpoint_voxels"
  )
  if (mode %in% names(aliases)) mode <- unname(aliases[[mode]])
  if (!mode %in% choices) {
    stop(
      "target_intersection_mode must be one of: ",
      paste(choices, collapse = ", "),
      ". The paper-aligned aliases 'all_voxels' and 'end_voxels' are also accepted.",
      call. = FALSE
    )
  }
  mode
}

ftq_validate_local_search_radius <- function(radius_mm) {
  radius_mm <- as.numeric(radius_mm)
  if (length(radius_mm) != 1L || is.na(radius_mm) ||
      !is.finite(radius_mm) || radius_mm <= 0) {
    stop(
      "local_search_radius_mm must be one finite number greater than zero.",
      call. = FALSE
    )
  }
  radius_mm
}

ftq_local_search_unit_info <- function(source) {
  units <- tryCatch(
    as.character(RNifti::pixunits(source))[[1L]],
    error = function(e) "Unknown"
  )
  if (!length(units) || is.na(units) || !nzchar(trimws(units))) {
    units <- "Unknown"
  }
  normalized <- tolower(trimws(units))
  scale_to_mm <- switch(
    normalized,
    "m" = 1000,
    "mm" = 1,
    "um" = 0.001,
    1
  )
  known <- normalized %in% c("m", "mm", "um")
  interpretation <- switch(
    normalized,
    "m" = "converted_m_to_mm",
    "mm" = "header_mm",
    "um" = "converted_um_to_mm",
    "unknown_assumed_mm"
  )
  list(
    spatial_units = units,
    spatial_units_known = known,
    scale_to_mm = scale_to_mm,
    units_interpretation = interpretation
  )
}

ftq_local_search_xform_mm <- function(source) {
  xform <- unclass(RNifti::xform(source))
  if (!is.matrix(xform) || !identical(dim(xform), c(4L, 4L)) ||
      any(!is.finite(xform))) {
    stop("The reference-image affine is unavailable or invalid.", call. = FALSE)
  }
  unit_info <- ftq_local_search_unit_info(source)
  xform[1:3, ] <- xform[1:3, , drop = FALSE] * unit_info$scale_to_mm
  xform
}

ftq_local_search_affine_equal <- function(a, b, tolerance = 1e-6) {
  tryCatch(
    isTRUE(all.equal(
      ftq_local_search_xform_mm(a),
      ftq_local_search_xform_mm(b),
      tolerance = tolerance,
      check.attributes = FALSE
    )),
    error = function(e) FALSE
  )
}

ftq_world_mm_to_voxel <- function(coords_mm, reference) {
  coords_mm <- as.matrix(coords_mm)[, 1:3, drop = FALSE]
  scale_to_mm <- ftq_local_search_unit_info(reference)$scale_to_mm
  world_native <- coords_mm / scale_to_mm
  matrix(
    RNifti::worldToVoxel(world_native, reference, simple = FALSE),
    ncol = 3L
  )
}

ftq_voxel_to_world_mm <- function(voxel, reference) {
  voxel <- as.matrix(voxel)[, 1:3, drop = FALSE]
  scale_to_mm <- ftq_local_search_unit_info(reference)$scale_to_mm
  matrix(
    RNifti::voxelToWorld(voxel, reference, simple = FALSE),
    ncol = 3L
  ) * scale_to_mm
}

ftq_local_search_geometry <- function(fNode, ref_image = NULL) {
  reference_image <- NULL
  if (!is.null(ref_image)) {
    reference_image <- tryCatch(
      {
        if (inherits(ref_image, "niftiImage")) {
          ref_image
        } else {
          RNifti::readNifti(ref_image, internal = FALSE)
        }
      },
      error = function(e) NULL
    )
    focal_dims <- dim(fNode)
    reference_dims <- tryCatch(
      {
        dims <- dim(reference_image)
        as.integer(dims[1:3])
      },
      error = function(e) NULL
    )
    if (!is.null(focal_dims) && length(focal_dims) >= 3L &&
        !is.null(reference_dims) &&
        !identical(as.integer(focal_dims[1:3]), reference_dims)) {
      stop(
        "For local search, ref_image dimensions must match the fNode grid.",
        call. = FALSE
      )
    }
    if (inherits(fNode, "niftiImage") &&
        inherits(reference_image, "niftiImage") &&
        !ftq_local_search_affine_equal(fNode, reference_image)) {
      stop(
        "For local search, ref_image must use the same affine/grid as fNode.",
        call. = FALSE
      )
    }
  }
  stored_basis <- attr(fNode, "fibrtraq_voxel_to_world_basis", exact = TRUE)
  unit_info <- list(
    spatial_units = "mm",
    spatial_units_known = TRUE,
    scale_to_mm = 1,
    units_interpretation = "stored_geometry_declared_mm"
  )
  if (!is.null(stored_basis)) {
    basis <- as.matrix(stored_basis)
  } else if (!is.null(ref_image) || inherits(fNode, "niftiImage")) {
    source <- if (is.null(ref_image)) fNode else reference_image
    basis_points <- rbind(
      c(0, 0, 0),
      c(1, 0, 0),
      c(0, 1, 0),
      c(0, 0, 1)
    )
    world <- tryCatch(
      matrix(
        RNifti::voxelToWorld(basis_points, source, simple = FALSE),
        ncol = 3L
      ),
      error = function(e) NULL
    )
    basis <- if (is.null(world)) {
      NULL
    } else {
      sweep(world[2:4, , drop = FALSE], 2L, world[1L, ], `-`)
    }
    if (!is.null(basis)) {
      unit_info <- ftq_local_search_unit_info(source)
      basis <- basis * unit_info$scale_to_mm
    }
  } else {
    voxel_size <- attr(fNode, "fibrtraq_voxel_size_mm", exact = TRUE)
    basis <- if (is.null(voxel_size)) {
      NULL
    } else {
      diag(as.numeric(voxel_size)[1:3], nrow = 3L)
    }
  }

  if (is.null(basis) || !identical(dim(basis), c(3L, 3L)) ||
      any(!is.finite(basis))) {
    stop(
      "target_intersection_mode = 'local_search' requires reference-image ",
      "spatial metadata. Supply ref_image (or an RNifti fNode) so distances ",
      "can be measured in scanner-space millimetres.",
      call. = FALSE
    )
  }
  storage.mode(basis) <- "double"
  gram <- basis %*% t(basis)
  inverse_gram <- tryCatch(solve(gram), error = function(e) NULL)
  if (is.null(inverse_gram) || any(!is.finite(inverse_gram)) ||
      any(diag(inverse_gram) <= 0)) {
    stop(
      "The reference-image voxel-to-world geometry is singular or invalid ",
      "for local-search distances.",
      call. = FALSE
    )
  }
  list(
    voxel_to_world = basis,
    voxel_radius_per_mm = sqrt(diag(inverse_gram)),
    voxel_size_mm = sqrt(rowSums(basis^2)),
    spatial_units = unit_info$spatial_units,
    spatial_units_known = unit_info$spatial_units_known,
    spatial_units_scale_to_mm = unit_info$scale_to_mm,
    spatial_units_interpretation = unit_info$units_interpretation
  )
}

ftq_local_search_voxel_radius <- function(radius_mm, geometry) {
  extent <- radius_mm * as.numeric(geometry$voxel_radius_per_mm)
  if (length(extent) != 3L || anyNA(extent) ||
      any(!is.finite(extent)) || any(extent < 0)) {
    stop(
      "local_search_radius_mm produces a non-finite voxel search extent ",
      "under the reference-image geometry.",
      call. = FALSE
    )
  }
  extent
}

ftq_node_labels_at_key <- function(key, lut_f, target_lut) {
  if (length(key) != 1L || is.na(key) || key < 1L || key > length(lut_f)) {
    return(integer())
  }
  labels <- if (isTRUE(lut_f[[key]])) 1L else integer()
  if (ncol(target_lut)) {
    labels <- c(labels, which(target_lut[key, , drop = TRUE]) + 1L)
  }
  unique(as.integer(labels))
}

ftq_unique_node_label_at_key <- function(key, lut_f, target_lut) {
  labels <- ftq_node_labels_at_key(key, lut_f, target_lut)
  if (length(labels) == 1L) labels[[1L]] else NA_integer_
}

ftq_paired_target_hits <- function(start_label, end_label, n_targets) {
  out <- rep(FALSE, n_targets)
  if (is.na(start_label) || is.na(end_label) || !n_targets) return(out)
  if (start_label == 1L && end_label > 1L && end_label <= n_targets + 1L) {
    out[[end_label - 1L]] <- TRUE
  } else if (end_label == 1L &&
             start_label > 1L && start_label <= n_targets + 1L) {
    out[[start_label - 1L]] <- TRUE
  }
  out
}

ftq_endpoint_voxel_target_hits_r <- function(chosen_keys,
                                              two_endpoints,
                                              lut_f,
                                              target_lut) {
  n_targets <- ncol(target_lut)
  if (!isTRUE(two_endpoints) || length(chosen_keys) < 2L) {
    return(rep(FALSE, n_targets))
  }
  start_label <- ftq_unique_node_label_at_key(
    chosen_keys[[1L]], lut_f, target_lut
  )
  end_label <- ftq_unique_node_label_at_key(
    chosen_keys[[length(chosen_keys)]], lut_f, target_lut
  )
  ftq_paired_target_hits(start_label, end_label, n_targets)
}

ftq_apply_continuous_flip <- function(coords, axis_flip, vol_dim) {
  out <- as.matrix(coords)[, 1:3, drop = FALSE]
  for (axis in seq_len(3L)) {
    if (isTRUE(axis_flip[[axis]])) {
      out[, axis] <- vol_dim[[axis]] - out[, axis] + 1
    }
  }
  out
}

ftq_local_search_endpoint_label_r <- function(point,
                                               lut_f,
                                               target_lut,
                                               vol_dim,
                                               radius_mm,
                                               geometry) {
  empty <- list(
    label = NA_integer_,
    ambiguous = FALSE,
    distance_mm = NA_real_,
    voxel_key = NA_integer_,
    voxel = rep(NA_real_, 3L)
  )
  point <- as.numeric(point)[1:3]
  if (length(point) != 3L || any(!is.finite(point))) return(empty)

  voxel_extent <- radius_mm * geometry$voxel_radius_per_mm
  lower <- pmax(1L, ceiling(point - voxel_extent))
  upper <- pmin(as.integer(vol_dim), floor(point + voxel_extent))
  if (any(lower > upper)) return(empty)

  radius_squared <- radius_mm^2
  best_squared <- Inf
  best_labels <- integer()
  best_key <- NA_integer_
  best_voxel <- rep(NA_real_, 3L)

  for (z in seq.int(lower[[3L]], upper[[3L]])) {
    for (y in seq.int(lower[[2L]], upper[[2L]])) {
      for (x in seq.int(lower[[1L]], upper[[1L]])) {
        key <- as.integer(
          x + (y - 1L) * vol_dim[[1L]] +
            (z - 1L) * vol_dim[[1L]] * vol_dim[[2L]]
        )
        labels <- ftq_node_labels_at_key(key, lut_f, target_lut)
        if (!length(labels)) next
        delta_world <- as.numeric(
          (c(x, y, z) - point) %*% geometry$voxel_to_world
        )
        distance_squared <- sum(delta_world^2)
        radius_tolerance <- 1e-12 * max(
          1,
          abs(distance_squared),
          abs(radius_squared)
        )
        if (distance_squared > radius_squared + radius_tolerance) next
        if (!is.finite(best_squared)) {
          best_squared <- distance_squared
          best_labels <- labels
          best_key <- key
          best_voxel <- c(x, y, z)
          next
        }
        nearest_tolerance <- 1e-12 * max(
          1,
          abs(distance_squared),
          abs(best_squared)
        )
        if (distance_squared < best_squared - nearest_tolerance) {
          best_squared <- distance_squared
          best_labels <- labels
          best_key <- key
          best_voxel <- c(x, y, z)
        } else if (
          abs(distance_squared - best_squared) <= nearest_tolerance
        ) {
          best_labels <- union(best_labels, labels)
        }
      }
    }
  }
  if (!is.finite(best_squared)) return(empty)
  list(
    label = if (length(best_labels) == 1L) best_labels[[1L]] else NA_integer_,
    ambiguous = length(best_labels) > 1L,
    distance_mm = sqrt(best_squared),
    voxel_key = best_key,
    voxel = best_voxel
  )
}

ftq_local_search_target_hits_r <- function(mapped,
                                            two_endpoints,
                                            lut_f,
                                            target_lut,
                                            axis_flip,
                                            vol_dim,
                                            radius_mm,
                                            geometry) {
  n_targets <- ncol(target_lut)
  if (!isTRUE(two_endpoints) || nrow(mapped) < 2L) {
    return(list(
      target_hit = rep(FALSE, n_targets),
      endpoints = list()
    ))
  }
  endpoints <- ftq_apply_continuous_flip(
    mapped[c(1L, nrow(mapped)), 1:3, drop = FALSE],
    axis_flip,
    vol_dim
  )
  assignments <- lapply(seq_len(2L), function(i) {
    ftq_local_search_endpoint_label_r(
      endpoints[i, ],
      lut_f = lut_f,
      target_lut = target_lut,
      vol_dim = vol_dim,
      radius_mm = radius_mm,
      geometry = geometry
    )
  })
  list(
    target_hit = ftq_paired_target_hits(
      assignments[[1L]]$label,
      assignments[[2L]]$label,
      n_targets
    ),
    endpoints = assignments
  )
}

ftq_primary_target_hits_r <- function(mode,
                                      full_path_target_hit,
                                      chosen_keys,
                                      mapped,
                                      two_endpoints,
                                      lut_f,
                                      target_lut,
                                      axis_flip,
                                      vol_dim,
                                      local_search_radius_mm,
                                      local_search_geometry = NULL) {
  mode <- ftq_normalize_target_intersection_mode(mode)
  if (identical(mode, "full_streamline")) {
    return(list(
      target_hit = as.logical(full_path_target_hit),
      endpoints = list()
    ))
  }
  if (identical(mode, "endpoint_voxels")) {
    return(list(
      target_hit = ftq_endpoint_voxel_target_hits_r(
        chosen_keys,
        two_endpoints,
        lut_f,
        target_lut
      ),
      endpoints = list()
    ))
  }
  if (is.null(local_search_geometry)) {
    stop("Internal error: local-search geometry was not prepared.", call. = FALSE)
  }
  ftq_local_search_target_hits_r(
    mapped = mapped,
    two_endpoints = two_endpoints,
    lut_f = lut_f,
    target_lut = target_lut,
    axis_flip = axis_flip,
    vol_dim = vol_dim,
    radius_mm = local_search_radius_mm,
    geometry = local_search_geometry
  )
}
