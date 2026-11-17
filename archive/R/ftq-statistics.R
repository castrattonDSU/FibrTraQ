#' Summarize a FibrTraQ result for subject-level statistics
#'
#' Collapse a per-voxel FibrTraQ table into ROI-level and target-level metrics
#' suitable for publication analyses. Use this before group comparisons so each
#' subject contributes one independent row per target.
#'
#' @param result_csv FibrTraQ per-voxel CSV path or data.frame.
#' @param subject Optional subject identifier.
#' @param group Optional group/condition label.
#' @param roi Optional focal ROI label.
#' @param target_labels Optional target columns to include.
#' @param endpoint_summary Optional endpoint-summary CSV/data frame. When
#'   supplied, an additive `endpoints` table is returned.
#' @param endpoint_penetration Optional endpoint-penetration CSV/data frame.
#'   When supplied with `endpoint_summary`, an additive `penetration` table is
#'   returned.
#'
#' @return A list with `roi`, `targets`, and `voxels` data.frames. When endpoint
#'   inputs are supplied, the list also contains additive `endpoints` and
#'   `penetration` members.
#' @export
ftq_stat_subject_summary <- function(result_csv,
                                     subject = NA_character_,
                                     group = NA_character_,
                                     roi = NA_character_,
                                     target_labels = NULL,
                                     endpoint_summary = NULL,
                                     endpoint_penetration = NULL) {
  result <- ftq_read_table(result_csv)
  targets <- ftq_result_target_columns(result, target_labels)
  target_mat <- as.matrix(result[targets])
  storage.mode(target_mat) <- "numeric"
  target_mat[!is.finite(target_mat)] <- 0

  tracks <- as.numeric(result$tracks)
  local <- as.numeric(result$local)
  leave_no_target <- as.numeric(result$leave_no_target)
  leave_with_target <- tracks - local - leave_no_target
  target_hits <- rowSums(target_mat)
  target_diversity <- rowSums(target_mat > 0)
  entropy <- ftq_target_entropy(target_mat, target_hits)
  entropy_fraction <- ifelse(target_diversity > 1, entropy / log2(target_diversity), 0)
  entropy_fraction[!is.finite(entropy_fraction)] <- 0

  focal_voxel_streamline_associations <- sum(tracks)
  total_target_hits <- sum(target_mat)
  roi_summary <- data.frame(
    subject = subject,
    group = group,
    roi = roi,
    focal_voxels = nrow(result),
    connected_focal_voxels = sum(tracks > 0),
    focal_voxel_streamline_associations = focal_voxel_streamline_associations,
    local_associations = sum(local),
    leave_no_target_associations = sum(leave_no_target),
    leave_with_target_associations = sum(leave_with_target),
    focal_voxel_streamline_target_associations = total_target_hits,
    leave_no_target_fraction = if (focal_voxel_streamline_associations > 0) sum(leave_no_target) / focal_voxel_streamline_associations else NA_real_,
    target_hits_per_focal_voxel_streamline_association = if (focal_voxel_streamline_associations > 0) total_target_hits / focal_voxel_streamline_associations else NA_real_,
    median_target_diversity = stats::median(target_diversity[tracks > 0]),
    median_entropy_fraction = stats::median(entropy_fraction[tracks > 0]),
    stringsAsFactors = FALSE
  )

  target_totals <- colSums(target_mat)
  target_summary <- data.frame(
    subject = subject,
    group = group,
    roi = roi,
    target = targets,
    focal_voxel_streamline_target_associations = as.numeric(target_totals),
    connected_focal_voxels = as.integer(colSums(target_mat > 0)),
    target_fraction_of_focal_voxel_streamline_associations = if (focal_voxel_streamline_associations > 0) as.numeric(target_totals) / focal_voxel_streamline_associations else NA_real_,
    target_fraction_of_target_hits = if (total_target_hits > 0) as.numeric(target_totals) / total_target_hits else NA_real_,
    stringsAsFactors = FALSE
  )
  target_summary <- target_summary[order(-target_summary$focal_voxel_streamline_target_associations, target_summary$target), , drop = FALSE]
  rownames(target_summary) <- NULL

  coords <- as.data.frame(ftq_coord_matrix(result$fNode_Coords))
  voxel_summary <- data.frame(
    subject = subject,
    group = group,
    roi = roi,
    fNode_Coords = result$fNode_Coords,
    x = coords$x,
    y = coords$y,
    z = coords$z,
    tracks = tracks,
    leave_no_target_fraction = ifelse(tracks > 0, leave_no_target / tracks, NA_real_),
    target_hit_events = target_hits,
    target_diversity = target_diversity,
    entropy_fraction = entropy_fraction,
    stringsAsFactors = FALSE
  )

  out <- list(roi = roi_summary, targets = target_summary, voxels = voxel_summary)
  if (!is.null(endpoint_summary)) {
    endpoints <- ftq_read_table(endpoint_summary)
    if (!is.null(target_labels)) endpoints <- endpoints[endpoints$target %in% target_labels, , drop = FALSE]
    endpoints <- cbind(
      data.frame(subject = subject, group = group, roi = roi, stringsAsFactors = FALSE),
      endpoints
    )
    out$endpoints <- endpoints
    if (!is.null(endpoint_penetration)) {
      penetration <- ftq_read_table(endpoint_penetration)
      if (!is.null(target_labels)) penetration <- penetration[penetration$target %in% target_labels, , drop = FALSE]
      penetration_summary <- ftq_endpoint_penetration_summary(
        penetration,
        target_labels = unique(as.character(endpoints$target))
      )
      out$penetration <- cbind(
        data.frame(subject = subject, group = group, roi = roi, stringsAsFactors = FALSE),
        penetration_summary
      )
    }
  }
  out
}

#' Build separate unique-streamline and voxel-association target summaries
#'
#' FibrTraQ reports two different target quantities that must not be mixed:
#' ROI-level unique streamline counts and focal-voxel/streamline/target
#' association counts. This helper creates both tables with explicit
#' denominators for publication and validation reports. Both primary tables
#' follow the target-assignment mode recorded by the source run. When supplied,
#' `endpoint_summary` remains a separate additive full-path endpoint audit and
#' does not replace the selected-mode primary counts.
#'
#' @param result_csv FibrTraQ per-voxel CSV path or data.frame.
#' @param roi_summary ROI summary data.frame or CSV path from a freshly computed
#'   FibrTraQ run.
#' @param target_family Target family label recycled across targets, or a named
#'   vector/data.frame with per-target values.
#' @param target_mask_voxels Optional named vector or data.frame containing
#'   target-mask voxel counts.
#' @param target_checksums Optional named vector or data.frame containing
#'   target-mask checksums.
#' @param target_labels Optional target columns to include.
#' @param endpoint_summary Optional endpoint-summary CSV/data frame. When
#'   supplied, an additive `endpoint_streamlines` table is returned.
#'
#' @return A list with selected-mode `unique_streamlines` and
#'   `voxel_associations` tables. When `endpoint_summary` is supplied, the list
#'   also includes the additive full-path `endpoint_streamlines` audit table.
#' @export
ftq_target_summary_tables <- function(result_csv,
                                      roi_summary,
                                      target_family = "unspecified",
                                      target_mask_voxels = NULL,
                                      target_checksums = NULL,
                                      target_labels = NULL,
                                      endpoint_summary = NULL) {
  result <- ftq_read_table(result_csv)
  if (!is.data.frame(roi_summary)) roi_summary <- ftq_read_table(roi_summary)
  if (!nrow(roi_summary)) stop("roi_summary must contain at least one row.")
  targets <- ftq_result_target_columns(result, target_labels)
  target_mat <- as.matrix(result[targets])
  storage.mode(target_mat) <- "numeric"
  target_mat[!is.finite(target_mat)] <- 0

  family <- ftq_target_metadata_value(targets, target_family, default = "unspecified")
  mask_voxels <- ftq_target_metadata_value(targets, target_mask_voxels, default = NA_real_)
  checksum <- ftq_target_metadata_value(targets, target_checksums, default = NA_character_)

  unique_cols <- paste0("unique_target_streamlines:", targets)
  unique_counts <- vapply(unique_cols, function(nm) {
    if (!nm %in% names(roi_summary)) return(NA_real_)
    as.numeric(roi_summary[[nm]][1])
  }, numeric(1))
  unique_denom <- as.numeric(roi_summary$unique_focal_intersecting_streamlines[1])

  association_counts <- colSums(target_mat)
  association_denom <- sum(as.numeric(result$tracks))

  unique_streamlines <- data.frame(
    target = targets,
    target_family = family,
    unique_target_positive_focal_intersecting_streamlines = unique_counts,
    percent_of_unique_focal_intersecting_streamlines = if (unique_denom > 0) 100 * unique_counts / unique_denom else NA_real_,
    target_mask_voxels = as.numeric(mask_voxels),
    checksum = as.character(checksum),
    stringsAsFactors = FALSE
  )
  unique_streamlines <- unique_streamlines[order(
    -unique_streamlines$unique_target_positive_focal_intersecting_streamlines,
    unique_streamlines$target
  ), , drop = FALSE]
  rownames(unique_streamlines) <- NULL

  voxel_associations <- data.frame(
    target = targets,
    target_family = family,
    focal_voxel_streamline_target_associations = as.numeric(association_counts),
    percent_of_all_focal_voxel_streamline_associations = if (association_denom > 0) 100 * association_counts / association_denom else NA_real_,
    focal_voxels_with_at_least_one_association = as.integer(colSums(target_mat > 0)),
    target_mask_voxels = as.numeric(mask_voxels),
    stringsAsFactors = FALSE
  )
  voxel_associations <- voxel_associations[order(
    -voxel_associations$focal_voxel_streamline_target_associations,
    voxel_associations$target
  ), , drop = FALSE]
  rownames(voxel_associations) <- NULL

  out <- list(
    unique_streamlines = unique_streamlines,
    voxel_associations = voxel_associations
  )
  if (!is.null(endpoint_summary)) {
    endpoint_streamlines <- ftq_read_table(endpoint_summary)
    endpoint_streamlines <- endpoint_streamlines[
      endpoint_streamlines$target %in% targets,
      ,
      drop = FALSE
    ]
    endpoint_streamlines$target_family <- family[
      match(endpoint_streamlines$target, targets)
    ]
    endpoint_streamlines$checksum <- checksum[
      match(endpoint_streamlines$target, targets)
    ]
    out$endpoint_streamlines <- endpoint_streamlines
  }
  out
}

ftq_target_metadata_value <- function(targets, value, default = NA) {
  if (is.null(value)) return(rep(default, length(targets)))
  if (is.data.frame(value)) {
    name_col <- intersect(c("target", "label", "name"), names(value))[1]
    val_col <- setdiff(names(value), name_col)[1]
    if (is.na(name_col) || is.na(val_col)) return(rep(default, length(targets)))
    idx <- match(targets, as.character(value[[name_col]]))
    out <- value[[val_col]][idx]
    out[is.na(idx)] <- default
    return(out)
  }
  if (length(value) == 1L && is.null(names(value))) return(rep(value, length(targets)))
  if (!is.null(names(value))) {
    out <- unname(value[targets])
    out[is.na(out)] <- default
    return(out)
  }
  rep(value, length.out = length(targets))
}

#' Estimate focal-voxel spatial gradients
#'
#' Fit target-fraction or entropy gradients against focal voxel x/y/z
#' coordinates. For publication inference, run this per subject and test the
#' resulting subject-level slopes.
#'
#' @param result_csv FibrTraQ per-voxel CSV path or data.frame.
#' @param target_labels Optional target columns to include.
#' @param min_connected_voxels Minimum nonzero voxels required per target.
#'
#' @return A data.frame of target-by-axis slopes.
#' @export
ftq_stat_spatial_gradients <- function(result_csv,
                                       target_labels = NULL,
                                       min_connected_voxels = 5L) {
  result <- ftq_read_table(result_csv)
  targets <- ftq_result_target_columns(result, target_labels)
  coords <- as.data.frame(ftq_coord_matrix(result$fNode_Coords))
  tracks <- as.numeric(result$tracks)
  connected <- tracks > 0
  rows <- list()

  for (target in targets) {
    hits <- as.numeric(result[[target]])
    y <- ifelse(tracks > 0, hits / tracks, NA_real_)
    keep <- connected & is.finite(y)
    if (sum(keep & hits > 0) < min_connected_voxels || length(unique(y[keep])) < 2L) next
    dat <- data.frame(y = y[keep], x = coords$x[keep], ycoord = coords$y[keep], z = coords$z[keep])
    fit <- stats::lm(y ~ x + ycoord + z, data = dat)
    coef_tab <- summary(fit)$coefficients
    for (axis in c("x", "ycoord", "z")) {
      axis_present <- axis %in% rownames(coef_tab)
      rows[[length(rows) + 1L]] <- data.frame(
        target = target,
        axis = if (axis == "ycoord") "y" else axis,
        slope = if (axis_present) unname(coef_tab[axis, "Estimate"]) else NA_real_,
        p_value = if (axis_present) unname(coef_tab[axis, "Pr(>|t|)"]) else NA_real_,
        connected_voxels = sum(keep),
        target_connected_voxels = sum(keep & hits > 0),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame(
    target = character(), axis = character(), slope = numeric(), p_value = numeric(),
    connected_voxels = integer(), target_connected_voxels = integer()
  )
  out$q_value <- stats::p.adjust(out$p_value, method = "BH")
  out
}

#' Compare correct runs against negative controls
#'
#' @param correct Data frame of subject/run summaries.
#' @param control Data frame of matched negative-control summaries.
#' @param by Column(s) used to pair rows.
#' @param metric Metric column to compare, or one of `"target_hit_rate"` and
#'   `"leave_no_target_fraction"` if raw count columns are present.
#' @param alternative Alternative hypothesis passed to [stats::wilcox.test()].
#'
#' @return A list with paired rows, effect size, and paired Wilcoxon result.
#' @export
ftq_stat_negative_control <- function(correct,
                                      control,
                                      by,
                                      metric = "target_hit_rate",
                                      alternative = c("greater", "two.sided", "less")) {
  alternative <- match.arg(alternative)
  correct_metric <- ftq_stat_metric(correct, metric)
  control_metric <- ftq_stat_metric(control, metric)
  cdat <- cbind(correct[by], correct_metric = correct_metric)
  kdat <- cbind(control[by], control_metric = control_metric)
  paired <- merge(cdat, kdat, by = by, all = FALSE)
  if (!nrow(paired)) stop("No matched rows found for negative-control comparison.")
  paired$difference <- paired$correct_metric - paired$control_metric
  test <- if (nrow(paired) >= 2L && any(paired$difference != 0)) {
    stats::wilcox.test(paired$correct_metric, paired$control_metric, paired = TRUE, alternative = alternative, exact = FALSE)
  } else {
    NULL
  }
  list(
    metric = metric,
    paired = paired,
    median_difference = stats::median(paired$difference, na.rm = TRUE),
    mean_difference = mean(paired$difference, na.rm = TRUE),
    test = test
  )
}

#' Compare subject-level metrics between groups
#'
#' @param data Subject-level summary data.
#' @param metric Metric column to compare.
#' @param group_col Grouping column.
#' @param subject_col Subject column.
#' @param method `"wilcox"` or `"t"`.
#'
#' @return A data.frame with one row per metric comparison.
#' @export
ftq_stat_compare_groups <- function(data,
                                    metric,
                                    group_col = "group",
                                    subject_col = "subject",
                                    method = c("wilcox", "t")) {
  method <- match.arg(method)
  if (!all(c(metric, group_col, subject_col) %in% names(data))) {
    stop("data must contain metric, group_col, and subject_col.")
  }
  groups <- unique(stats::na.omit(data[[group_col]]))
  if (length(groups) != 2L) stop("Exactly two groups are required.")
  x <- data[data[[group_col]] == groups[1], metric]
  y <- data[data[[group_col]] == groups[2], metric]
  test <- if (identical(method, "wilcox")) stats::wilcox.test(x, y, exact = FALSE) else stats::t.test(x, y)
  data.frame(
    metric = metric,
    group_a = groups[1],
    group_b = groups[2],
    n_a = length(stats::na.omit(x)),
    n_b = length(stats::na.omit(y)),
    mean_a = mean(x, na.rm = TRUE),
    mean_b = mean(y, na.rm = TRUE),
    median_a = stats::median(x, na.rm = TRUE),
    median_b = stats::median(y, na.rm = TRUE),
    estimate_difference = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
    p_value = test$p.value,
    method = method,
    stringsAsFactors = FALSE
  )
}

#' Write shifted ROI masks for spatial-null controls
#'
#' @param mask Path to a NIfTI ROI mask.
#' @param output_dir Directory for shifted masks.
#' @param shifts Matrix/data.frame with columns `x`, `y`, `z`, or a list of
#'   numeric length-3 shifts in voxels.
#' @param prefix Output prefix.
#' @param overwrite Whether existing shifted masks may be replaced.
#'
#' @return A data.frame of written mask paths and shifts.
#' @export
ftq_null_shifted_rois <- function(mask,
                                  output_dir,
                                  shifts,
                                  prefix = tools::file_path_sans_ext(basename(mask)),
                                  overwrite = FALSE) {
  img <- RNifti::readNifti(mask, internal = FALSE)
  arr <- as.array(img)
  shift_df <- ftq_shift_table(shifts)
  output_dir <- normalizePath(path.expand(output_dir), mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- lapply(seq_len(nrow(shift_df)), function(i) {
    s <- as.integer(shift_df[i, c("x", "y", "z")])
    shifted <- ftq_shift_array(arr, s)
    label <- paste0("shift-x", s[1], "_y", s[2], "_z", s[3])
    out <- file.path(output_dir, paste0(prefix, "_", label, ".nii.gz"))
    if (file.exists(out) && !isTRUE(overwrite)) {
      stop("Shifted ROI exists. Use overwrite = TRUE to replace: ", out)
    }
    RNifti::writeNifti(shifted, out, template = img)
    data.frame(x = s[1], y = s[2], z = s[3], path = out, nonzero_voxels = sum(shifted != 0), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

ftq_stat_metric <- function(data, metric) {
  if (metric %in% names(data)) return(as.numeric(data[[metric]]))
  if (identical(metric, "target_hit_rate")) {
    numerator <- if ("focal_voxel_streamline_target_associations" %in% names(data)) {
      data$focal_voxel_streamline_target_associations
    } else if ("target_hit_events" %in% names(data)) {
      data$target_hit_events
    } else {
      NULL
    }
    denominator <- if ("focal_voxel_streamline_associations" %in% names(data)) {
      data$focal_voxel_streamline_associations
    } else if ("assigned_tracks" %in% names(data)) {
      data$assigned_tracks
    } else {
      NULL
    }
    if (is.null(numerator) || is.null(denominator)) {
      stop("target_hit_rate requires target and focal-voxel association count columns.")
    }
    return(ifelse(denominator > 0, numerator / denominator, NA_real_))
  }
  if (identical(metric, "leave_no_target_fraction")) {
    if ("leave_no_target_fraction" %in% names(data)) return(as.numeric(data$leave_no_target_fraction))
    numerator <- if ("leave_no_target_associations" %in% names(data)) {
      data$leave_no_target_associations
    } else if ("leave_no_target_tracks" %in% names(data)) {
      data$leave_no_target_tracks
    } else {
      NULL
    }
    denominator <- if ("focal_voxel_streamline_associations" %in% names(data)) {
      data$focal_voxel_streamline_associations
    } else if ("assigned_tracks" %in% names(data)) {
      data$assigned_tracks
    } else {
      NULL
    }
    if (is.null(numerator) || is.null(denominator)) {
      stop("leave_no_target_fraction requires leave-no-target and focal-voxel association count columns.")
    }
    return(ifelse(denominator > 0, numerator / denominator, NA_real_))
  }
  stop("Unknown metric column: ", metric)
}

ftq_shift_table <- function(shifts) {
  if (is.data.frame(shifts) || is.matrix(shifts)) {
    out <- as.data.frame(shifts)
    if (!all(c("x", "y", "z") %in% names(out))) names(out)[1:3] <- c("x", "y", "z")
    out <- out[, c("x", "y", "z"), drop = FALSE]
  } else if (is.list(shifts)) {
    out <- do.call(rbind, lapply(shifts, function(x) as.integer(x[1:3])))
    out <- as.data.frame(out)
    names(out) <- c("x", "y", "z")
  } else {
    stop("shifts must be a data.frame, matrix, or list of length-3 shifts.")
  }
  out[] <- lapply(out, as.integer)
  if (anyNA(out)) stop("shifts must contain integer x/y/z values.")
  out
}

ftq_shift_array <- function(arr, shift) {
  out <- array(0, dim = dim(arr))
  dims <- dim(arr)[1:3]
  src <- lapply(seq_along(dims), function(i) {
    if (shift[i] >= 0L) seq_len(dims[i] - shift[i]) else seq.int(1L - shift[i], dims[i])
  })
  dst <- lapply(seq_along(dims), function(i) src[[i]] + shift[i])
  if (any(vapply(src, length, integer(1)) <= 0L)) return(out)
  out[dst[[1]], dst[[2]], dst[[3]]] <- arr[src[[1]], src[[2]], src[[3]]]
  out
}
