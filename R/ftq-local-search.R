#' Retrieve local-search end-voxel audit tables
#'
#' Local-search runs retain one audit row for every focal-intersecting
#' streamline. The audit records the nearest-node status and scanner-space
#' distance for each end voxel, the paired target when one was assigned, and
#' whether that assigned streamline also intersected the target along its full
#' voxelized path.
#'
#' @param result A data frame returned directly by [FibrTraQ()].
#' @return A list with `assignments`, `summary`, `endpoint_status`, and
#'   `metadata`. For results not computed with `assignment_mode =
#'   "local_search"`, the tables are empty and `metadata$available` is false.
#' @export
ftq_local_search_tables <- function(result) {
  if (!is.data.frame(result)) {
    stop("result must be a data frame returned directly by FibrTraQ().")
  }
  assignments <- attr(result, "fibrtraq_local_search_assignments")
  if (is.null(assignments)) {
    return(list(
      assignments = data.frame(),
      summary = data.frame(),
      endpoint_status = data.frame(),
      metadata = list(available = FALSE)
    ))
  }
  list(
    assignments = assignments,
    summary = ftq_local_search_assignment_summary(assignments),
    endpoint_status = ftq_local_search_endpoint_status(assignments),
    metadata = list(
      available = TRUE,
      search_radius_mm = attr(result, "fibrtraq_local_search_radius_mm"),
      spatial_units = attr(result, "fibrtraq_local_search_spatial_units"),
      spatial_units_scale_to_mm = attr(
        result,
        "fibrtraq_local_search_spatial_units_scale_to_mm"
      ),
      spatial_units_interpretation = attr(
        result,
        "fibrtraq_local_search_spatial_units_interpretation"
      ),
      focal_intersecting_streamlines = nrow(assignments)
    )
  )
}

ftq_local_search_assignment_summary <- function(assignments) {
  required <- c(
    "paired_target", "paired_focal_distance_mm",
    "paired_target_distance_mm", "paired_max_distance_mm",
    "paired_target_intersects_full_path", "search_radius_mm"
  )
  missing <- setdiff(required, names(assignments))
  if (length(missing)) {
    stop(
      "Local-search assignment table is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  paired <- assignments[
    !is.na(assignments$paired_target) & nzchar(assignments$paired_target),
    ,
    drop = FALSE
  ]
  if (!nrow(paired)) {
    return(data.frame(
      target = character(),
      unique_connections = integer(),
      connections_with_full_path_intersection = integer(),
      connections_without_full_path_intersection = integer(),
      full_path_intersection_percent_of_connections = numeric(),
      focal_distance_median_mm = numeric(),
      target_distance_median_mm = numeric(),
      paired_max_distance_q25_mm = numeric(),
      paired_max_distance_median_mm = numeric(),
      paired_max_distance_q75_mm = numeric(),
      paired_max_distance_q95_mm = numeric(),
      paired_max_distance_max_mm = numeric(),
      full_path_target_distance_median_mm = numeric(),
      full_path_max_distance_median_mm = numeric(),
      full_path_max_distance_q95_mm = numeric(),
      proximity_only_target_distance_median_mm = numeric(),
      proximity_only_max_distance_median_mm = numeric(),
      proximity_only_max_distance_q95_mm = numeric(),
      search_radius_mm = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  quantile_value <- function(x, probability) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    as.numeric(stats::quantile(
      x,
      probs = probability,
      names = FALSE,
      type = 7,
      na.rm = TRUE
    ))
  }
  targets <- unique(as.character(paired$paired_target))
  rows <- lapply(targets, function(target) {
    target_rows <- paired[paired$paired_target == target, , drop = FALSE]
    with_full_path <- sum(
      target_rows$paired_target_intersects_full_path %in% TRUE
    )
    without_full_path <- sum(
      target_rows$paired_target_intersects_full_path %in% FALSE
    )
    full_path_rows <- target_rows[
      target_rows$paired_target_intersects_full_path %in% TRUE,
      ,
      drop = FALSE
    ]
    proximity_only_rows <- target_rows[
      target_rows$paired_target_intersects_full_path %in% FALSE,
      ,
      drop = FALSE
    ]
    total <- nrow(target_rows)
    data.frame(
      target = target,
      unique_connections = total,
      connections_with_full_path_intersection = with_full_path,
      connections_without_full_path_intersection = without_full_path,
      full_path_intersection_percent_of_connections = if (total > 0L) {
        100 * with_full_path / total
      } else {
        NA_real_
      },
      focal_distance_median_mm = quantile_value(
        target_rows$paired_focal_distance_mm,
        0.5
      ),
      target_distance_median_mm = quantile_value(
        target_rows$paired_target_distance_mm,
        0.5
      ),
      paired_max_distance_q25_mm = quantile_value(
        target_rows$paired_max_distance_mm,
        0.25
      ),
      paired_max_distance_median_mm = quantile_value(
        target_rows$paired_max_distance_mm,
        0.5
      ),
      paired_max_distance_q75_mm = quantile_value(
        target_rows$paired_max_distance_mm,
        0.75
      ),
      paired_max_distance_q95_mm = quantile_value(
        target_rows$paired_max_distance_mm,
        0.95
      ),
      paired_max_distance_max_mm = quantile_value(
        target_rows$paired_max_distance_mm,
        1
      ),
      full_path_target_distance_median_mm = quantile_value(
        full_path_rows$paired_target_distance_mm,
        0.5
      ),
      full_path_max_distance_median_mm = quantile_value(
        full_path_rows$paired_max_distance_mm,
        0.5
      ),
      full_path_max_distance_q95_mm = quantile_value(
        full_path_rows$paired_max_distance_mm,
        0.95
      ),
      proximity_only_target_distance_median_mm = quantile_value(
        proximity_only_rows$paired_target_distance_mm,
        0.5
      ),
      proximity_only_max_distance_median_mm = quantile_value(
        proximity_only_rows$paired_max_distance_mm,
        0.5
      ),
      proximity_only_max_distance_q95_mm = quantile_value(
        proximity_only_rows$paired_max_distance_mm,
        0.95
      ),
      search_radius_mm = unique(as.numeric(target_rows$search_radius_mm))[[1L]],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

ftq_local_search_endpoint_status <- function(assignments) {
  required <- c("start_status", "end_status")
  missing <- setdiff(required, names(assignments))
  if (length(missing)) {
    stop(
      "Local-search assignment table is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  denominator <- nrow(assignments)
  rows <- lapply(c("start", "end"), function(endpoint) {
    values <- as.character(assignments[[paste0(endpoint, "_status")]])
    counts <- table(values, useNA = "ifany")
    data.frame(
      endpoint = endpoint,
      status = names(counts),
      streamline_count = as.integer(counts),
      percent_of_focal_intersecting_streamlines = if (denominator > 0L) {
        100 * as.integer(counts) / denominator
      } else {
        numeric(length(counts))
      },
      denominator_focal_intersecting_streamlines = denominator,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

ftq_validate_local_search_assignments <- function(assignments,
                                                   roi_summary = NULL,
                                                   target_labels = NULL,
                                                   tolerance = 1e-8) {
  assignments <- ftq_read_table(assignments)
  checks <- list()
  add <- function(check, passed, detail = "") {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = check,
      passed = isTRUE(passed),
      detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  }
  required <- c(
    "streamline_index", "start_status", "start_distance_mm",
    "start_ambiguous", "end_status", "end_distance_mm", "end_ambiguous",
    "paired_target", "paired_focal_distance_mm",
    "paired_target_distance_mm", "paired_max_distance_mm",
    "paired_target_intersects_full_path", "search_radius_mm"
  )
  missing <- setdiff(required, names(assignments))
  add(
    "local_search_required_columns",
    !length(missing),
    paste(missing, collapse = ", ")
  )
  if (length(missing)) return(do.call(rbind, checks))

  indices <- as.numeric(assignments$streamline_index)
  add(
    "local_search_streamline_indices_valid_unique",
    all(is.finite(indices) & indices >= 1 & indices == floor(indices)) &&
      !anyDuplicated(indices),
    paste("rows =", nrow(assignments))
  )
  radii <- unique(as.numeric(assignments$search_radius_mm))
  add(
    "local_search_radius_single_finite_positive_value",
    length(radii) == 1L && is.finite(radii) && radii > 0,
    paste(radii, collapse = ", ")
  )
  radius <- if (length(radii) == 1L) radii[[1L]] else NA_real_

  allowed_status <- c(
    "endpoint_unavailable", "ambiguous_nearest_node",
    "no_node_within_radius", "focal", target_labels
  )
  observed_status <- unique(c(
    as.character(assignments$start_status),
    as.character(assignments$end_status)
  ))
  add(
    "local_search_endpoint_status_values_supported",
    all(observed_status %in% allowed_status),
    paste(setdiff(observed_status, allowed_status), collapse = ", ")
  )

  check_endpoint_distance <- function(endpoint) {
    status <- as.character(assignments[[paste0(endpoint, "_status")]])
    distance <- as.numeric(assignments[[paste0(endpoint, "_distance_mm")]])
    assigned <- status %in% c("ambiguous_nearest_node", "focal", target_labels)
    unassigned <- status %in% c(
      "endpoint_unavailable",
      "no_node_within_radius"
    )
    all(is.finite(distance[assigned]) & distance[assigned] >= 0 &
        distance[assigned] <= radius + tolerance) &&
      all(is.na(distance[unassigned]))
  }
  add(
    "local_search_start_distances_match_status_and_radius",
    is.finite(radius) && check_endpoint_distance("start")
  )
  add(
    "local_search_end_distances_match_status_and_radius",
    is.finite(radius) && check_endpoint_distance("end")
  )

  paired <- !is.na(assignments$paired_target) &
    nzchar(as.character(assignments$paired_target))
  paired_targets <- as.character(assignments$paired_target[paired])
  add(
    "local_search_paired_targets_supported",
    all(paired_targets %in% target_labels),
    paste(setdiff(paired_targets, target_labels), collapse = ", ")
  )
  paired_focal <- as.numeric(assignments$paired_focal_distance_mm[paired])
  paired_target <- as.numeric(assignments$paired_target_distance_mm[paired])
  paired_max <- as.numeric(assignments$paired_max_distance_mm[paired])
  add(
    "local_search_paired_distances_finite_and_within_radius",
    all(is.finite(paired_focal) & paired_focal >= 0 &
          paired_focal <= radius + tolerance) &&
      all(is.finite(paired_target) & paired_target >= 0 &
            paired_target <= radius + tolerance) &&
      all(is.finite(paired_max) & paired_max >= 0 &
            paired_max <= radius + tolerance)
  )
  add(
    "local_search_paired_max_distance_reconciles",
    all(abs(paired_max - pmax(paired_focal, paired_target)) <= tolerance)
  )
  add(
    "local_search_unpaired_rows_have_no_paired_distances",
    all(is.na(assignments$paired_focal_distance_mm[!paired])) &&
      all(is.na(assignments$paired_target_distance_mm[!paired])) &&
      all(is.na(assignments$paired_max_distance_mm[!paired]))
  )
  add(
    "local_search_paired_full_path_status_complete",
    all(!is.na(assignments$paired_target_intersects_full_path[paired])) &&
      all(is.na(assignments$paired_target_intersects_full_path[!paired]))
  )

  if (!is.null(roi_summary)) {
    roi <- ftq_read_table(roi_summary)
    roi <- roi[1L, , drop = FALSE]
    focal_column <- "unique_focal_intersecting_streamlines"
    add(
      "local_search_rows_match_unique_focal_streamlines",
      focal_column %in% names(roi) &&
        nrow(assignments) == as.numeric(roi[[focal_column]][[1L]]),
      if (focal_column %in% names(roi)) {
        paste(
          "assignment rows =", nrow(assignments),
          "; unique focal streamlines =", roi[[focal_column]][[1L]]
        )
      } else {
        paste("missing ROI-summary column", focal_column)
      }
    )
    for (target in target_labels) {
      target_rows <- paired & assignments$paired_target == target
      observed <- sum(target_rows)
      expected_column <- paste0("unique_target_streamlines:", target)
      expected <- if (expected_column %in% names(roi)) {
        as.numeric(roi[[expected_column]][[1L]])
      } else {
        NA_real_
      }
      add(
        paste0("local_search_target_count_reconciles:", target),
        is.finite(expected) && observed == expected,
        paste("assignment rows =", observed, "; ROI summary =", expected)
      )

      with_full_path <- sum(
        target_rows & assignments$paired_target_intersects_full_path %in% TRUE
      )
      without_full_path <- sum(
        target_rows & assignments$paired_target_intersects_full_path %in% FALSE
      )
      add(
        paste0("local_search_full_path_partition_reconciles:", target),
        with_full_path + without_full_path == observed,
        paste(
          "with full-path intersection =", with_full_path,
          "; proximity only =", without_full_path,
          "; primary connections =", observed
        )
      )
      full_path_column <- paste0(
        "unique_full_path_target_streamlines:",
        target
      )
      available_full_path <- if (full_path_column %in% names(roi)) {
        as.numeric(roi[[full_path_column]][[1L]])
      } else {
        NA_real_
      }
      add(
        paste0("local_search_intersection_subset_of_full_path:", target),
        is.finite(available_full_path) &&
          with_full_path <= available_full_path,
        paste(
          "paired connections also intersecting target =", with_full_path,
          "; all full-path target intersections =", available_full_path
        )
      )
    }
  }
  out <- do.call(rbind, checks)
  rownames(out) <- NULL
  out
}

ftq_plot_local_search_distances <- function(assignments, path) {
  summary <- ftq_local_search_assignment_summary(assignments)
  paired <- assignments[
    !is.na(assignments$paired_target) & nzchar(assignments$paired_target),
    ,
    drop = FALSE
  ]
  grDevices::png(ftq_png_filename(path), width = 2200, height = 1000, res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(matrix(seq_len(2L), nrow = 1L), widths = c(1.3, 1))

  if (!nrow(paired)) {
    graphics::par(mar = c(4, 4, 3, 1) + 0.1)
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No paired local-search connections")
    graphics::plot.new()
    return(invisible(path))
  }

  targets <- unique(as.character(paired$paired_target))
  colors <- grDevices::hcl.colors(length(targets), palette = "Dark 3")
  radius <- max(as.numeric(paired$search_radius_mm), na.rm = TRUE)
  graphics::par(mar = c(4.5, 5, 3, 2.5) + 0.1)
  graphics::plot(
    NA,
    xlim = c(0, radius),
    ylim = c(0, 100),
    xlab = "Larger of the two end-voxel distances (mm)",
    ylab = "Connections at or below distance (%)",
    main = "Paired-connection distance distribution"
  )
  graphics::grid(col = "#E5E7EB")
  for (index in seq_along(targets)) {
    target <- targets[[index]]
    values <- sort(as.numeric(
      paired$paired_max_distance_mm[paired$paired_target == target]
    ))
    values <- values[is.finite(values)]
    if (!length(values)) next
    graphics::lines(
      c(0, values),
      c(0, 100 * seq_along(values) / length(values)),
      type = "s",
      lwd = 2.5,
      col = colors[[index]]
    )
  }
  graphics::abline(v = radius, lty = 3, col = "#475569")
  graphics::legend(
    "bottomright",
    legend = targets,
    col = colors,
    lwd = 2.5,
    bty = "n",
    cex = 0.75,
    inset = c(0.02, 0.02)
  )

  partition <- rbind(
    "Also intersects cNode" = summary$connections_with_full_path_intersection,
    "Proximity-only assignment" = summary$connections_without_full_path_intersection
  )
  colnames(partition) <- summary$target
  graphics::par(mar = c(4.5, 5, 3, 1) + 0.1)
  graphics::barplot(
    partition,
    col = c("#2A9D8F", "#E9C46A"),
    border = NA,
    las = 1,
    ylab = "Unique connections",
    main = "Full-path reconciliation"
  )
  graphics::legend(
    "topright",
    legend = rownames(partition),
    fill = c("#2A9D8F", "#E9C46A"),
    bty = "n",
    cex = 0.9
  )
  invisible(path)
}
