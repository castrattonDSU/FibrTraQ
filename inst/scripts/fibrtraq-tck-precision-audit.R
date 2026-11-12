#!/usr/bin/env Rscript
# Exact accounting of representation differences, never a relaxed parity test.
# Rscript fibrtraq-tck-precision-audit.R <prepared_dir> <results_dir> <new_audit_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || any(args %in% c("--help", "-h"))) {
  cat("Usage: Rscript fibrtraq-tck-precision-audit.R <prepared_dir> <results_dir> <new_audit_dir>\n")
  quit(status = if (any(args %in% c("--help", "-h"))) 0L else 1L)
}
suppressPackageStartupMessages(library(FibrTraQ))
prepared <- normalizePath(args[[1L]], mustWork = TRUE)
results <- normalizePath(args[[2L]], mustWork = TRUE)
audit <- normalizePath(args[[3L]], mustWork = FALSE)
if (file.exists(audit)) stop("Choose a NEW audit directory; existing evidence is preserved.")
if (!dir.create(audit, recursive = TRUE, mode = "0700")) stop("Cannot create audit directory.")
audit <- normalizePath(audit, mustWork = TRUE)
get_ftq <- function(name) getFromNamespace(name, "FibrTraQ")
checks <- data.frame(check = character(), passed = logical(), detail = character())
state <- list(status = "running", exact_format_parity = FALSE,
              count_tolerance = 0,
              world_coordinate_tolerance_mm = 1e-4,
              scope = paste("Causal accounting of Float32 representation differences for the pinned public example.",
                            "Independent voxel-path reconstruction and affected-streamline core reruns are used;",
                            "this is not an independent counting oracle or an anatomical accuracy claim."),
              prepared_directory = prepared, results_directory = results,
              package_version = as.character(packageVersion("FibrTraQ")),
              package_description = as.list(utils::packageDescription("FibrTraQ")[c("Version", "Built", "RemoteSha")]),
              invocation = commandArgs(),
              checks = checks)
write_state <- function() {
  state$checks <- checks
  path <- file.path(audit, "precision-audit.json")
  jsonlite::write_json(state, path, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA)
  Sys.chmod(path, "0600")
  utils::write.csv(checks, file.path(audit, "audit-checks.csv"), row.names = FALSE)
  Sys.chmod(file.path(audit, "audit-checks.csv"), "0600")
}
check <- function(name, passed, detail = "") {
  checks <<- rbind(checks, data.frame(check = name, passed = isTRUE(passed), detail = detail))
  write_state()
  if (!isTRUE(passed)) stop("Unaccounted discrepancy: ", name, if (nzchar(detail)) paste0("; ", detail))
}
csv <- function(path) utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
one_file <- function(root, suffix) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE)
  files <- files[endsWith(files, suffix)]
  if (length(files) != 1L) stop("Expected one ", suffix, " in ", root)
  files[[1L]]
}
resolve_bundle <- function(kind) {
  root <- file.path(results, paste0(kind, "-stream-cpp"), "sub-CON02", "negative_x_interface")
  metadata_path <- one_file(root, "_fibrtraq_run-metadata.json")
  metadata <- jsonlite::fromJSON(metadata_path, simplifyVector = FALSE)
  local_file <- function(path) {
    matches <- list.files(root, recursive = TRUE, full.names = TRUE)
    matches <- matches[basename(matches) == basename(path)]
    if (length(matches) != 1L) stop("Cannot resolve result artifact: ", path)
    matches[[1L]]
  }
  metadata$files$map_paths <- lapply(metadata$files$map_paths, local_file)
  for (name in setdiff(names(metadata$files), c("map_paths", "context"))) {
    if (!is.null(metadata$files[[name]])) metadata$files[[name]] <- local_file(metadata$files[[name]])
  }
  list(root = root, metadata_path = metadata_path, metadata = metadata)
}
serialize_table <- function(table, name) {
  path <- file.path(audit, paste0(name, ".csv"))
  utils::write.csv(table, path, row.names = FALSE, na = "NA")
  Sys.chmod(path, "0600")
  csv(path)
}
compare_table <- function(expected, observed, name) {
  # Match the package's CSV serialization of derived percentages. All counts
  # remain integers and are compared exactly; no numerical tolerance is used.
  expected <- serialize_table(expected, paste0("reconstructed-", name))
  rownames(expected) <- rownames(observed) <- NULL
  difference <- all.equal(expected, observed, tolerance = 0, check.attributes = TRUE)
  check(name, isTRUE(difference), if (isTRUE(difference)) "Exact serialized table match" else paste(difference, collapse = "; "))
}
linear_key <- function(points, dims) {
  points <- round(points)
  valid <- rowSums(points >= 1 & sweep(points, 2L, dims, "<=")) == 3L
  kept <- points[valid, , drop = FALSE]
  list(keys = unique(kept[, 1L] + dims[[1L]] * (kept[, 2L] - 1) +
                       prod(dims[1:2]) * (kept[, 3L] - 1)), outside = sum(!valid))
}
map_points <- function(points, affine) {
  # This operation order reproduces the documented compiled affine mapping.
  matrix(vapply(seq_len(3L), function(axis) points[, 1L] * affine[axis, 1L] +
           points[, 2L] * affine[axis, 2L] + points[, 3L] * affine[axis, 3L] +
           affine[axis, 4L], numeric(nrow(points))), ncol = 3L)
}
step_count <- function(points) {
  if (nrow(points) < 2L) return(integer())
  pmax(1L, ceiling(apply(abs(points[-1L, , drop = FALSE] -
                              points[-nrow(points), , drop = FALSE]), 1L, max) / 0.5))
}
dense_path <- function(points, steps) {
  if (!length(steps)) return(points)
  pieces <- lapply(seq_along(steps), function(segment) {
    fraction <- seq.int(0L, steps[[segment]] - 1L) / steps[[segment]]
    sweep(fraction %o% (points[segment + 1L, ] - points[segment, ]), 2L, points[segment, ], "+")
  })
  rbind(do.call(rbind, pieces), points[nrow(points), , drop = FALSE])
}

main <- function() {
  a <- resolve_bundle("tck"); b <- resolve_bundle("trk")
  expected_contract <- list(discretization = "round", focal_assignment_mode = "all_focal_voxels",
                            streamline_voxelization_mode = "segment_interpolation",
                            target_intersection_mode = "full_streamline", endpoint_analysis = "exact_and_dilated1")
  for (bundle in list(a, b)) for (name in names(expected_contract)) {
    check(paste0(basename(dirname(dirname(bundle$root))), ":contract:", name),
          identical(bundle$metadata$manifest[[name]], expected_contract[[name]]))
  }
  paths <- c(tck = file.path(prepared, "tracks_10k_6fed19354195.tck"),
             trk = file.path(prepared, "native-equivalent.trk"),
             reference = file.path(prepared, "reference.nii.gz"),
             focal = file.path(prepared, "focal-left-interface.nii.gz"),
             positive_x_cortical_GM = file.path(prepared, "target-right-gm.nii.gz"),
             negative_x_cortical_GM = file.path(prepared, "target-left-gm.nii.gz"))
  input_hashes <- stats::setNames(unname(tools::md5sum(paths)), names(paths))
  state$input_md5 <<- as.list(stats::setNames(unname(input_hashes), names(paths)))
  if (!exists("sha256sum", envir = asNamespace("tools"), inherits = FALSE)) {
    stop("This optional provenance audit requires an R version providing tools::sha256sum.")
  }
  sha256 <- getFromNamespace("sha256sum", "tools")
  script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  if (length(script_path) != 1L || !file.exists(script_path)) stop("Cannot identify the audit script for provenance.")
  state$audit_script_sha256 <<- unname(sha256(script_path))
  provenance_path <- file.path(prepared, "provenance.json")
  preparation <- jsonlite::fromJSON(provenance_path, simplifyVector = FALSE)
  input_sha256 <- stats::setNames(unname(sha256(paths)), names(paths))
  state$input_sha256 <<- as.list(stats::setNames(unname(input_sha256), names(paths)))
  state$preparation_provenance_sha256 <<- unname(sha256(provenance_path))
  check("independent_preparation_completed", identical(preparation$status, "prepared_and_independently_checked"))
  check("independent_conversion_within_fixed_bound", is.finite(preparation$independent_tractogram_check$trk_coordinate_max_abs_difference_mm) &&
          preparation$independent_tractogram_check$trk_coordinate_max_abs_difference_mm <= 1e-4 &&
          isTRUE(preparation$independent_tractogram_check$streamline_and_point_order_preserved))
  for (name in names(paths)) {
    recorded <- preparation$derived_files[[basename(paths[[name]])]]
    check(paste0("prepared_sha256:", name), !is.null(recorded) &&
            identical(unname(input_sha256[[name]]), recorded$sha256) &&
            identical(as.numeric(file.info(paths[[name]])$size), as.numeric(recorded$bytes)))
  }
  source_tck <- Filter(function(x) identical(x$role, "native_tractogram"), preparation$source_specification$sources)
  check("pinned_native_tck_source", length(source_tck) == 1L &&
          identical(source_tck[[1L]]$sha256, unname(input_sha256[["tck"]])))
  for (j in 1:2) {
    bundle <- list(a, b)[[j]]; kind <- c("tck", "trk")[[j]]
    manifest <- bundle$metadata$manifest
    check(paste0(kind, ":coordinate_contract"),
          identical(manifest$track_space, c("world", "trackvis")[[j]]) &&
            identical(as.logical(unlist(manifest$axis_flip)), rep(FALSE, 3L)))
    roles <- c(fnode = "focal", tractogram = kind, ref_image = "reference", background_image = "reference",
               "target:positive_x_cortical_GM" = "positive_x_cortical_GM",
               "target:negative_x_cortical_GM" = "negative_x_cortical_GM")
    for (role in names(roles)) {
      recorded <- Filter(function(x) identical(x$name, role), bundle$metadata$inputs)
      check(paste0(kind, ":run_input_binding:", role), length(recorded) == 1L &&
              identical(recorded[[1L]]$checksum_md5, unname(input_hashes[[roles[[role]]]])))
    }
  }
  result_files <- unique(c(a$metadata_path, b$metadata_path,
    unlist(a$metadata$files), unlist(b$metadata$files),
    list.files(file.path(results, "tck-density-context"), recursive = TRUE, full.names = TRUE),
    list.files(file.path(results, "trk-density-context"), recursive = TRUE, full.names = TRUE)))
  result_hashes <- tools::md5sum(result_files)
  state$result_evidence_md5 <<- as.list(stats::setNames(unname(result_hashes), substring(result_files, nchar(results) + 2L)))
  reference <- RNifti::readNifti(paths[["reference"]], internal = FALSE)
  focal <- RNifti::readNifti(paths[["focal"]], internal = FALSE)
  targets <- lapply(paths[c("positive_x_cortical_GM", "negative_x_cortical_GM")], RNifti::readNifti, internal = FALSE)
  dims <- dim(reference)
  native <- get_ftq("ftq_read_tractogram")(paths[["tck"]], "memory")
  equivalent <- get_ftq("ftq_read_tractogram")(paths[["trk"]], "memory")
  n <- length(native$tracks)
  check("paired_streamline_counts", n == length(equivalent$tracks) && n > 0)
  track_affine <- get_ftq("ftq_trackvis_affine")(equivalent$header)
  affines <- list(
    get_ftq("ftq_qc_raw_to_voxel_affine")("world", reference),
    get_ftq("ftq_qc_raw_to_voxel_affine")("trackvis", reference, track_affine = track_affine)
  )
  densities <- list(numeric(prod(dims)), numeric(prod(dims)))
  statistics <- list(c(stored = 0, densified = 0, outside = 0), c(stored = 0, densified = 0, outside = 0))
  path_rows <- boundary_rows <- step_rows <- list()
  affected <- integer(); path_affected <- integer(); endpoint_only <- integer()
  max_world_delta <- max_voxel_delta <- 0
  for (i in seq_len(n)) {
    raw_a <- native$tracks[[i]]$coords; raw_b <- equivalent$tracks[[i]]$coords
    if (!identical(dim(raw_a), dim(raw_b))) stop("Point counts changed at streamline ", i)
    p <- map_points(raw_a, affines[[1L]]); q <- map_points(raw_b, affines[[2L]])
    max_world_delta <- max(max_world_delta, abs(raw_a - map_points(raw_b, track_affine)))
    point_delta <- apply(abs(p - q), 2L, max)
    max_voxel_delta <- max(max_voxel_delta, point_delta)
    steps_a <- step_count(p); steps_b <- step_count(q)
    pd <- dense_path(p, steps_a); qd <- dense_path(q, steps_b)
    paths_pair <- list(linear_key(pd, dims), linear_key(qd, dims))
    for (j in 1:2) {
      densities[[j]][paths_pair[[j]]$keys] <- densities[[j]][paths_pair[[j]]$keys] + 1
      statistics[[j]] <- statistics[[j]] + c(stored = nrow(p), densified = if (j == 1L) nrow(pd) else nrow(qd),
                                             outside = paths_pair[[j]]$outside)
    }
    different_path <- !setequal(paths_pair[[1L]]$keys, paths_pair[[2L]]$keys)
    different_endpoint <- any(round(p[c(1L, nrow(p)), , drop = FALSE]) != round(q[c(1L, nrow(q)), , drop = FALSE]))
    if (different_path || different_endpoint) affected <- c(affected, i)
    if (different_endpoint && !different_path) endpoint_only <- c(endpoint_only, i)
    if (different_path) {
      path_affected <- c(path_affected, i)
      for (j in 1:2) for (key in setdiff(paths_pair[[j]]$keys, paths_pair[[3L - j]]$keys)) {
        voxel <- arrayInd(key, dims)
        path_rows[[length(path_rows) + 1L]] <- data.frame(streamline = i, x = voxel[1L], y = voxel[2L], z = voxel[3L],
                                                        delta_trk_minus_tck = if (j == 1L) -1L else 1L)
      }
    }
    changed_steps <- which(steps_a != steps_b)
    for (segment in changed_steps) {
      span_a <- max(abs(p[segment + 1L, ] - p[segment, ]))
      span_b <- max(abs(q[segment + 1L, ] - q[segment, ]))
      threshold <- min(steps_a[[segment]], steps_b[[segment]]) * 0.5
      span_bound <- max(abs(p[segment + 1L, ] - q[segment + 1L, ]) + abs(p[segment, ] - q[segment, ]))
      accounted <- abs(steps_a[[segment]] - steps_b[[segment]]) == 1 &&
        threshold >= min(span_a, span_b) && threshold <= max(span_a, span_b) &&
        abs(span_a - span_b) <= span_bound + 64 * .Machine$double.eps * max(1, abs(c(p, q)))
      step_rows[[length(step_rows) + 1L]] <- data.frame(streamline = i, segment = segment,
        tck_steps = steps_a[[segment]], trk_steps = steps_b[[segment]], tck_max_axis_span = span_a,
        trk_max_axis_span = span_b, threshold_voxels = threshold, endpoint_perturbation_bound = span_bound,
        occupancy_changed = different_path, accounted = accounted)
    }
    if (different_path && length(changed_steps)) stop("Changed interpolation grid with changed occupancy requires additional accounting at streamline ", i)
    if (different_path && !length(changed_steps)) {
      changed <- which(round(pd) != round(qd), arr.ind = TRUE)
      for (k in seq_len(nrow(changed))) {
        row <- changed[k, 1L]; axis <- changed[k, 2L]
        va <- pd[row, axis]; vb <- qd[row, axis]
        threshold <- (round(va) + round(vb)) / 2
        bound <- point_delta[[axis]]
        accounted <- abs(round(va) - round(vb)) == 1 && threshold >= min(va, vb) &&
          threshold <= max(va, vb) && abs(va - vb) <= bound + 64 * .Machine$double.eps * max(1, abs(c(va, vb)))
        boundary_rows[[length(boundary_rows) + 1L]] <- data.frame(streamline = i, densified_point = row,
          axis = c("x", "y", "z")[[axis]], tck_coordinate = va, trk_coordinate = vb,
          rounding_boundary = threshold, tck_signed_distance = va - threshold,
          trk_signed_distance = vb - threshold, coordinate_difference = abs(va - vb),
          measured_endpoint_perturbation_bound = bound, accounted = accounted)
      }
    }
    if (i %% 1000L == 0L) message("Audited paired streamline geometry: ", i, "/", n)
  }
  state$geometry <<- list(streamlines = n, path_difference_indices = path_affected,
                          affected_core_indices = affected, max_world_difference_mm = max_world_delta,
                          endpoint_only_difference_indices = endpoint_only,
                          max_voxel_difference = max_voxel_delta,
                          interpolation_step_change_indices = unique(vapply(step_rows, function(x) x$streamline, integer(1))))
  boundaries <- if (length(boundary_rows)) do.call(rbind, boundary_rows) else data.frame()
  steps <- if (length(step_rows)) do.call(rbind, step_rows) else data.frame()
  path_deltas <- if (length(path_rows)) do.call(rbind, path_rows) else data.frame()
  serialize_table(boundaries, "rounding-boundary-crossings")
  serialize_table(steps, "interpolation-step-crossings")
  serialize_table(path_deltas, "streamline-occupancy-deltas")
  check("world_coordinate_difference_within_fixed_bound", is.finite(max_world_delta) && max_world_delta <= 1e-4,
        "Fixed 0.0001 mm bound from the independent preparation contract; never applied to counts.")
  check("no_unexamined_endpoint_only_crossings", !length(endpoint_only),
        "This fixed public example has no changed endpoint voxels without changed path occupancy.")
  check("all_rounding_crossings_bounded", nrow(boundaries) > 0 && all(boundaries$accounted))
  check("all_interpolation_step_changes_accounted", !nrow(steps) || all(steps$accounted))
  check("nonzero_representation_difference_retained", length(path_affected) > 0)
  density_summary <- function(values, stats) {
    occupied <- which(values > 0)
    coordinates <- arrayInd(occupied, dims)
    lo <- apply(coordinates, 2L, min); hi <- apply(coordinates, 2L, max)
    world <- get_ftq("ftq_bbox_world")(focal, lo, hi)
    data.frame(processed_streamlines = n, declared_streamlines = n,
      nonzero_tractogram_occupied_voxels = length(occupied),
      tractogram_occupancy_volume_mm3 = length(occupied) * get_ftq("ftq_voxel_volume_mm3")(focal),
      total_streamline_voxel_occupancy_events = sum(values), max_density = max(values),
      total_sample_points = unname(stats[["stored"]]),
      total_voxelization_points = unname(stats[["densified"]]),
      out_of_bounds_sample_points = unname(stats[["outside"]]), mapping_error_streamlines = 0,
      voxel_bbox_x_min = lo[[1L]], voxel_bbox_x_max = hi[[1L]],
      voxel_bbox_y_min = lo[[2L]], voxel_bbox_y_max = hi[[2L]],
      voxel_bbox_z_min = lo[[3L]], voxel_bbox_z_max = hi[[3L]],
      world_bbox_x_min = world[[1L]], world_bbox_x_max = world[[2L]],
      world_bbox_y_min = world[[3L]], world_bbox_y_max = world[[4L]],
      world_bbox_z_min = world[[5L]], world_bbox_z_max = world[[6L]],
      note = "Tractogram-derived density context is not anatomical brain volume.")
  }
  observed_density <- list()
  for (j in 1:2) {
    kind <- c("tck", "trk")[[j]]
    context <- file.path(results, paste0(kind, "-density-context"), "context")
    density_path <- one_file(context, "_tractogram-density.nii.gz")
    image <- RNifti::readNifti(density_path, internal = FALSE)
    observed_density[[j]] <- as.numeric(image)
    check(paste0(kind, ":full_density_reconstruction"), identical(dim(image), dims) &&
            identical(observed_density[[j]], densities[[j]]))
    compare_table(density_summary(densities[[j]], statistics[[j]]),
                  csv(one_file(context, "_tractogram-context-summary.csv")), paste0(kind, "-density-summary"))
  }
  delta <- observed_density[[2L]] - observed_density[[1L]]
  indices <- which(delta != 0)
  coordinates <- arrayInd(indices, dims)
  serialize_table(data.frame(x = coordinates[, 1L], y = coordinates[, 2L], z = coordinates[, 3L],
                             native_tck = observed_density[[1L]][indices],
                             equivalent_trk = observed_density[[2L]][indices],
                             delta_trk_minus_tck = delta[indices]), "observed-density-deltas")
  check("all_density_deltas_accounted", identical(delta, densities[[2L]] - densities[[1L]]))

  # Re-run only paths/endpoints whose discrete representation changed. Every
  # unaffected streamline has identical path occupancy and endpoint voxels.
  run_subset <- function(object, space) {
    object$tracks <- object$tracks[affected]
    FibrTraQ(fNode = focal, cNode_list = targets, cNode_labels = names(targets),
              all_tracks = object, ref_image = reference, track_space = space,
              axis_flip = rep(FALSE, 3L), discretization = "round",
              focal_assignment_mode = "all_focal_voxels",
              streamline_voxelization_mode = "segment_interpolation",
              assignment_mode = "full_streamline", endpoint_analysis = "exact_and_dilated1",
              backend = "cpp", verbose = FALSE, progress = FALSE)
  }
  subsets <- list(run_subset(native, "world"), run_subset(equivalent, "trackvis"))
  core_a <- csv(a$metadata$files$result_csv); core_b <- csv(b$metadata$files$result_csv)
  check("core_coordinate_row_identity", identical(core_a$fNode_Coords, core_b$fNode_Coords) &&
          identical(core_a$fNode_Coords, subsets[[1L]]$fNode_Coords) &&
          identical(core_a$fNode_Coords, subsets[[2L]]$fNode_Coords))
  counts <- setdiff(names(core_a), c("fNode_Coords", "none"))
  predicted_core <- core_a
  for (name in counts) predicted_core[[name]] <- core_a[[name]] + subsets[[2L]][[name]] - subsets[[1L]][[name]]
  predicted_core$none <- as.integer(predicted_core$tracks == 0)
  compare_table(predicted_core, core_b, "per-voxel-hits")
  changed_core <- list()
  for (name in c(counts, "none")) for (row in which(core_a[[name]] != core_b[[name]])) {
    changed_core[[length(changed_core) + 1L]] <- data.frame(fNode_Coords = core_a$fNode_Coords[[row]],
      field = name, native_tck = core_a[[name]][[row]], equivalent_trk = core_b[[name]][[row]],
      predicted_trk = predicted_core[[name]][[row]])
  }
  serialize_table(do.call(rbind, changed_core), "core-count-deltas")
  map_values <- function(core) c(list(tracks = core$tracks, local = core$local,
    "leave-no-target" = core$leave_no_target,
    "leave-with-target" = core$tracks - core$local - core$leave_no_target,
    "tracks-binary" = as.integer(core$tracks > 0)),
    stats::setNames(as.list(core[names(targets)]), paste0("target-", names(targets))))
  core_coordinates <- get_ftq("ftq_coord_matrix")(core_a$fNode_Coords)
  for (j in 1:2) {
    bundle <- list(a, b)[[j]]
    values <- map_values(list(core_a, predicted_core)[[j]])
    check(paste0(c("tck", "trk")[[j]], ":map_names"), setequal(names(values), names(bundle$metadata$files$map_paths)))
    for (name in names(values)) {
      expected <- array(0, dims); expected[core_coordinates] <- values[[name]]
      observed <- RNifti::readNifti(bundle$metadata$files$map_paths[[name]], internal = FALSE)
      check(paste0(c("tck", "trk")[[j]], ":map:", name), identical(dim(observed), dims) &&
              identical(as.numeric(observed), as.numeric(expected)))
    }
  }
  roi_a <- csv(a$metadata$files$roi_summary_csv); roi_b <- csv(b$metadata$files$roi_summary_csv)
  roi_sub <- lapply(subsets, get_ftq("ftq_roi_summary"))
  predicted_roi <- roi_a
  numeric_roi <- names(roi_a)[vapply(roi_a, is.numeric, logical(1))]
  for (name in numeric_roi) predicted_roi[[name]] <- roi_a[[name]] + roi_sub[[2L]][[name]] - roi_sub[[1L]][[name]]
  predicted_roi$unique_exclusive_fate_sum_matches_focal <-
    predicted_roi$unique_exclusive_fate_sum == predicted_roi$unique_focal_intersecting_streamlines
  provenance_columns <- c("selected_track_space", "track_voxel_size_mm", "mapping_selection_rule")
  state$excluded_format_provenance_columns <<- provenance_columns
  compare_table(predicted_roi[setdiff(names(predicted_roi), provenance_columns)],
                roi_b[setdiff(names(roi_b), provenance_columns)], "roi-summary")
  compare_table(rsc_summary(predicted_core, percent = TRUE, warn_legacy = FALSE),
                csv(b$metadata$files$rsc_csv), "rsc-summary")
  compare_table(get_ftq("ftq_fate_summary")(predicted_core), csv(b$metadata$files$fate_csv), "fate-summary")

  endpoints <- lapply(subsets, ftq_endpoint_tables)
  endpoint_a <- csv(a$metadata$files$endpoint_summary_csv)
  endpoint_b <- csv(b$metadata$files$endpoint_summary_csv)
  check("endpoint_row_identity", identical(endpoint_a[c("target", "tolerance")], endpoint_b[c("target", "tolerance")]) &&
          identical(endpoint_a[c("target", "tolerance")], endpoints[[1L]]$summary[c("target", "tolerance")]) &&
          identical(endpoint_a[c("target", "tolerance")], endpoints[[2L]]$summary[c("target", "tolerance")]))
  predicted_endpoint <- endpoint_a
  numeric_endpoint <- names(endpoint_a)[vapply(endpoint_a, is.numeric, logical(1))]
  additive_endpoint <- setdiff(numeric_endpoint, grep("percent|^voxel_size_|_mask_voxels$|overlap_voxels$", numeric_endpoint, value = TRUE))
  for (name in additive_endpoint) predicted_endpoint[[name]] <- endpoint_a[[name]] +
    endpoints[[2L]]$summary[[name]] - endpoints[[1L]]$summary[[name]]
  ratios <- list(
    target_intersection_percent_of_unique_focal_intersecting_streamlines = c("target_intersection_streamlines", "unique_focal_intersecting_streamlines"),
    paired_percent_of_target_intersections = c("target_intersection_with_paired_endpoints", "target_intersection_streamlines"),
    paired_percent_of_focal_endpoint_streamlines = c("paired_endpoint_streamlines", "focal_endpoint_streamlines"),
    target_endpoint_percent_of_target_intersections = c("target_intersection_with_any_target_endpoint", "target_intersection_streamlines"),
    focal_endpoint_percent_of_target_intersections = c("target_intersection_with_any_focal_endpoint", "target_intersection_streamlines"),
    target_intersection_percent_of_focal_endpoint_streamlines = c("focal_endpoint_with_target_intersection_streamlines", "focal_endpoint_streamlines"),
    shared_mask_ambiguous_percent_of_target_intersections = c("shared_mask_ambiguous", "target_intersection_streamlines"),
    target_endpoint_only_percent_of_target_intersections = c("target_endpoint_only", "target_intersection_streamlines"),
    focal_endpoint_only_percent_of_target_intersections = c("focal_endpoint_only", "target_intersection_streamlines"),
    neither_relevant_endpoint_percent_of_target_intersections = c("neither_relevant_endpoint", "target_intersection_streamlines"))
  for (name in names(ratios)) predicted_endpoint[[name]] <- get_ftq("ftq_endpoint_percent")(
    predicted_endpoint[[ratios[[name]][[1L]]]], predicted_endpoint[[ratios[[name]][[2L]]]])
  compare_table(predicted_endpoint, endpoint_b, "endpoint-summary")

  penetration_a <- csv(a$metadata$files$endpoint_penetration_csv)
  penetration_b <- csv(b$metadata$files$endpoint_penetration_csv)
  penetration_rows <- list()
  for (target in names(targets)) {
    tables <- list(penetration_a, endpoints[[1L]]$penetration, endpoints[[2L]]$penetration)
    bins <- sort(unique(unlist(lapply(tables, function(x) x$unique_target_voxels[x$target == target]))))
    count_bin <- function(table, bin) sum(table$streamline_count[table$target == target & table$unique_target_voxels == bin])
    for (bin in bins) {
      count <- count_bin(tables[[1L]], bin) - count_bin(tables[[2L]], bin) + count_bin(tables[[3L]], bin)
      if (count == 0) next
      denominator <- predicted_roi[[paste0("unique_target_streamlines:", target)]]
      penetration_rows[[length(penetration_rows) + 1L]] <- data.frame(target = target,
        unique_target_voxels = bin, streamline_count = count,
        target_intersection_streamlines = denominator,
        percent_of_target_intersections = get_ftq("ftq_endpoint_percent")(count, denominator),
        streamline_voxelization_mode = "segment_interpolation",
        penetration_interpretation = "segment-interpolated target penetration")
    }
  }
  compare_table(do.call(rbind, penetration_rows), penetration_b, "endpoint-penetration")
  context_a <- csv(a$metadata$files$endpoint_context_csv); context_b <- csv(b$metadata$files$endpoint_context_csv)
  predicted_context <- context_a
  context_counts <- names(context_a)[vapply(context_a, is.numeric, logical(1))]
  context_counts <- setdiff(context_counts, grep("percent|mask_voxels", context_counts, value = TRUE))
  for (name in context_counts) predicted_context[[name]] <- context_a[[name]] +
    endpoints[[2L]]$context[[name]] - endpoints[[1L]]$context[[name]]
  predicted_context$endpoints_inside_percent <- get_ftq("ftq_endpoint_percent")(predicted_context$endpoints_inside,
                                                                              predicted_context$endpoint_count)
  predicted_context$streamlines_with_both_endpoints_inside_percent <- get_ftq("ftq_endpoint_percent")(
    predicted_context$streamlines_with_both_endpoints_inside, predicted_context$valid_streamlines)
  compare_table(predicted_context, context_b, "endpoint-context")
  check("immutable_prepared_inputs", identical(unname(tools::md5sum(paths)), unname(input_hashes)))
  check("immutable_result_evidence", identical(unname(tools::md5sum(result_files)), unname(result_hashes)))
  check("immutable_preparation_provenance", identical(unname(sha256(provenance_path)), state$preparation_provenance_sha256))
  state$accounting <<- list(density_voxels_changed = length(indices),
    density_net_occupancy_change = sum(delta),
    focal_voxel_rows_changed = length(unique(vapply(changed_core, function(x) x$fNode_Coords, character(1)))),
    all_primary_count_tables_and_maps_reconciled = TRUE,
    exact_format_parity = FALSE,
    derived_percentage_comparison = "Exact equality after the same CSV serialization used by package outputs; no count tolerance.")
}

write_state()
tryCatch({
  main()
  state$status <- "passed_causal_accounting"
  state$exact_format_parity <- FALSE
  state$completed_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  write_state()
  writeLines(capture.output(sessionInfo()), file.path(audit, "session-info.txt"))
  Sys.chmod(file.path(audit, "session-info.txt"), "0600")
  message("All observed differences accounted for; exact format parity remains FALSE: ", audit)
}, error = function(e) {
  state$status <<- "failed_causal_accounting"
  state$error <<- conditionMessage(e)
  write_state()
  stop(conditionMessage(e), call. = FALSE)
})
