ftq_endpoint_fixture <- function(overlap = FALSE) {
  blank <- array(0, dim = c(8, 4, 4))
  focal <- blank
  target <- blank
  context <- array(1, dim = dim(blank))
  focal[2, 2, 2] <- 1
  target[6, 2, 2] <- 1
  if (overlap) target[2, 2, 2] <- 1
  tracks <- list(
    rbind(c(2, 2, 2), c(6, 2, 2)),
    rbind(c(1, 2, 2), c(2, 2, 2), c(6, 2, 2), c(7, 2, 2)),
    rbind(c(1, 2, 2), c(2, 2, 2), c(6, 2, 2)),
    rbind(c(2, 2, 2), c(6, 2, 2), c(7, 2, 2)),
    matrix(c(2, 2, 2), nrow = 1L)
  )
  list(focal = focal, target = target, context = context, tracks = tracks)
}

ftq_endpoint_write_trk <- function(path, tracks, dims) {
  write_fixed <- function(con, text, n) {
    bytes <- rep(as.raw(0), n)
    raw <- charToRaw(text)
    if (length(raw)) bytes[seq_len(min(length(raw), n))] <- raw[seq_len(min(length(raw), n))]
    writeBin(bytes, con)
  }
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  endian <- "little"
  write_fixed(con, "TRACK", 6L)
  writeBin(as.integer(dims), con, size = 2L, endian = endian)
  writeBin(as.numeric(c(1, 1, 1)), con, size = 4L, endian = endian)
  writeBin(as.numeric(c(0, 0, 0)), con, size = 4L, endian = endian)
  writeBin(as.integer(0), con, size = 2L, endian = endian)
  write_fixed(con, "", 200L)
  writeBin(as.integer(0), con, size = 2L, endian = endian)
  write_fixed(con, "", 200L)
  writeBin(as.numeric(t(diag(4))), con, size = 4L, endian = endian)
  write_fixed(con, "", 444L)
  write_fixed(con, "RAS", 4L)
  write_fixed(con, "", 4L)
  writeBin(as.numeric(c(1, 0, 0, 0, 1, 0)), con, size = 4L, endian = endian)
  write_fixed(con, "", 2L)
  writeBin(as.integer(rep(0, 6)), con, size = 1L, endian = endian)
  writeBin(as.integer(length(tracks)), con, size = 4L, endian = endian)
  writeBin(as.integer(2), con, size = 4L, endian = endian)
  writeBin(as.integer(1000), con, size = 4L, endian = endian)
  for (track in tracks) {
    coords <- as.matrix(track)[, 1:3, drop = FALSE]
    writeBin(as.integer(nrow(coords)), con, size = 4L, endian = endian)
    writeBin(as.numeric(t(coords)), con, size = 4L, endian = endian)
  }
  invisible(path)
}

ftq_run_endpoint_fixture <- function(backend = "r", reverse = FALSE, overlap = FALSE,
                                     endpoint_analysis = "exact_and_dilated1") {
  fixture <- ftq_endpoint_fixture(overlap = overlap)
  FibrTraQ(
    fNode = fixture$focal,
    localROI = NULL,
    cNode_list = list(target = fixture$target),
    cNode_labels = "target",
    all_tracks = fixture$tracks,
    discretization = "round",
    focal_assignment_mode = "all_focal_voxels",
    streamline_voxelization_mode = "sample_points",
    reverse_streamlines = reverse,
    track_space = "voxel",
    axis_flip = c(FALSE, FALSE, FALSE),
    backend = backend,
    endpoint_analysis = endpoint_analysis,
    endpoint_context_list = if (identical(endpoint_analysis, "none")) NULL else list(brain = fixture$context)
  )
}

test_that("endpoint analysis partitions target intersections and matches both backends", {
  r_result <- ftq_run_endpoint_fixture("r")
  cpp_result <- ftq_run_endpoint_fixture("cpp")
  r_endpoint <- ftq_endpoint_tables(r_result)
  cpp_endpoint <- ftq_endpoint_tables(cpp_result)

  expect_equal(r_result, cpp_result, ignore_attr = TRUE)
  expect_equal(r_endpoint, cpp_endpoint, ignore_attr = TRUE)
  expect_equal(r_endpoint$summary$tolerance, c("exact", "dilated1"))

  exact <- r_endpoint$summary[r_endpoint$summary$tolerance == "exact", ]
  expect_equal(exact$target_intersection_streamlines, 4)
  expect_equal(exact$paired_focal_target_endpoints, 1)
  expect_equal(exact$target_endpoint_only, 1)
  expect_equal(exact$focal_endpoint_only, 1)
  expect_equal(exact$neither_relevant_endpoint, 1)
  expect_equal(exact$shared_mask_ambiguous, 0)
  expect_equal(
    sum(exact[c(
      "shared_mask_ambiguous", "paired_focal_target_endpoints",
      "target_endpoint_only", "focal_endpoint_only", "neither_relevant_endpoint"
    )]),
    exact$target_intersection_streamlines
  )
  expect_equal(sum(r_endpoint$penetration$streamline_count), 4)
  expect_equal(r_endpoint$context$context, c("reference_grid", "brain"))
  expect_equal(r_endpoint$context$endpoints_outside, c(0, 0))
})

test_that("endpoint classifications are reversal invariant and shared masks are explicit", {
  forward <- ftq_endpoint_tables(ftq_run_endpoint_fixture("cpp", reverse = FALSE))
  reverse <- ftq_endpoint_tables(ftq_run_endpoint_fixture("cpp", reverse = TRUE))
  expect_equal(forward$summary, reverse$summary)
  expect_equal(forward$penetration, reverse$penetration)
  expect_equal(forward$context, reverse$context)

  overlap <- ftq_endpoint_tables(ftq_run_endpoint_fixture("cpp", overlap = TRUE))$summary
  exact <- overlap[overlap$tolerance == "exact", ]
  expect_gt(exact$shared_mask_ambiguous, 0)
  expect_gte(
    exact$target_intersection_with_paired_endpoints,
    exact$paired_focal_target_endpoints
  )
  expect_lte(
    exact$target_intersection_with_paired_endpoints -
      exact$paired_focal_target_endpoints,
    exact$shared_mask_ambiguous
  )
  expect_equal(
    exact$paired_focal_target_endpoints + exact$shared_mask_ambiguous +
      exact$target_endpoint_only + exact$focal_endpoint_only +
      exact$neither_relevant_endpoint,
    exact$target_intersection_streamlines
  )
})

test_that("all paired endpoint memberships are monotonic across dilation", {
  summary <- ftq_endpoint_tables(ftq_run_endpoint_fixture("cpp", overlap = TRUE))$summary
  exact <- summary[summary$tolerance == "exact", ]
  dilated1 <- summary[summary$tolerance == "dilated1", ]
  expect_lte(exact$paired_endpoint_streamlines, dilated1$paired_endpoint_streamlines)
  expect_lte(
    exact$target_intersection_with_paired_endpoints,
    dilated1$target_intersection_with_paired_endpoints
  )
})

test_that("compiled endpoint dilation is clipped 26-neighbor dilation", {
  dims <- c(4L, 3L, 2L)
  mask <- array(FALSE, dim = dims)
  mask[1, 1, 1] <- TRUE
  mask[4, 3, 2] <- TRUE
  observed <- array(
    FibrTraQ:::ftq_endpoint_dilate_lut(as.vector(mask), dims),
    dim = dims
  )
  expected <- array(FALSE, dim = dims)
  expected[1:2, 1:2, 1:2] <- TRUE
  expected[3:4, 2:3, 1:2] <- TRUE
  expect_identical(observed, expected)
})

test_that("endpoint R and C++ parity covers discretization and voxelization modes", {
  fixture <- ftq_endpoint_fixture()
  for (discretization in c("round", "floor", "ceil", "any")) {
    for (voxelization in c("sample_points", "segment_interpolation")) {
      args <- list(
        fNode = fixture$focal,
        localROI = NULL,
        cNode_list = list(target = fixture$target),
        cNode_labels = "target",
        all_tracks = fixture$tracks,
        discretization = discretization,
        streamline_voxelization_mode = voxelization,
        track_space = "voxel",
        axis_flip = c(FALSE, FALSE, FALSE),
        endpoint_context_list = list(brain = fixture$context)
      )
      r_result <- do.call(FibrTraQ, c(args, list(backend = "r")))
      cpp_result <- do.call(FibrTraQ, c(args, list(backend = "cpp")))
      expect_equal(r_result, cpp_result, ignore_attr = TRUE)
      expect_equal(
        ftq_endpoint_tables(r_result),
        ftq_endpoint_tables(cpp_result),
        ignore_attr = TRUE
      )
    }
  }

  flipped_focal <- array(0, dim = c(8, 4, 4))
  flipped_target <- flipped_focal
  flipped_focal[2, 2, 2] <- 1
  flipped_target[6, 2, 2] <- 1
  flipped_tracks <- list(rbind(c(7, 2, 2), c(3, 2, 2)))
  flipped_args <- list(
    fNode = flipped_focal,
    localROI = NULL,
    cNode_list = list(target = flipped_target),
    cNode_labels = "target",
    all_tracks = flipped_tracks,
    discretization = "round",
    track_space = "voxel",
    axis_flip = c(TRUE, FALSE, FALSE),
    endpoint_context_list = list(brain = array(1, dim = c(8, 4, 4)))
  )
  r_flipped <- do.call(FibrTraQ, c(flipped_args, list(backend = "r")))
  cpp_flipped <- do.call(FibrTraQ, c(flipped_args, list(backend = "cpp")))
  expect_equal(
    ftq_endpoint_tables(r_flipped),
    ftq_endpoint_tables(cpp_flipped),
    ignore_attr = TRUE
  )
  expect_equal(
    ftq_endpoint_tables(cpp_flipped)$summary$paired_focal_target_endpoints[1],
    1
  )
})

test_that("segment interpolation changes penetration without changing endpoint membership", {
  focal <- array(0, dim = c(8, 3, 3))
  target <- focal
  focal[2, 2, 2] <- 1
  target[4:6, 2, 2] <- 1
  track <- list(rbind(c(2, 2, 2), c(6, 2, 2)))
  sample <- FibrTraQ(
    focal, NULL, list(target), "target", track,
    discretization = "round", streamline_voxelization_mode = "sample_points",
    track_space = "voxel", axis_flip = c(FALSE, FALSE, FALSE), backend = "cpp"
  )
  interpolated <- FibrTraQ(
    focal, NULL, list(target), "target", track,
    discretization = "round", streamline_voxelization_mode = "segment_interpolation",
    track_space = "voxel", axis_flip = c(FALSE, FALSE,FALSE), backend = "cpp"
  )
  expect_equal(
    ftq_endpoint_tables(sample)$summary$paired_focal_target_endpoints,
    ftq_endpoint_tables(interpolated)$summary$paired_focal_target_endpoints
  )
  expect_equal(ftq_endpoint_tables(sample)$penetration$unique_target_voxels, 1)
  expect_equal(ftq_endpoint_tables(interpolated)$penetration$unique_target_voxels, 3)
})

test_that("streamed endpoint tables are chunk-size invariant", {
  fixture <- ftq_endpoint_fixture()
  path <- tempfile(fileext = ".trk")
  ftq_endpoint_write_trk(path, fixture$tracks, dim(fixture$focal))
  run_stream <- function(chunk_size) {
    FibrTraQ(
      fixture$focal,
      NULL,
      list(fixture$target),
      "target",
      ftq_trk_stream(path, chunk_size = chunk_size),
      discretization = "round",
      track_space = "voxel",
      axis_flip = c(FALSE, FALSE, FALSE),
      backend = "cpp",
      endpoint_context_list = list(brain = fixture$context)
    )
  }
  chunk_one <- run_stream(1L)
  chunk_three <- run_stream(3L)
  expect_equal(chunk_one, chunk_three, ignore_attr = TRUE)
  expect_equal(
    ftq_endpoint_tables(chunk_one),
    ftq_endpoint_tables(chunk_three),
    ignore_attr = TRUE
  )
})

test_that("streamed local-search audits are chunk-size invariant", {
  dims <- c(12L, 2L, 2L)
  focal <- array(0, dim = dims)
  focal[2, 1, 1] <- 1
  target <- array(0, dim = dims)
  target[10, 1, 1] <- 1
  attr(focal, "fibrtraq_voxel_to_world_basis") <- diag(3)
  tracks <- list(
    rbind(c(2, 1, 1), c(10, 1, 1)),
    rbind(c(1.4, 1, 1), c(9.4, 1, 1)),
    rbind(c(2, 1, 1), c(7, 1, 1))
  )
  path <- tempfile(fileext = ".trk")
  ftq_endpoint_write_trk(path, tracks, dims)

  run_stream <- function(chunk_size) {
    FibrTraQ(
      focal,
      array(0, dim = dims),
      list(target),
      "target",
      ftq_trk_stream(path, chunk_size = chunk_size),
      discretization = "round",
      streamline_voxelization_mode = "segment_interpolation",
      track_space = "voxel",
      axis_flip = c(FALSE, FALSE, FALSE),
      assignment_mode = "local_search",
      search_radius_mm = 0.75,
      endpoint_analysis = "none",
      backend = "cpp"
    )
  }
  chunk_one <- run_stream(1L)
  chunk_three <- run_stream(3L)
  roi_one <- attr(chunk_one, "fibrtraq_roi_summary")
  roi_three <- attr(chunk_three, "fibrtraq_roi_summary")
  roi_one$backend <- NULL
  roi_three$backend <- NULL

  expect_equal(chunk_one, chunk_three, ignore_attr = TRUE)
  expect_equal(roi_one, roi_three)
  expect_equal(
    ftq_local_search_tables(chunk_one),
    ftq_local_search_tables(chunk_three),
    tolerance = 1e-10
  )
})

test_that("manifest endpoint contexts are parsed, preflighted, and written", {
  skip_if_not_installed("RNifti")
  fixture <- ftq_endpoint_fixture()
  root <- file.path(tempdir(), paste0("endpoint-manifest-", sample.int(1e8, 1)))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  fnode <- file.path(root, "focal.nii.gz")
  target <- file.path(root, "target.nii.gz")
  context <- file.path(root, "brain.nii.gz")
  tractogram <- file.path(root, "tracks.trk")
  RNifti::writeNifti(fixture$focal, fnode)
  RNifti::writeNifti(fixture$target, target)
  RNifti::writeNifti(fixture$context, context)
  ftq_endpoint_write_trk(tractogram, fixture$tracks, dim(fixture$focal))
  manifest <- list(
    dataset_root = root,
    subject = "sub-endpoint",
    description = "endpoint-manifest",
    fnode = basename(fnode),
    tractogram = basename(tractogram),
    targets = list(target = basename(target)),
    endpoint_contexts = list(brain_and_cord = basename(context)),
    coordinate = list(
      track_space = "voxel",
      discretization = "round",
      axis_flip = c(FALSE, FALSE, FALSE)
    ),
    analysis = list(endpoint_analysis = "exact_and_dilated1"),
    execution = list(track_reader = "stream", track_chunk_size = 2L),
    outputs = list(
      output_dir = file.path(root, "outputs"),
      write_maps = FALSE,
      write_report = FALSE,
      overwrite = TRUE
    )
  )
  parsed <- ftq_read_manifest(manifest)
  expect_equal(parsed$endpoint_context_labels, "brain_and_cord")
  expect_equal(parsed$endpoint_analysis, "exact_and_dilated1")
  expect_true(ftq_preflight(parsed)$passed)
  run <- ftq_run(parsed, backend = "cpp")
  expect_true(all(file.exists(c(
    run$outputs$endpoint_summary_csv,
    run$outputs$endpoint_penetration_csv,
    run$outputs$endpoint_context_csv
  ))))
  expect_equal(run$validation$biological_qc_status, "pass")

  manifest$analysis$endpoint_analysis <- "none"
  expect_error(ftq_read_manifest(manifest), "endpoint_contexts require")
})

test_that("endpoint-off mode preserves existing outputs and attributes", {
  enabled <- ftq_run_endpoint_fixture("cpp")
  disabled <- ftq_run_endpoint_fixture("cpp", endpoint_analysis = "none")
  expect_equal(enabled, disabled, ignore_attr = TRUE)

  legacy_attributes <- c(
    "fibrtraq_backend", "fibrtraq_focal_assignment_mode",
    "fibrtraq_streamline_voxelization_mode", "fibrtraq_target_intersection_mode",
    "fibrtraq_local_roi_defaulted_to_focal", "fibrtraq_track_space",
    "fibrtraq_track_voxel_size", "fibrtraq_axis_flip", "fibrtraq_discretization",
    "fibrtraq_reverse_streamlines", "fibrtraq_roi_summary",
    "fibrtraq_mapping_scores", "fibrtraq_mapping_selection_rule"
  )
  for (attribute in legacy_attributes) {
    expect_equal(attr(enabled, attribute), attr(disabled, attribute))
  }
  expect_null(attr(disabled, "fibrtraq_endpoint_summary"))
  expect_error(ftq_endpoint_tables(disabled), "unavailable")
})

test_that("endpoint context masks are strict and boundary distances are classified", {
  fixture <- ftq_endpoint_fixture()
  bad_values <- fixture$context
  bad_values[1] <- 2
  expect_error(
    FibrTraQ(
      fixture$focal, NULL, list(fixture$target), "target", fixture$tracks,
      discretization = "round", track_space = "voxel",
      axis_flip = c(FALSE, FALSE, FALSE), endpoint_context_list = list(bad = bad_values)
    ),
    "binary"
  )
  expect_error(
    FibrTraQ(
      fixture$focal, NULL, list(fixture$target), "target", fixture$tracks,
      discretization = "round", track_space = "voxel",
      axis_flip = c(FALSE, FALSE, FALSE),
      endpoint_context_list = list(bad = array(1, c(2, 2, 2)))
    ),
    "dimensions"
  )

  outside_tracks <- list(
    rbind(c(0, 2, 2), c(2, 2, 2), c(6, 2, 2)),
    rbind(c(-2, 2, 2), c(2, 2, 2), c(6, 2, 2))
  )
  result <- FibrTraQ(
    fixture$focal, NULL, list(fixture$target), "target", outside_tracks,
    discretization = "round", track_space = "voxel",
    axis_flip = c(FALSE, FALSE, FALSE), backend = "cpp",
    endpoint_context_list = list(brain = fixture$context)
  )
  grid <- ftq_endpoint_tables(result)$context[1, ]
  expect_equal(grid$endpoints_one_voxel_outside, 1)
  expect_equal(grid$endpoints_more_than_one_voxel_outside, 1)
  expect_equal(grid$x_below_grid, 2)
})

test_that("endpoint context records all boundary directions and one-point tracks", {
  dims <- c(5, 5, 5)
  focal <- array(0, dim = dims)
  target <- focal
  focal[3, 3, 3] <- 1
  target[4, 3, 3] <- 1
  tracks <- list(
    matrix(c(0, 3, 3), nrow = 1),
    matrix(c(6, 3, 3), nrow = 1),
    matrix(c(3, 0, 3), nrow = 1),
    matrix(c(3, 6, 3), nrow = 1),
    matrix(c(3, 3, 0), nrow = 1),
    matrix(c(3, 3, 6), nrow = 1),
    matrix(c(3, 3, 3), nrow = 1)
  )
  result <- FibrTraQ(
    focal, NULL, list(target), "target", tracks,
    discretization = "round", track_space = "voxel",
    axis_flip = c(FALSE, FALSE, FALSE), backend = "cpp"
  )
  endpoint <- ftq_endpoint_tables(result)
  context <- endpoint$context
  grid <- context[context$context == "reference_grid", , drop = FALSE]
  expect_equal(grid$endpoint_count, 7)
  expect_equal(grid$endpoints_one_voxel_outside, 6)
  expect_equal(grid$endpoints_more_than_one_voxel_outside, 0)
  expect_equal(
    unlist(grid[c(
      "x_below_grid", "x_above_grid", "y_below_grid",
      "y_above_grid", "z_below_grid", "z_above_grid"
    )], use.names = FALSE),
    rep(1, 6)
  )
  exact <- endpoint$summary[endpoint$summary$tolerance == "exact", ]
  expect_equal(exact$target_intersection_streamlines, 0)
  expect_equal(exact$focal_endpoint_streamlines, 1)
  expect_equal(exact$paired_endpoint_streamlines, 0)
  expect_equal(nrow(endpoint$penetration), 0)
})

test_that("endpoint analysis keeps target-specific classes for multiple targets", {
  dims <- c(8, 3, 3)
  focal <- array(0, dim = dims)
  target_a <- focal
  target_b <- focal
  focal[2, 2, 2] <- 1
  target_a[5, 2, 2] <- 1
  target_b[7, 2, 2] <- 1
  track <- list(rbind(c(2, 2, 2), c(5, 2, 2), c(7, 2, 2)))
  endpoint <- ftq_endpoint_tables(FibrTraQ(
    focal, NULL, list(target_a, target_b), c("target_a", "target_b"), track,
    discretization = "round", streamline_voxelization_mode = "segment_interpolation",
    track_space = "voxel", axis_flip = c(FALSE, FALSE, FALSE), backend = "cpp"
  ))
  exact <- endpoint$summary[endpoint$summary$tolerance == "exact", ]
  expect_equal(exact$target_intersection_streamlines, c(1, 1))
  expect_equal(exact$focal_endpoint_only, c(1, 0))
  expect_equal(exact$paired_focal_target_endpoints, c(0, 1))
  expect_setequal(endpoint$penetration$target, c("target_a", "target_b"))
})

test_that("endpoint R and C++ parity covers supported coordinate spaces", {
  fixture <- ftq_endpoint_fixture()
  reference <- RNifti::asNifti(fixture$focal)
  desired <- fixture$tracks
  voxel_size <- c(2, 3, 4)
  cases <- list(
    voxel = list(tracks = desired),
    vox10 = list(tracks = lapply(desired, function(track) track / 10)),
    voxelmm = list(
      tracks = lapply(desired, function(track) sweep(track - 0.5, 2, voxel_size, "*")),
      track_voxel_size = voxel_size
    ),
    world = list(
      tracks = lapply(desired, function(track) {
        matrix(RNifti::voxelToWorld(track, reference, simple = FALSE), ncol = 3)
      })
    ),
    trackvis = list(
      tracks = lapply(desired, function(track) {
        matrix(RNifti::voxelToWorld(track, reference, simple = FALSE), ncol = 3)
      }),
      track_affine = diag(4)
    )
  )
  for (space in names(cases)) {
    args <- list(
      fNode = fixture$focal,
      localROI = NULL,
      cNode_list = list(target = fixture$target),
      cNode_labels = "target",
      all_tracks = cases[[space]]$tracks,
      discretization = "round",
      track_space = space,
      ref_image = reference,
      axis_flip = c(FALSE, FALSE, FALSE),
      endpoint_context_list = list(brain = fixture$context)
    )
    if (!is.null(cases[[space]]$track_voxel_size)) {
      args$track_voxel_size <- cases[[space]]$track_voxel_size
    }
    if (!is.null(cases[[space]]$track_affine)) {
      args$track_affine <- cases[[space]]$track_affine
    }
    r_result <- do.call(FibrTraQ, c(args, list(backend = "r")))
    cpp_result <- do.call(FibrTraQ, c(args, list(backend = "cpp")))
    expect_equal(r_result, cpp_result, ignore_attr = TRUE, info = space)
    expect_equal(
      ftq_endpoint_tables(r_result),
      ftq_endpoint_tables(cpp_result),
      ignore_attr = TRUE,
      info = space
    )
  }
})

test_that("endpoint files, validation, reports, and statistics are additive", {
  result <- ftq_run_endpoint_fixture("cpp")
  fixture <- ftq_endpoint_fixture()
  output_dir <- file.path(tempdir(), paste0("endpoint-output-", sample.int(1e8, 1)))
  report_dir <- file.path(output_dir, "report")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rsc <- rsc_summary(result, percent = TRUE, warn_legacy = FALSE)
  outputs <- ftq_write_outputs(
    result,
    rsc,
    template_image = fixture$focal,
    output_dir = output_dir,
    prefix = "endpoint-test",
    write_maps = FALSE
  )
  expect_true(all(file.exists(c(
    outputs$endpoint_summary_csv,
    outputs$endpoint_penetration_csv,
    outputs$endpoint_context_csv
  ))))

  validation <- ftq_validate_run(
    outputs$result_csv,
    outputs$rsc_csv,
    roi_summary_csv = outputs$roi_summary_csv,
    endpoint_summary = outputs$endpoint_summary_csv,
    endpoint_penetration = outputs$endpoint_penetration_csv,
    endpoint_context = outputs$endpoint_context_csv
  )
  expect_true(validation$passed)
  expect_equal(validation$biological_qc_status, "pass")

  report <- ftq_report(
    outputs$result_csv,
    outputs$rsc_csv,
    output_dir = report_dir,
    endpoint_summary = outputs$endpoint_summary_csv,
    endpoint_penetration = outputs$endpoint_penetration_csv,
    endpoint_context = outputs$endpoint_context_csv
  )
  html <- paste(readLines(report$html, warn = FALSE), collapse = "\n")
  expect_match(html, "Full-Path Intersections vs End-Voxel Connectivity", fixed = TRUE)
  expect_match(html, "Different estimands", fixed = TRUE)
  expect_match(
    html,
    "<div class='table-scroll' role='region' tabindex='0'>",
    fixed = TRUE
  )
  expect_match(
    html,
    "<meta name='viewport' content='width=device-width,initial-scale=1'>",
    fixed = TRUE
  )
  expect_true(all(file.exists(unlist(report$figures[grepl("endpoint", names(report$figures))]))))

  legacy_report <- ftq_report(
    outputs$result_csv,
    outputs$rsc_csv,
    output_dir = file.path(output_dir, "legacy-report")
  )
  legacy_html <- paste(readLines(legacy_report$html, warn = FALSE), collapse = "\n")
  expect_false(grepl("Full-Path Intersections vs End-Voxel Connectivity", legacy_html, fixed = TRUE))

  old_subject <- ftq_stat_subject_summary(result)
  expect_identical(names(old_subject), c("roi", "targets", "voxels"))
  endpoint_subject <- ftq_stat_subject_summary(
    result,
    endpoint_summary = outputs$endpoint_summary_csv,
    endpoint_penetration = outputs$endpoint_penetration_csv
  )
  expect_true(all(c("endpoints", "penetration") %in% names(endpoint_subject)))

  old_targets <- ftq_target_summary_tables(result, attr(result, "fibrtraq_roi_summary"))
  expect_identical(names(old_targets), c("unique_streamlines", "voxel_associations"))
  endpoint_targets <- ftq_target_summary_tables(
    result,
    attr(result, "fibrtraq_roi_summary"),
    endpoint_summary = outputs$endpoint_summary_csv
  )
  expect_true("endpoint_streamlines" %in% names(endpoint_targets))
})

test_that("penetration percentages are calculated per sparse histogram bin", {
  expect_equal(
    FibrTraQ:::ftq_endpoint_percent(c(1, 2, 7), 10),
    c(10, 20, 70)
  )
  expect_equal(
    FibrTraQ:::ftq_endpoint_percent(c(1, 2), c(4, 8)),
    c(25, 25)
  )
  expect_true(all(is.na(FibrTraQ:::ftq_endpoint_percent(c(1, 2), 0))))

  penetration <- data.frame(
    target = rep("T", 5),
    unique_target_voxels = c(1, 2, 3, 7, 10),
    streamline_count = c(2, 3, 4, 5, 6),
    target_intersection_streamlines = rep(20, 5),
    stringsAsFactors = FALSE
  )
  summary <- FibrTraQ:::ftq_endpoint_penetration_summary(penetration)
  exclusive_counts <- unlist(summary[c(
    "one_voxel_touch_streamlines",
    "two_voxel_only_streamlines",
    "three_to_four_voxels_streamlines",
    "five_to_nine_voxels_streamlines",
    "ten_or_more_voxels_streamlines"
  )], use.names = FALSE)
  exclusive_percent <- unlist(summary[c(
    "one_voxel_touch_percent_of_target_intersections",
    "two_voxel_only_percent_of_target_intersections",
    "three_to_four_voxels_percent_of_target_intersections",
    "five_to_nine_voxels_percent_of_target_intersections",
    "ten_or_more_voxels_percent_of_target_intersections"
  )], use.names = FALSE)
  expect_equal(exclusive_counts, c(2, 3, 4, 5, 6))
  expect_equal(sum(exclusive_counts), 20)
  expect_equal(sum(exclusive_percent), 100)
})
