library(testthat)
library(FibrTraQ)

ftq_write_test_trk <- function(path,
                               tracks,
                               dims = c(4L, 4L, 4L),
                               voxel_order = "RAS",
                               voxel_size = c(1, 1, 1)) {
  write_fixed <- function(con, text, n) {
    bytes <- rep(as.raw(0), n)
    text_bytes <- charToRaw(text)
    keep <- seq_len(min(length(text_bytes), n))
    bytes[keep] <- text_bytes[keep]
    writeBin(bytes, con)
  }

  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  endian <- "little"
  write_fixed(con, "TRACK", 6L)
  writeBin(as.integer(dims), con, size = 2L, endian = endian)
  writeBin(as.numeric(voxel_size), con, size = 4L, endian = endian)
  writeBin(as.numeric(c(0, 0, 0)), con, size = 4L, endian = endian)
  writeBin(as.integer(0), con, size = 2L, endian = endian)
  write_fixed(con, "", 200L)
  writeBin(as.integer(0), con, size = 2L, endian = endian)
  write_fixed(con, "", 200L)
  writeBin(as.numeric(t(diag(4))), con, size = 4L, endian = endian)
  write_fixed(con, "", 444L)
  write_fixed(con, voxel_order, 4L)
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

ftq_test_path <- function(path) {
  gsub("\\\\", "/", path)
}

test_that("workflow output writer validates CSVs, maps, and report", {
  skip_if_not_installed("RNifti")

  tmp <- tempdir()
  template <- file.path(tmp, "template.nii.gz")
  RNifti::writeNifti(array(0, dim = c(3, 3, 3)), template)
  background <- file.path(tmp, "background.nii.gz")
  RNifti::writeNifti(array(seq_len(27), dim = c(3, 3, 3)), background)

  result <- data.frame(
    fNode_Coords = c("1_1_1", "2_2_2"),
    tracks = c(3L, 0L),
    local = c(1L, 0L),
    leave_no_target = c(1L, 0L),
    T = c(1L, 0L),
    none = c(0L, 1L),
    stringsAsFactors = FALSE
  )
  rsc <- rsc_summary(result, percent = TRUE, warn_legacy = FALSE)
  roi_summary <- data.frame(
    unique_focal_intersecting_streamlines = 3,
    unique_local_streamlines = 1,
    unique_leave_no_target_streamlines = 1,
    unique_leave_with_target_streamlines = 1,
    `unique_target_streamlines:T` = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  out <- ftq_write_outputs(
    result = result,
    rsc = rsc,
    template_image = template,
    output_dir = file.path(tmp, "out"),
    prefix = "sub-001_ses-01_desc-test",
    target_labels = "T",
    write_maps = TRUE,
    overwrite = TRUE
  )

  expect_true(file.exists(out$result_csv))
  expect_true(file.exists(out$rsc_csv))
  expect_true(file.exists(out$roi_summary_csv))
  expect_true(file.exists(out$map_paths$tracks))
  expect_true(file.exists(out$map_paths[["leave-no-target"]]))
  expect_true(file.exists(out$map_paths[["target-T"]]))

  validation <- ftq_validate_run(
    out$result_csv,
    out$rsc_csv,
    out$map_paths,
    expected = list(
      result_csv = out$result_csv,
      rsc_csv = out$rsc_csv,
      map_paths = out$map_paths[c("tracks", "leave-no-target")]
    )
  )
  expect_true(validation$passed)
  expect_true("expected_result_values_match" %in% validation$checks$check)

  report <- ftq_report(
    result_csv = out$result_csv,
    rsc_csv = out$rsc_csv,
    manifest = list(
      prefix = "sub-001_ses-01_desc-test",
      subject = "sub-001",
      description = "test",
      focal_assignment_mode = "all_focal_voxels",
      target_intersection_mode = "full_streamline",
      background_image = background
    ),
    output_dir = file.path(tmp, "report"),
    map_paths = out$map_paths,
    roi_summary = roi_summary,
    validation = validation,
    preflight = list(
      passed = TRUE,
      checks = data.frame(
        check = "background_image_dims_match_fnode",
        severity = "warning",
        passed = TRUE,
        detail = "background = 3x3x3 fnode = 3x3x3",
        stringsAsFactors = FALSE
      )
    )
  )
  expect_true(file.exists(report$html))
  expect_true(any(file.exists(unlist(report$figures))))
  expect_true(any(file.exists(unlist(report$tables))))
  connection_table <- utils::read.csv(
    report$tables$connection_summary,
    check.names = FALSE
  )
  expect_identical(connection_table$cNode, c("local", "T"))
  expect_equal(connection_table$unique_streamline_connections, c(1, 1))
  expect_equal(connection_table$connected_fNode_voxels, c(1, 1))
  fate_table <- utils::read.csv(report$tables$fate_summary)
  expect_equal(sum(fate_table$count), 3)
  expect_true(all(
    fate_table$summary_basis == "unique_focal_intersecting_streamlines"
  ))
  report_html <- paste(readLines(report$html, warn = FALSE), collapse = "\n")
  expect_match(report_html, "Unique streamline connections", fixed = TRUE)
  expect_match(report_html, "fNode streamline density heatmap", fixed = TRUE)
  expect_match(report_html, "Full preflight details", fixed = TRUE)
  expect_match(report_html, "automated computational checks passed", fixed = TRUE)
  expect_match(report_html, "full-streamline intersections", fixed = TRUE)
  expect_false(grepl("This is not a 1 mm local-search run", report_html, fixed = TRUE))
  expect_false(grepl("Local-Search Connection Audit", report_html, fixed = TRUE))

  nested_report <- ftq_report(
    result_csv = out$result_csv,
    rsc_csv = out$rsc_csv,
    manifest = list(
      prefix = "sub-001_ses-01_desc-nested-report",
      subject = "sub-001",
      description = "nested-report",
      analysis = list(
        focal_assignment_mode = "all_focal_voxels",
        streamline_voxelization_mode = "segment_interpolation",
        target_intersection_mode = "local_search",
        local_search_radius_mm = 9,
        endpoint_analysis = "none"
      )
    ),
    output_dir = file.path(tmp, "nested-report"),
    map_paths = out$map_paths,
    overwrite = TRUE
  )
  nested_html <- paste(
    readLines(nested_report$html, warn = FALSE),
    collapse = "\n"
  )
  expect_match(nested_html, "Target assignment mode:</b> local_search",
               fixed = TRUE)
  expect_match(nested_html, "Local search radius (mm):</b> 9",
               fixed = TRUE)

  relative_manifest <- ftq_read_manifest(list(
    dataset_root = tmp,
    subject = "sub-001",
    description = "relative-background",
    fnode = basename(template),
    tractogram = "unused.trk",
    targets = list(T = basename(template)),
    report = list(background_image = basename(background)),
    outputs = list(output_dir = file.path(tmp, "relative-background"))
  ))
  expect_equal(
    relative_manifest$background_image,
    normalizePath(background, mustWork = TRUE)
  )
  relative_report <- ftq_report(
    result_csv = out$result_csv,
    rsc_csv = out$rsc_csv,
    manifest = relative_manifest,
    output_dir = file.path(tmp, "relative-background-report"),
    map_paths = out$map_paths,
    overwrite = TRUE
  )
  relative_html <- paste(
    readLines(relative_report$html, warn = FALSE),
    collapse = "\n"
  )
  expect_match(
    relative_html,
    normalizePath(background, mustWork = TRUE),
    fixed = TRUE
  )

  bad_result <- result
  bad_result$T[1] <- bad_result$T[1] + 1L
  bad_result_csv <- file.path(tmp, "bad-result.csv")
  utils::write.csv(bad_result, bad_result_csv, row.names = FALSE)
  bad_validation <- ftq_validate_run(
    out$result_csv,
    out$rsc_csv,
    expected = list(result_csv = bad_result_csv)
  )
  expect_false(bad_validation$passed)
  expect_false(bad_validation$checks$passed[bad_validation$checks$check == "expected_result_values_match"])
})

test_that("local-search audit is written, validated, and explained in reports", {
  skip_if_not_installed("RNifti")

  tmp <- file.path(tempdir(), paste0("local-search-audit-", sample.int(1e8, 1)))
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  dims <- c(12L, 2L, 2L)
  focal <- array(0, dim = dims)
  focal[2, 1, 1] <- 1
  target <- array(0, dim = dims)
  target[10, 1, 1] <- 1
  attr(focal, "fibrtraq_voxel_to_world_basis") <- diag(3)
  result <- FibrTraQ(
    fNode = focal,
    localROI = array(0, dim = dims),
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = list(
      rbind(c(2, 1, 1), c(10, 1, 1)),
      rbind(c(1.4, 1, 1), c(9.4, 1, 1)),
      rbind(c(2, 1, 1), c(7, 1, 1))
    ),
    discretization = "round",
    streamline_voxelization_mode = "segment_interpolation",
    track_space = "voxel",
    axis_flip = c(FALSE, FALSE, FALSE),
    assignment_mode = "local_search",
    search_radius_mm = 0.75,
    endpoint_analysis = "none",
    backend = "cpp"
  )
  rsc <- rsc_summary(result, percent = TRUE, warn_legacy = FALSE)
  template <- file.path(tmp, "template.nii.gz")
  RNifti::writeNifti(array(0, dim = dims), template)
  out <- ftq_write_outputs(
    result = result,
    rsc = rsc,
    template_image = template,
    output_dir = file.path(tmp, "outputs"),
    prefix = "sub-test_desc-local-search",
    target_labels = "T",
    write_maps = FALSE,
    overwrite = TRUE
  )

  expect_true(file.exists(out$local_search_assignments_csv))
  assignments <- utils::read.csv(
    out$local_search_assignments_csv,
    check.names = FALSE
  )
  expect_equal(nrow(assignments), 3L)
  validation <- ftq_validate_run(
    out$result_csv,
    out$rsc_csv,
    roi_summary_csv = out$roi_summary_csv,
    local_search_assignments = out$local_search_assignments_csv
  )
  expect_true(validation$passed)
  expect_true(any(grepl(
    "local_search_full_path_partition_reconciles",
    validation$checks$check,
    fixed = TRUE
  )))

  report <- ftq_report(
    result_csv = out$result_csv,
    rsc_csv = out$rsc_csv,
    manifest = list(
      prefix = "sub-test_desc-local-search",
      subject = "sub-test",
      description = "local-search-audit",
      focal_assignment_mode = "all_focal_voxels",
      streamline_voxelization_mode = "segment_interpolation",
      target_intersection_mode = "local_search",
      local_search_radius_mm = 0.75,
      endpoint_analysis = "none",
      output_dir = file.path(tmp, "outputs")
    ),
    output_dir = file.path(tmp, "report"),
    roi_summary = out$roi_summary_csv,
    local_search_assignments = out$local_search_assignments_csv,
    validation = validation,
    overwrite = TRUE
  )
  html <- paste(readLines(report$html, warn = FALSE), collapse = "\n")

  audit_link <- regmatches(html, regexpr(
    "href='[^']+'[^>]*>Download streamline-level local-search audit CSV", html
  ))
  audit_href <- sub("^href='([^']+)'.*$", "\\1", audit_link)
  expect_identical(
    normalizePath(file.path(dirname(report$html), utils::URLdecode(audit_href)),
                  winslash = "/", mustWork = TRUE),
    normalizePath(out$local_search_assignments_csv, winslash = "/", mustWork = TRUE)
  )

  expect_true(file.exists(report$tables$local_search_summary))
  expect_true(file.exists(report$tables$local_search_endpoint_status))
  expect_true(file.exists(report$figures$local_search_distances))
  expect_match(html, "Local-Search Connection Audit", fixed = TRUE)
  expect_match(html, "Proximity-only connection", fixed = TRUE)
  expect_match(
    html,
    "Download streamline-level local-search audit CSV",
    fixed = TRUE
  )
})

test_that("publication summaries keep unique streamlines separate from voxel associations", {
  result <- data.frame(
    fNode_Coords = c("1_1_1", "2_1_1"),
    tracks = c(4L, 3L),
    local = c(2L, 1L),
    leave_no_target = c(1L, 1L),
    T = c(3L, 2L),
    none = c(0L, 0L),
    stringsAsFactors = FALSE
  )
  rsc <- rsc_summary(result, percent = TRUE, warn_legacy = FALSE)
  roi <- data.frame(
    unique_focal_intersecting_streamlines = 5,
    unique_local_streamlines = 2,
    unique_leave_no_target_streamlines = 1,
    unique_leave_with_target_streamlines = 2,
    `unique_target_streamlines:T` = 3,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  fate <- FibrTraQ:::ftq_report_fate_summary(result, roi)
  expect_equal(fate$count, c(2, 1, 2))
  expect_equal(sum(fate$count), 5)
  expect_equal(sum(fate$percent_of_denominator), 100)
  expect_false(identical(sum(fate$count), sum(result$tracks)))

  connections <- FibrTraQ:::ftq_report_connection_summary(result, rsc, roi)
  expect_identical(connections$cNode, c("local", "T"))
  expect_equal(connections$unique_streamline_connections, c(2, 3))
  expect_equal(connections$connected_fNode_voxels, c(2, 2))
  expect_equal(connections$rsc_percent, c(100, 100))

  bad_roi <- roi
  bad_roi$unique_leave_with_target_streamlines <- 3
  expect_error(
    FibrTraQ:::ftq_report_fate_summary(result, bad_roi),
    "do not partition"
  )
})

test_that("report heatmap scaling is robust and display-only", {
  values <- c(rep(1, 1000), 10000)
  expect_equal(FibrTraQ:::ftq_positive_zlim(values), c(0, 10000))
  expect_equal(FibrTraQ:::ftq_positive_zlim(values, quantile = 0.99), c(0, 1))
  expect_error(
    FibrTraQ:::ftq_positive_zlim(values, quantile = 0),
    "quantile must"
  )
})

test_that("HTML reports escape metadata and encode artifact link paths", {
  tmp <- file.path(tempdir(), "fibrtraq-report-untrusted-html")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  payload <- "unsafe' onclick='globalThis.FIBRTRAQ_XSS=1"
  table_path <- file.path(tmp, paste0(payload, ".csv"))
  endpoint_path <- file.path(tmp, paste0("endpoint-", payload, ".csv"))
  writeLines("value\n1", table_path)
  writeLines("value\n1", endpoint_path)
  report_path <- file.path(tmp, "report.html")
  result <- data.frame(
    fNode_Coords = "1_1_1",
    tracks = 1L,
    local = 0L,
    leave_no_target = 0L,
    T = 1L,
    none = 0L,
    stringsAsFactors = FALSE
  )
  rsc <- rsc_summary(result, percent = TRUE, warn_legacy = FALSE)
  manifest <- list(
    subject = payload,
    session = "",
    description = payload,
    focal_assignment_mode = "all_focal_voxels",
    streamline_voxelization_mode = "sample_points",
    target_intersection_mode = "full_streamline",
    endpoint_analysis = "none",
    endpoint_context_labels = character(),
    reverse_streamlines = FALSE,
    discretization = "round",
    track_space = "voxel",
    output_dir = tmp,
    background_image = NULL
  )
  endpoint_data <- list(
    missing = FALSE,
    summary = data.frame(target = "T"),
    penetration = data.frame(
      target = "T",
      streamline_voxelization_mode = "segment_interpolation"
    ),
    penetration_summary = data.frame(target = "T"),
    context = data.frame(context = "focal"),
    source_paths = stats::setNames(list(endpoint_path), payload)
  )

  FibrTraQ:::ftq_write_html_report(
    path = report_path,
    result = result,
    rsc = rsc,
    manifest = manifest,
    validation = NULL,
    preflight = NULL,
    figs = list(),
    top_tables = stats::setNames(list(table_path), payload),
    endpoint_data = endpoint_data
  )

  html <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
  expect_false(grepl(" onclick='", html, fixed = TRUE))
  expect_true(grepl("&#39;", html, fixed = TRUE))
  expect_true(grepl(
    utils::URLencode(basename(table_path), reserved = TRUE),
    html,
    fixed = TRUE
  ))
  expect_true(grepl(
    utils::URLencode(basename(endpoint_path), reserved = TRUE),
    html,
    fixed = TRUE
  ))
})

test_that("manifest stream reader matches in-memory TRK reader", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("freesurferformats")

  tmp <- file.path(tempdir(), "stream-reader-equivalence")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  fnode <- file.path(tmp, "fnode.nii.gz")
  target <- file.path(tmp, "target.nii.gz")
  tractogram <- file.path(tmp, "tracts.trk")

  fnode_arr <- array(0, dim = c(4, 4, 4))
  fnode_arr[2, 2, 2] <- 1
  target_arr <- array(0, dim = c(4, 4, 4))
  target_arr[4, 2, 2] <- 1
  RNifti::writeNifti(fnode_arr, fnode)
  RNifti::writeNifti(target_arr, target)
  ftq_write_test_trk(
    tractogram,
    list(
      rbind(c(2, 2, 2), c(4, 2, 2)),
      rbind(c(2, 2, 2), c(2, 3, 2)),
      rbind(c(1, 1, 1), c(1, 1, 2))
    )
  )

  base_manifest <- list(
    dataset_root = tmp,
    subject = "sub-test",
    description = "stream-equivalence",
    fnode = basename(fnode),
    tractogram = basename(tractogram),
    targets = list(T = basename(target)),
    coordinate = list(
      track_space = "voxel",
      discretization = "round",
      axis_flip = c(FALSE, FALSE, FALSE)
    ),
    outputs = list(write_maps = FALSE, write_report = FALSE, overwrite = TRUE)
  )

  stream_manifest <- base_manifest
  stream_manifest$execution <- list(track_reader = "stream", track_chunk_size = 2L)
  stream_manifest$outputs$output_dir <- file.path(tmp, "stream-out")
  memory_manifest <- base_manifest
  memory_manifest$execution <- list(track_reader = "memory")
  memory_manifest$outputs$output_dir <- file.path(tmp, "memory-out")

  stream_run <- ftq_run(stream_manifest)
  memory_run <- ftq_run(memory_manifest)

  expect_equal(stream_run$result, memory_run$result)
  expect_true(attr(stream_run$result, "fibrtraq_local_roi_defaulted_to_focal"))
  expect_true(attr(memory_run$result, "fibrtraq_local_roi_defaulted_to_focal"))
  expect_true(stream_run$roi_summary$local_roi_defaulted_to_focal)
  expect_equal(stream_run$result$tracks, 2L)
  expect_equal(stream_run$result$T, 1L)
  expect_equal(stream_run$result$leave_no_target, 1L)
})

test_that("tractogram context builds discretized streamline occupancy density", {
  skip_if_not_installed("RNifti")

  tmp <- file.path(tempdir(), "tractogram-context-density")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  template <- file.path(tmp, "template.nii.gz")
  tractogram <- file.path(tmp, "tracts.trk")
  RNifti::writeNifti(array(0, dim = c(5, 5, 5)), template)
  ftq_write_test_trk(
    tractogram,
    list(
      rbind(c(2, 2, 2), c(2, 2, 2), c(3, 2, 2), c(6, 2, 2)),
      rbind(c(2, 2, 2), c(4, 4, 4))
    ),
    dims = c(5L, 5L, 5L)
  )

  context <- ftq_tractogram_context(
    tractogram = tractogram,
    template_image = template,
    output_dir = tmp,
    prefix = "sub-test_desc-context",
    track_space = "voxel",
    axis_flip = c(FALSE, FALSE, FALSE),
    track_chunk_size = 1L,
    overwrite = TRUE
  )
  density <- as.array(RNifti::readNifti(context$density_image, internal = FALSE))
  expect_equal(density[2, 2, 2], 2)
  expect_equal(density[3, 2, 2], 1)
  expect_equal(density[4, 4, 4], 1)
  expect_equal(sum(density > 0), 3)
  expect_equal(context$summary$total_streamline_voxel_occupancy_events, 4)
  expect_equal(context$summary$out_of_bounds_sample_points, 1)
  expect_match(context$summary$note, "not anatomical brain volume")
})

test_that("tractogram context applies axis flip on the template grid", {
  skip_if_not_installed("RNifti")

  tmp <- file.path(tempdir(), "tractogram-context-axis-flip")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  template <- file.path(tmp, "template.nii.gz")
  tractogram <- file.path(tmp, "tracts.trk")
  RNifti::writeNifti(array(0, dim = c(5, 5, 5)), template)
  ftq_write_test_trk(tractogram, list(rbind(c(1, 2, 2), c(1, 2, 2))), dims = c(5L, 5L, 5L))

  context <- ftq_tractogram_context(
    tractogram = tractogram,
    template_image = template,
    output_dir = tmp,
    prefix = "sub-test_desc-flip",
    track_space = "voxel",
    axis_flip = c(TRUE, FALSE, FALSE),
    overwrite = TRUE
  )
  density <- as.array(RNifti::readNifti(context$density_image, internal = FALSE))
  expect_equal(density[5, 2, 2], 1)
  expect_equal(density[1, 2, 2], 0)
})

test_that("tractogram context uses TrackVis voxel-mm header scaling", {
  skip_if_not_installed("RNifti")

  tmp <- file.path(tempdir(), "tractogram-context-voxelmm")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  template <- file.path(tmp, "template.nii.gz")
  tractogram <- file.path(tmp, "tracts.trk")
  RNifti::writeNifti(array(0, dim = c(5, 5, 5)), template)
  ftq_write_test_trk(
    tractogram,
    list(rbind(c(0.05, 0.1, 0.25), c(0.35, 0.1, 0.25))),
    dims = c(5L, 5L, 5L),
    voxel_size = c(0.1, 0.2, 0.5)
  )

  context <- ftq_tractogram_context(
    tractogram = tractogram,
    template_image = template,
    output_dir = tmp,
    prefix = "sub-test_desc-voxelmm",
    track_space = "auto",
    axis_flip = c(FALSE, FALSE, FALSE),
    streamline_voxelization_mode = "segment_interpolation",
    overwrite = TRUE
  )
  density <- as.array(RNifti::readNifti(context$density_image, internal = FALSE))
  expect_equal(density[1, 1, 1], 1)
  expect_equal(density[2, 1, 1], 1)
  expect_equal(density[3, 1, 1], 1)
  expect_equal(density[4, 1, 1], 1)
  expect_identical(context$metadata$selected_track_space, "voxelmm")
  expect_equal(context$metadata$track_voxel_size_mm, c(0.1, 0.2, 0.5), tolerance = 1e-6)
})

test_that("ftq_run automatically writes tractogram context when anatomy is absent", {
  skip_if_not_installed("RNifti")

  tmp <- file.path(tempdir(), "workflow-context-fallback")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  fnode <- file.path(tmp, "fnode.nii.gz")
  target <- file.path(tmp, "target.nii.gz")
  tractogram <- file.path(tmp, "tracts.trk")

  fnode_arr <- array(0, dim = c(5, 5, 5))
  fnode_arr[2, 2, 2] <- 1
  target_arr <- array(0, dim = c(5, 5, 5))
  target_arr[4, 2, 2] <- 1
  RNifti::writeNifti(fnode_arr, fnode)
  RNifti::writeNifti(target_arr, target)
  ftq_write_test_trk(
    tractogram,
    list(
      rbind(c(2, 2, 2), c(4, 2, 2)),
      rbind(c(2, 2, 2), c(3, 2, 2))
    ),
    dims = c(5L, 5L, 5L)
  )

  run <- ftq_run(list(
    dataset_root = tmp,
    subject = "sub-test",
    description = "context-fallback",
    fnode = basename(fnode),
    tractogram = basename(tractogram),
    targets = list(T = basename(target)),
    coordinate = list(track_space = "voxel", discretization = "round", axis_flip = c(FALSE, FALSE, FALSE)),
    execution = list(track_reader = "stream", track_chunk_size = 1L),
    outputs = list(output_dir = file.path(tmp, "out"), write_maps = TRUE, write_report = TRUE, overwrite = TRUE)
  ))

  expect_true(file.exists(run$outputs$context$density_image))
  expect_null(run$outputs$context$summary_csv)
  expect_identical(basename(dirname(run$outputs$context$density_image)), "maps")
  expect_false(dir.exists(file.path(run$manifest$output_dir, "context")))
  html <- paste(readLines(run$report$html, warn = FALSE), collapse = "\n")
  expect_match(html, "Tractogram-derived density context", fixed = TRUE)
})

test_that("ftq_run prefers supplied anatomical background over automatic context", {
  skip_if_not_installed("RNifti")

  tmp <- file.path(tempdir(), "workflow-context-prefers-anatomy")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  fnode <- file.path(tmp, "fnode.nii.gz")
  target <- file.path(tmp, "target.nii.gz")
  background <- file.path(tmp, "background.nii.gz")
  tractogram <- file.path(tmp, "tracts.trk")
  fnode_arr <- array(0, dim = c(4, 4, 4))
  fnode_arr[2, 2, 2] <- 1
  target_arr <- array(0, dim = c(4, 4, 4))
  target_arr[4, 2, 2] <- 1
  RNifti::writeNifti(fnode_arr, fnode)
  RNifti::writeNifti(target_arr, target)
  RNifti::writeNifti(array(seq_len(64), dim = c(4, 4, 4)), background)
  ftq_write_test_trk(tractogram, list(rbind(c(2, 2, 2), c(4, 2, 2))), dims = c(4L, 4L, 4L))

  run <- ftq_run(list(
    dataset_root = tmp,
    subject = "sub-test",
    description = "anatomy-preferred",
    fnode = basename(fnode),
    tractogram = basename(tractogram),
    targets = list(T = basename(target)),
    report = list(background_image = basename(background)),
    coordinate = list(track_space = "voxel", discretization = "round", axis_flip = c(FALSE, FALSE, FALSE)),
    outputs = list(output_dir = file.path(tmp, "out"), write_maps = TRUE, write_report = TRUE, overwrite = TRUE)
  ))

  expect_null(run$outputs$context)
  html <- paste(readLines(run$report$html, warn = FALSE), collapse = "\n")
  expect_false(grepl("Tractogram-derived density context", html, fixed = TRUE))
})

test_that("target report derives primary targets without research files", {
  skip_if_not_installed("RNifti")

  tmp <- file.path(tempdir(), "exploratory-summary")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  template <- file.path(tmp, "template.nii.gz")
  background <- file.path(tmp, "background.nii.gz")
  RNifti::writeNifti(array(0, dim = c(4, 4, 4)), template)
  RNifti::writeNifti(array(seq_len(64), dim = c(4, 4, 4)), background)

  result <- data.frame(
    fNode_Coords = c("1_1_1", "2_2_2", "3_3_3"),
    tracks = c(10L, 5L, 0L),
    local = c(1L, 0L, 0L),
    leave_no_target = c(2L, 1L, 0L),
    A = c(6L, 1L, 0L),
    B = c(1L, 3L, 0L),
    none = c(0L, 0L, 1L),
    stringsAsFactors = FALSE
  )

  report <- ftq_exploratory_report(
    result_csv = result,
    output_dir = file.path(tmp, "out"),
    template_image = template,
    background_image = background,
    prefix = "sub-test_desc-explore",
    overwrite = TRUE
  )

  expect_equal(report$summary$data$primary_target, c("A", "B", "none"))
  expect_equal(report$summary$data$target_diversity, c(2, 2, 0))
  expect_equal(report$summary$target_totals$target[1], "A")
  expect_null(report$summary$outputs$summary_csv)
  expect_length(report$summary$outputs$maps, 0L)
  expect_null(report$summary$outputs$provenance_csv)
  expect_true(file.exists(report$summary$outputs$lut_csv))
  expect_equal(
    report$summary$provenance$target_intersection_mode,
    "full_streamline"
  )
  expect_true(file.exists(report$html))
  expect_true(all(file.exists(unlist(report$figures))))
  report_html <- paste(readLines(report$html, warn = FALSE), collapse = "\n")
  expect_match(
    report_html,
    "Target assignment mode:</b> full_streamline",
    fixed = TRUE
  )

  local_result <- result
  attr(local_result, "fibrtraq_target_intersection_mode") <- "local_search"
  attr(local_result, "fibrtraq_local_search_radius_mm") <- 3.25
  local_summary <- ftq_exploratory_summary(
    local_result,
    write_maps = FALSE
  )
  expect_equal(
    local_summary$provenance$target_intersection_mode,
    "local_search"
  )
  expect_equal(local_summary$provenance$local_search_radius_mm, 3.25)

  context <- file.path(tmp, "context.nii.gz")
  context_arr <- array(0, dim = c(4, 4, 4))
  context_arr[1:3, 1:3, 1:3] <- 1
  RNifti::writeNifti(context_arr, context)
  context_report <- ftq_exploratory_report(
    result_csv = result,
    output_dir = file.path(tmp, "out-context"),
    template_image = template,
    prefix = "sub-test_desc-explore-context",
    overwrite = TRUE,
    tractogram_context = context
  )
  context_html <- paste(readLines(context_report$html, warn = FALSE), collapse = "\n")
  expect_match(context_html, "Tractogram-derived density context", fixed = TRUE)

  unsafe_prefix <- "unsafe' onerror='globalThis.FIBRTRAQXSS=1"
  unsafe_report <- ftq_exploratory_report(
    result_csv = result,
    output_dir = file.path(tmp, "out-untrusted-prefix"),
    template_image = template,
    prefix = unsafe_prefix,
    overwrite = TRUE
  )
  unsafe_html <- paste(
    readLines(unsafe_report$html, warn = FALSE),
    collapse = "\n"
  )
  expect_false(grepl(" onerror='", unsafe_html, fixed = TRUE))
  expect_true(grepl("%27", unsafe_html, fixed = TRUE))
  # Prefixes now appear only in URL-encoded links, not visible filename labels.
})

test_that("preflight catches blocking dataset issues and non-blocking recommendations", {
  skip_if_not_installed("RNifti")

  tmp <- tempdir()
  fnode <- file.path(tmp, "fnode.nii.gz")
  target <- file.path(tmp, "target.nii.gz")
  bad_target <- file.path(tmp, "bad-target.nii.gz")
  tractogram <- file.path(tmp, "tracts.trk")
  RNifti::writeNifti(array(c(1, rep(0, 26)), dim = c(3, 3, 3)), fnode)
  RNifti::writeNifti(array(c(0, 1, rep(0, 25)), dim = c(3, 3, 3)), target)
  RNifti::writeNifti(array(0, dim = c(4, 3, 3)), bad_target)
  file.create(tractogram)

  manifest <- list(
    dataset_root = tmp,
    subject = "sub-001",
    description = "preflight",
    fnode = "fnode.nii.gz",
    tractogram = "tracts.trk",
    targets = list(T = "target.nii.gz"),
    outputs = list(output_dir = file.path(tmp, "preflight-out"), overwrite = FALSE)
  )
  preflight <- ftq_preflight(manifest)
  expect_true(preflight$passed)
  expect_true("spatial_context_fallback_available" %in% preflight$checks$check)
  expect_true("focal_local_overlap" %in% preflight$checks$check)
  expect_true(preflight$checks$passed[preflight$checks$check == "spatial_context_fallback_available"])
  expect_equal(preflight$checks$severity[preflight$checks$check == "spatial_context_fallback_available"], "recommendation")
  expect_false(
    preflight$checks$passed[
      preflight$checks$check == "ref_image_recommended_for_auto_mapping"
    ]
  )

  auto_with_reference <- manifest
  auto_with_reference$ref_image <- "fnode.nii.gz"
  reference_preflight <- ftq_preflight(auto_with_reference)
  expect_true(reference_preflight$passed)
  expect_true(
    reference_preflight$checks$passed[
      reference_preflight$checks$check ==
        "ref_image_recommended_for_auto_mapping"
    ]
  )

  world_without_reference <- manifest
  world_without_reference$coordinate <- list(track_space = "world")
  world_preflight <- ftq_preflight(world_without_reference)
  expect_false(world_preflight$passed)
  expect_false(
    world_preflight$checks$passed[
      world_preflight$checks$check ==
        "ref_image_required_for_world_mapping"
    ]
  )

  no_map_manifest <- manifest
  no_map_manifest$outputs$write_maps <- FALSE
  no_map_preflight <- ftq_preflight(no_map_manifest)
  expect_true(no_map_preflight$passed)
  expect_true(no_map_preflight$checks$passed[no_map_preflight$checks$check == "map_outputs_available"])

  bad_manifest <- manifest
  bad_manifest$targets <- list(T = "bad-target.nii.gz")
  bad_preflight <- ftq_preflight(bad_manifest)
  expect_false(bad_preflight$passed)
  expect_false(bad_preflight$checks$passed[bad_preflight$checks$check == "target_dims_match_fnode:T"])

  planned <- file.path(tmp, "preflight-out", "sub-001_desc-preflight_fibrtraq_per-voxel-hits.csv")
  dir.create(dirname(planned), recursive = TRUE, showWarnings = FALSE)
  file.create(planned)
  collision_preflight <- ftq_preflight(manifest)
  expect_false(collision_preflight$passed)
  expect_false(collision_preflight$checks$passed[collision_preflight$checks$check == "core_outputs_available"])
})

test_that("manifest parser resolves generic dataset paths", {
  tmp <- tempdir()
  manifest <- list(
    dataset_root = tmp,
    subject = "sub-001",
    session = "ses-01",
    description = "demo",
    fnode = "rois/fnode.nii.gz",
    tractogram = "trk/test.trk",
    targets = list(A = "rois/a.nii.gz"),
    report = list(background_image = "anat/t1w.nii.gz"),
    coordinate = list(track_space = "auto", discretization = "any", axis_flip = "auto"),
    outputs = list(
      derivative_name = "fibrtraq-test",
      write_maps = TRUE,
      expected = list(
        result_csv = "expected/per-voxel.csv",
        rsc_csv = "expected/rsc.csv",
        maps = list(tracks = "expected/tracks.nii.gz")
      )
    )
  )

  parsed <- ftq_read_manifest(manifest)
  expect_equal(parsed$prefix, "sub-001_ses-01_desc-demo")
  expect_equal(parsed$target_labels, "A")
  expect_null(parsed$axis_flip)
  expect_true(grepl("fibrtraq-test$", parsed$output_dir))
  expect_true(grepl("anat/t1w.nii.gz$", ftq_test_path(parsed$background_image)))
  expect_true(grepl("expected/per-voxel.csv$", ftq_test_path(parsed$expected$result_csv)))
  expect_true(grepl("expected/tracks.nii.gz$", ftq_test_path(parsed$expected$map_paths$tracks)))
})

test_that("manifest parser preserves Windows absolute paths", {
  manifest <- list(
    dataset_root = "D:/dataset-root",
    subject = "sub-001",
    description = "windows-absolute",
    fnode = "C:/data/fnode.nii.gz",
    tractogram = "C:/data/tracts.trk",
    targets = list(A = "C:/data/target.nii.gz"),
    report = list(background_image = "C:/data/t1w.nii.gz"),
    outputs = list(output_dir = "C:/data/out")
  )

  parsed <- ftq_read_manifest(manifest)
  expect_false(grepl("dataset-root.*C:/data", ftq_test_path(parsed$fnode), fixed = FALSE))
  expect_true(grepl("^C:/data/fnode.nii.gz$", ftq_test_path(parsed$fnode)))
  expect_true(grepl("^C:/data/tracts.trk$", ftq_test_path(parsed$tractogram)))
  expect_true(grepl("^C:/data/target.nii.gz$", ftq_test_path(parsed$target_paths[1])))
  expect_true(grepl("^C:/data/t1w.nii.gz$", ftq_test_path(parsed$background_image)))
  expect_true(grepl("^C:/data/out$", ftq_test_path(parsed$output_dir)))
})

test_that("installed validation manifest templates parse", {
  manifest_dir <- system.file("validation", "manifests", package = "FibrTraQ")
  skip_if(!nzchar(manifest_dir), "validation templates not installed")

  manifests <- list.files(manifest_dir, pattern = "\\.yml$", full.names = TRUE)
  expect_gt(length(manifests), 0)
  parsed <- lapply(manifests, ftq_read_manifest)

  expect_true(all(vapply(parsed, function(x) length(x$target_labels) > 0L, logical(1))))
  expect_true(all(vapply(parsed, function(x) identical(x$track_space, "auto"), logical(1))))
  expect_true(all(vapply(parsed, function(x) isTRUE(x$write_report), logical(1))))
})

test_that("validation suite summarizes manifest preflight rows", {
  skip_if_not_installed("RNifti")

  tmp <- tempdir()
  fnode <- file.path(tmp, "suite-fnode.nii.gz")
  target <- file.path(tmp, "suite-target.nii.gz")
  tractogram <- file.path(tmp, "suite-tracts.trk")
  RNifti::writeNifti(array(c(1, rep(0, 26)), dim = c(3, 3, 3)), fnode)
  RNifti::writeNifti(array(c(0, 1, rep(0, 25)), dim = c(3, 3, 3)), target)
  file.create(tractogram)

  manifest <- list(
    dataset_root = tmp,
    subject = "sub-suite",
    description = "validation-suite",
    fnode = basename(fnode),
    tractogram = basename(tractogram),
    targets = list(T = basename(target)),
    outputs = list(output_dir = file.path(tmp, "suite-out"), overwrite = FALSE)
  )
  manifest_path <- file.path(tmp, "suite-manifest.yml")
  yaml::write_yaml(manifest, manifest_path)
  summary_csv <- file.path(tmp, "suite-summary.csv")

  summary <- ftq_run_validation_suite(
    manifests = manifest_path,
    output_csv = summary_csv,
    preflight_only = TRUE
  )

  expect_true(file.exists(summary_csv))
  expect_equal(nrow(summary), 1L)
  expect_equal(summary$run_status, "preflight_only")
  expect_true(summary$preflight_passed)
  expect_equal(summary$preflight_recommendations, 0L)
  expect_true(is.na(summary$validation_passed))
})

test_that("ISMRM 2015 manifest generator maps configured endpoints", {
  tmp <- file.path(tempdir(), "ismrm-generator")
  scoring <- file.path(tmp, "scoring")
  manifest_dir <- file.path(tmp, "manifests")
  dir.create(file.path(scoring, "bundles"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(scoring, "bundles", "sub_bundles"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(scoring, "ROI", "endpoints"), recursive = TRUE, showWarnings = FALSE)

  file.create(file.path(scoring, "bundles", "CP.trk"))
  file.create(file.path(scoring, "bundles", "sub_bundles", "CC_u_shaped.trk"))
  file.create(file.path(scoring, "ROI", "endpoints", "CP_head.nii.gz"))
  file.create(file.path(scoring, "ROI", "endpoints", "CP_tail.nii.gz"))
  file.create(file.path(scoring, "ROI", "endpoints", "CC_left.nii.gz"))
  file.create(file.path(scoring, "ROI", "endpoints", "CC_right.nii.gz"))
  file.create(file.path(scoring, "t1.nii.gz"))
  jsonlite::write_json(
    list(
      CP = list(
        all_mask = "ROI/all_masks/CP.nii.gz",
        head = "ROI/endpoints/CP_head.nii.gz",
        tail = "ROI/endpoints/CP_tail.nii.gz"
      ),
      CC_u_shaped = list(
        all_mask = "ROI/all_masks/CC_u_shaped.nii.gz",
        head = "ROI/endpoints/CC_left.nii.gz",
        tail = "ROI/endpoints/CC_right.nii.gz"
      )
    ),
    file.path(scoring, "config_file_segmentation.json"),
    auto_unbox = TRUE,
    pretty = TRUE
  )

  generated <- ftq_make_ismrm2015_manifests(
    scoring_dir = scoring,
    output_dir = manifest_dir,
    run_output_root = file.path(tmp, "outputs"),
    force = TRUE
  )

  expect_equal(generated$bundle, c("CP", "CC_u_shaped"))
  expect_true(all(file.exists(generated$manifest_path)))

  parsed <- ftq_read_manifest(generated$manifest_path[1])
  expect_equal(parsed$description, "ismrm2015-CP-groundtruth")
  expect_true(grepl("bundles/CP.trk$", ftq_test_path(parsed$tractogram)))
  expect_true(grepl("ROI/endpoints/CP_head.nii.gz$", ftq_test_path(parsed$fnode)))
  expect_equal(parsed$target_labels, "tail")
  expect_true(grepl("ROI/endpoints/CP_tail.nii.gz$", ftq_test_path(parsed$target_paths[1])))
  expect_true(grepl("t1.nii.gz$", ftq_test_path(parsed$background_image)))
})

test_that("ISMRM 2015 control manifest generators build specificity tests", {
  skip_if_not_installed("RNifti")

  tmp <- file.path(tempdir(), "ismrm-control-generator")
  scoring <- file.path(tmp, "scoring")
  dir.create(file.path(scoring, "bundles"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(scoring, "ROI", "endpoints"), recursive = TRUE, showWarnings = FALSE)

  file.create(file.path(scoring, "bundles", "CP.trk"))
  file.create(file.path(scoring, "bundles", "MCP.trk"))
  RNifti::writeNifti(array(c(1, rep(0, 26)), dim = c(3, 3, 3)),
                     file.path(scoring, "ROI", "endpoints", "CP_head.nii.gz"))
  RNifti::writeNifti(array(c(0, 1, rep(0, 25)), dim = c(3, 3, 3)),
                     file.path(scoring, "ROI", "endpoints", "CP_tail.nii.gz"))
  RNifti::writeNifti(array(c(0, 0, 1, rep(0, 24)), dim = c(3, 3, 3)),
                     file.path(scoring, "ROI", "endpoints", "MCP_head.nii.gz"))
  RNifti::writeNifti(array(c(0, 0, 0, 1, rep(0, 23)), dim = c(3, 3, 3)),
                     file.path(scoring, "ROI", "endpoints", "MCP_tail.nii.gz"))
  RNifti::writeNifti(array(seq_len(27), dim = c(3, 3, 3)),
                     file.path(scoring, "t1.nii.gz"))
  jsonlite::write_json(
    list(
      CP = list(head = "ROI/endpoints/CP_head.nii.gz", tail = "ROI/endpoints/CP_tail.nii.gz"),
      MCP = list(head = "ROI/endpoints/MCP_head.nii.gz", tail = "ROI/endpoints/MCP_tail.nii.gz")
    ),
    file.path(scoring, "config_file_segmentation.json"),
    auto_unbox = TRUE,
    pretty = TRUE
  )

  wrong <- ftq_make_ismrm2015_wrong_target_manifests(
    scoring_dir = scoring,
    output_dir = file.path(tmp, "wrong"),
    bundles = "CP",
    target_bundles = c("CP", "MCP"),
    run_output_root = file.path(tmp, "wrong-out"),
    force = TRUE
  )
  expect_equal(nrow(wrong), 1L)
  expect_equal(wrong$target_bundle, "MCP")
  wrong_parsed <- ftq_read_manifest(wrong$manifest_path)
  expect_equal(wrong_parsed$target_labels, "wrong_MCP")
  expect_false(wrong_parsed$write_maps)

  wrong_matrix <- ftq_make_ismrm2015_wrong_target_manifests(
    scoring_dir = scoring,
    output_dir = file.path(tmp, "wrong-matrix"),
    bundles = "CP",
    target_bundles = c("CP", "MCP"),
    run_output_root = file.path(tmp, "wrong-matrix-out"),
    combine_targets = TRUE,
    force = TRUE
  )
  expect_equal(nrow(wrong_matrix), 1L)
  expect_equal(wrong_matrix$target_count, 1L)
  wrong_matrix_parsed <- ftq_read_manifest(wrong_matrix$manifest_path)
  expect_equal(wrong_matrix_parsed$target_labels, "wrong_MCP")
  expect_true(grepl("all-wrong-targets$", wrong_matrix_parsed$output_dir))

  shifted <- ftq_make_ismrm2015_shifted_target_manifests(
    scoring_dir = scoring,
    output_dir = file.path(tmp, "shifted"),
    shifted_roi_dir = file.path(tmp, "shifted-rois"),
    shift = c(1L, 0L, 0L),
    bundles = "CP",
    run_output_root = file.path(tmp, "shifted-out"),
    force = TRUE
  )
  expect_equal(nrow(shifted), 1L)
  expect_true(file.exists(shifted$target))
  shifted_arr <- RNifti::readNifti(shifted$target, internal = FALSE)
  expect_equal(sum(as.array(shifted_arr) != 0), 1)
  shifted_parsed <- ftq_read_manifest(shifted$manifest_path)
  expect_equal(shifted_parsed$target_labels, "shifted_tail")
  expect_true(grepl("shifted-rois", shifted_parsed$target_paths[1]))

  submission <- file.path(tmp, "submission.trk")
  file.create(submission)
  sweep <- ftq_make_ismrm2015_submission_sweep_manifests(
    scoring_dir = scoring,
    tractograms = submission,
    output_dir = file.path(tmp, "sweep"),
    bundles = c("CP", "MCP"),
    run_output_root = file.path(tmp, "sweep-out"),
    force = TRUE
  )
  expect_equal(nrow(sweep), 2L)
  expect_equal(sweep$bundle, c("CP", "MCP"))
  sweep_parsed <- ftq_read_manifest(sweep$manifest_path[1])
  expect_true(grepl("submission.trk$", sweep_parsed$tractogram))
})
