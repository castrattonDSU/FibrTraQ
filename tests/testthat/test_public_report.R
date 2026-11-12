ftq_public_report_result <- function() {
  coords <- expand.grid(x = 1:4, y = 1:3, z = 1:5)
  data.frame(
    fNode_Coords = apply(coords, 1L, paste, collapse = "_"),
    tracks = rep(c(0, 4, 8), 20L),
    local = rep(c(0, 1, 8), 20L),
    leave_no_target = rep(c(0, 1, 0), 20L),
    x = rep(c(0, 2, 0), 20L),
    y = rep(c(0, 1, 0), 20L),
    z = rep(c(0, 0, 2), 20L),
    voxel_x = rep(c(0, 1, 0), 20L),
    voxel_x_2 = rep(c(0, 0, 1), 20L),
    none = rep(c(1, 0, 0), 20L),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

test_that("complete voxel exports preserve coordinate-named targets on readback", {
  root <- tempfile("public-voxel-view-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  result <- ftq_public_report_result()
  original <- result
  paths <- FibrTraQ:::ftq_report_write_voxel_exports(result, root, "fixture")
  view <- utils::read.csv(paths$tables$all_voxels, check.names = FALSE)
  roles <- utils::read.csv(paths$tables$all_voxel_columns, check.names = FALSE)

  expect_equal(nrow(view), 60L)
  expect_equal(anyDuplicated(names(view)), 0L)
  expect_identical(roles$column, names(view))
  expect_identical(roles$column[roles$role == "coordinate"],
                   c("voxel_x_3", "voxel_y", "voxel_z"))
  expect_setequal(roles$column[roles$role == "target"],
                  c("x", "y", "z", "voxel_x", "voxel_x_2"))
  expect_equal(view[names(result)], result)
  expect_equal(unname(as.matrix(view[roles$column[roles$role == "coordinate"]])),
               unname(as.matrix(expand.grid(x = 1:4, y = 1:3, z = 1:5))))
  expect_equal(FibrTraQ:::ftq_read_table(paths$tables$all_voxels), result)
  expect_identical(result, original)

  canonical <- file.path(root, "fixture_fibrtraq_per-voxel-hits.csv")
  utils::write.csv(result, canonical, row.names = FALSE)
  expect_equal(FibrTraQ:::ftq_read_table(canonical), result)
  expect_error(FibrTraQ:::ftq_report_write_voxel_exports(result, root, "fixture"),
               "already exist")

  view$voxel_x_3[1L] <- 99L
  utils::write.csv(view, paths$tables$all_voxels, row.names = FALSE)
  expect_error(FibrTraQ:::ftq_read_table(paths$tables$all_voxels),
               "coordinates do not reconstruct")
  unlink(paths$tables$all_voxel_columns)
  expect_error(FibrTraQ:::ftq_read_table(paths$tables$all_voxels),
               "companion all-voxels-columns")
})

test_that("ranked views include all voxels with deterministic ties and no repeated fields", {
  root <- tempfile("public-ranked-voxels-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  result <- ftq_public_report_result()
  paths <- FibrTraQ:::ftq_write_top_tables(result, root, "fixture")
  expect_setequal(names(paths), c("leave_no_target", "x", "y", "z", "voxel_x", "voxel_x_2"))
  no_target <- utils::read.csv(paths$leave_no_target, check.names = FALSE)
  expect_true(all(grepl("_ranked_.*_voxels[.]csv$", unlist(paths))))
  expect_false(any(grepl("_top50_", list.files(root))))
  coords <- FibrTraQ:::ftq_coord_matrix(result$fNode_Coords)
  for (metric in names(paths)) {
    ranked <- utils::read.csv(paths[[metric]], check.names = FALSE)
    ord <- order(-result[[metric]], coords[, 1L], coords[, 2L], coords[, 3L])
    expect_equal(nrow(ranked), nrow(result))
    expect_equal(ranked[[metric]], result[[metric]][ord])
    expect_identical(apply(ranked[c("voxel_x_3", "voxel_y", "voxel_z")],
                           1L, paste, collapse = "_"), result$fNode_Coords[ord])
    expect_equal(sum(ranked$tracks == 0), sum(result$tracks == 0))
    expect_equal(anyDuplicated(names(ranked)), 0L)
  }
  expect_equal(sum(names(no_target) == "leave_no_target"), 1L)
  expect_equal(anyDuplicated(names(no_target)), 0L)
  expect_true(all(diff(no_target$leave_no_target) <= 0))
  target_x <- utils::read.csv(paths$x, check.names = FALSE)
  expect_true(all(c("voxel_x_3", "x") %in% names(target_x)))
  expect_true(all(diff(target_x$x) <= 0))
  full <- FibrTraQ:::ftq_report_write_voxel_exports(result, root, "fixture")
  expect_equal(nrow(FibrTraQ:::ftq_read_table(full$tables$all_voxels)), nrow(result))
})

test_that("explicit limited ranked views remain available for research callers", {
  root <- tempfile("research-ranked-voxels-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  result <- ftq_public_report_result()
  paths <- FibrTraQ:::ftq_write_top_tables(result, root, "fixture", n = 50L)
  expect_true(all(grepl("_top50_", unlist(paths))))
  expect_true(all(vapply(paths, function(path) nrow(utils::read.csv(path)) == 50L, logical(1))))
  for (invalid in list(-1L, NA_real_, Inf, 1.5, c(1L, 2L), "50")) {
    expect_error(FibrTraQ:::ftq_write_top_tables(result, root, "invalid", n = invalid),
                 "single nonnegative integer")
  }
})

test_that("main reports use public RSC labels, complete exports and portable figure paths", {
  root <- tempfile("public report # 100%20 ")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  result <- ftq_public_report_result()
  original <- result
  rsc <- rsc_summary(result, percent = TRUE, warn_legacy = FALSE)
  expected_rsc <- FibrTraQ:::ftq_report_connection_summary(result, rsc)
  report <- ftq_report(result, rsc, output_dir = root,
                       manifest = list(prefix = "fixture", subject = "synthetic"))
  html <- paste(readLines(report$html, warn = FALSE), collapse = "\n")
  expect_match(html, "<h2>RSC summary</h2>", fixed = TRUE)
  expect_match(html, "Total number of fNode voxels: 60.</p>", fixed = TRUE)
  expect_match(html, "<th>Destination</th>", fixed = TRUE)
  expect_match(html, "RSC summary (CSV)", fixed = TRUE)
  expect_match(html, "All voxels, ranked by x (CSV)", fixed = TRUE)
  expect_false(grepl("Connection summary|RSC denominator|Top-50", html))
  expect_false(dir.exists(file.path(root, "figures")))
  expect_true(all(basename(dirname(unlist(report$figures))) == "report-figures"))
  expect_true(all(file.exists(unlist(report$figures))))
  expect_match(html, "src='report-figures/fixture_fig_rsc-summary.png'", fixed = TRUE)
  observed_rsc <- utils::read.csv(report$tables$connection_summary, check.names = FALSE)
  expect_identical(names(observed_rsc), names(expected_rsc))
  expect_identical(observed_rsc$cNode, expected_rsc$cNode)
  for (column in setdiff(names(expected_rsc), "cNode")) {
    expect_equal(as.numeric(observed_rsc[[column]]), as.numeric(expected_rsc[[column]]))
  }
  expect_identical(result, original)
})

test_that("context settings and summary survive without standalone sidecars", {
  root <- tempfile("public-context-metadata-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  result <- ftq_public_report_result()
  summary <- data.frame(streamline_count = 8L, occupied_voxels = 12L,
                         note = "Not anatomical brain volume.")
  context <- list(
    density_image = file.path(root, "maps", "fixture_fibrtraq_tractogram-density.nii.gz"),
    summary_csv = NULL, metadata_json = NULL, summary = summary,
    metadata = list(density_unit = "streamline occupancy per voxel",
                    track_chunk_size = 50000L, axis_flip = c(FALSE, TRUE, FALSE),
                    streamline_voxelization_mode = "segment_interpolation",
                    summary = summary)
  )
  original <- context
  outputs <- ftq_write_outputs(result, rsc_summary(result, percent = TRUE, warn_legacy = FALSE),
                               template_image = NULL, output_dir = root, prefix = "fixture",
                               context = context, write_maps = FALSE)
  metadata <- jsonlite::read_json(outputs$metadata_json, simplifyVector = TRUE)
  saved <- metadata$files$context
  expect_equal(saved$summary, summary)
  expect_equal(saved$metadata, context$metadata[setdiff(names(context$metadata), "summary")])
  expect_null(saved$summary_csv)
  expect_null(saved$metadata_json)
  expect_false("summary" %in% names(saved$metadata))
  methods <- FibrTraQ:::ftq_report_metadata_html(outputs$metadata_json, root)
  expect_match(methods, "tractogram context / summary", fixed = TRUE)
  expect_match(methods, "streamline occupancy per voxel", fixed = TRUE)
  expect_match(methods, "segment_interpolation", fixed = TRUE)
  expect_identical(context, original)
  expect_identical(FibrTraQ:::ftq_compact_context(context$density_image)$density_image,
                   context$density_image)
})

test_that("derived fate columns do not replace target names or their rankings", {
  root <- tempfile("public-derived-name-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  result <- data.frame(
    fNode_Coords = c("1_1_1", "2_1_1", "3_1_1"), tracks = c(10, 10, 10),
    local = c(1, 2, 3), leave_no_target = c(2, 2, 2), none = 0,
    leave_with_target = c(1, 9, 2), derived_leave_with_target = c(8, 1, 3)
  )
  paths <- FibrTraQ:::ftq_write_top_tables(result, root, "fixture")
  target <- utils::read.csv(paths$leave_with_target, check.names = FALSE)
  fallback_target <- utils::read.csv(paths$derived_leave_with_target, check.names = FALSE)
  expect_equal(target$leave_with_target, c(9, 2, 1))
  expect_equal(target$x, c(2, 3, 1))
  expect_equal(target$derived_derived_leave_with_target, c(6, 5, 7))
  expect_equal(fallback_target$derived_leave_with_target, c(8, 3, 1))
  expect_equal(fallback_target$derived_derived_leave_with_target, c(7, 5, 6))
  expect_equal(anyDuplicated(names(target)), 0L)
})

test_that("local fractions mark zero assignments missing and reject invalid ratios", {
  result <- data.frame(fNode_Coords = c("1_1_1", "2_2_2", "3_3_3", "4_4_4"),
                       tracks = c(0, 4, 8, 2), local = c(0, 1, 8, 0))
  fraction <- FibrTraQ:::ftq_report_local_fraction(result)
  expect_equal(fraction$local_fraction, c(NA_real_, 0.25, 1, 0))
  expect_identical(fraction$has_assigned_streamlines, c(FALSE, TRUE, TRUE, TRUE))
  expect_identical(fraction$fNode_Coords, result$fNode_Coords)
  for (invalid in c(-1, 5, NA_real_, Inf)) {
    bad <- result
    bad$local[2L] <- invalid
    expect_error(FibrTraQ:::ftq_report_local_fraction(bad),
                 "finite local counts between zero and tracks")
  }
})

test_that("fraction and validity NIfTI maps preserve the template grid and affine", {
  skip_if_not_installed("RNifti")
  root <- tempfile("public-fraction-maps-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  template_path <- file.path(root, "template.nii.gz")
  template <- RNifti::asNifti(array(0, c(3, 4, 5)), reference = list(
    sform_code = 2L, qform_code = 0L,
    srow_x = c(0, -2, 0, 12), srow_y = c(1.5, 0, 0, -4),
    srow_z = c(0, 0, 3, 7)
  ))
  RNifti::pixdim(template) <- c(1.5, 2, 3)
  RNifti::pixunits(template) <- "mm"
  RNifti::writeNifti(template, template_path)
  template <- RNifti::readNifti(template_path, internal = FALSE)
  template_hash <- tools::md5sum(template_path)
  result <- data.frame(
    fNode_Coords = c("1_1_1", "2_3_4", "3_4_5", "2_2_2"),
    tracks = c(4, 0, 8, 2), local = c(1, 0, 8, 0),
    leave_no_target = c(0, 0, 0, 2), T = c(3, 0, 0, 0), none = c(0, 1, 0, 0)
  )
  paths <- FibrTraQ:::ftq_report_write_voxel_exports(
    result, root, "fixture", template_image = template_path
  )
  maps <- lapply(paths$maps, RNifti::readNifti, internal = FALSE)
  for (map in maps) {
    expect_identical(dim(map), dim(template))
    expect_equal(RNifti::xform(map), RNifti::xform(template), tolerance = 1e-7)
    expect_equal(RNifti::pixdim(map), RNifti::pixdim(template))
    expect_equal(RNifti::pixunits(map), RNifti::pixunits(template))
  }
  expected_fraction <- expected_valid <- array(0, dim(template))
  expected_fraction[1, 1, 1] <- 0.25
  expected_fraction[3, 4, 5] <- 1
  expected_valid[1, 1, 1] <- expected_valid[3, 4, 5] <- expected_valid[2, 2, 2] <- 1
  expect_equal(as.numeric(maps$local_fraction), as.numeric(expected_fraction))
  expect_equal(as.numeric(maps$local_fraction_valid), as.numeric(expected_valid))
  fractions <- utils::read.csv(paths$tables$local_fraction)
  expect_true(is.na(fractions$local_fraction[2L]))
  expect_identical(tools::md5sum(template_path), template_hash)

  result$fNode_Coords[1L] <- "4_1_1"
  expect_error(FibrTraQ:::ftq_report_write_voxel_exports(
    result, root, "outside", template_image = template_path
  ), "outside the local-fraction template grid")
})

test_that("explicit report slices are distinct one-based indices within the grid", {
  expect_null(FibrTraQ:::ftq_report_slices(NULL, 5L))
  expect_identical(FibrTraQ:::ftq_report_slices(c(5, 1, 3), 5L), c(5L, 1L, 3L))
  for (invalid in list(0, -1, 1.5, NA_real_, Inf, numeric(), c(1, 1), "2", TRUE)) {
    expect_error(FibrTraQ:::ftq_report_slices(invalid, 5L), "positive one-based axial indices")
  }
  expect_error(FibrTraQ:::ftq_report_slices(c(1, 6), 5L), "exceeds.*axial dimension")
})

test_that("endpoint HTML uses the unambiguous paired class and its stated denominator", {
  root <- tempfile("public-endpoint-percent-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  path <- file.path(root, "report.html")
  result <- data.frame(fNode_Coords = "1_1_1", tracks = 4535, local = 0,
                       leave_no_target = 0, overlap = 4535, zero = 0, none = 0)
  endpoint <- list(
    summary = data.frame(
      target = c("overlap", "zero"), tolerance = "exact",
      target_intersection_streamlines = c(4535, 0),
      paired_focal_target_endpoints = c(1088, 0),
      shared_mask_ambiguous = c(3447, 0),
      paired_percent_of_target_intersections = c(55.3693, 99)
    ),
    context = data.frame(context_type = character(), endpoints_outside = numeric(),
                          endpoints_nonfinite = numeric()),
    penetration = data.frame(), penetration_summary = data.frame()
  )
  original <- endpoint$summary
  FibrTraQ:::ftq_write_html_report(
    path, result, rsc_summary(result, percent = TRUE, warn_legacy = FALSE),
    manifest = list(subject = "synthetic", description = "overlap percentage"),
    validation = NULL, preflight = NULL, figs = list(), top_tables = list(),
    endpoint_data = endpoint
  )
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  section <- strsplit(html, "<h3>End-Voxel Classification</h3>", fixed = TRUE)[[1L]][2L]
  table <- strsplit(section, "</table>", fixed = TRUE)[[1L]][1L]
  rows <- regmatches(table, gregexpr("<tr>.*?</tr>", table, perl = TRUE))[[1L]]
  cells <- lapply(rows, function(row) {
    values <- regmatches(row, gregexpr("<t[hd]>[^<]*</t[hd]>", row, perl = TRUE))[[1L]]
    gsub("</?t[hd]>", "", values)
  })
  percent <- match("Unambiguous paired (% of full-path intersections)", cells[[1L]])
  expect_false(is.na(percent))
  expect_equal(length(cells), 3L)
  expect_identical(cells[[2L]][percent], "23.99")
  expect_identical(cells[[3L]][percent], "Unavailable")
  expect_false(grepl("55.3693|55.37", table))
  expect_identical(endpoint$summary, original)
})

test_that("fraction montage retains defined zeros and masks unassigned voxels", {
  root <- tempfile("public-zero-fraction-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  values <- array(0, c(3, 3, 2))
  values[2, 2, 2] <- 0.5
  valid <- array(FALSE, dim(values))
  valid[2, 2, ] <- TRUE
  captured <- list()
  testthat::local_mocked_bindings(
    ftq_draw_overlay_slice = function(map_mat, x, y, zlim, pal, main, ..., show_zeros = FALSE) {
      captured[[length(captured) + 1L]] <<- list(values = map_mat, title = main, zeros = show_zeros)
      graphics::plot.new()
    },
    ftq_draw_colorbar = function(...) graphics::plot.new(),
    .package = "FibrTraQ"
  )
  FibrTraQ:::ftq_plot_montage(values, file.path(root, "fraction.png"), "Local fraction",
                             validity_mask = valid)
  expect_length(captured, 2L)
  expect_identical(vapply(captured, `[[`, "", "title"), c("z = 1", "z = 2"))
  expect_true(all(vapply(captured, `[[`, TRUE, "zeros")))
  expect_equal(captured[[1L]]$values[2, 2], 0)
  expect_equal(captured[[2L]]$values[2, 2], 0.5)
  expect_equal(sum(is.na(captured[[1L]]$values)), 8L)
  expect_equal(sum(is.na(captured[[2L]]$values)), 8L)
  expect_error(FibrTraQ:::ftq_plot_montage(values, "unused.png", "invalid",
                                         validity_mask = array(TRUE, c(3, 3, 1))),
               "validity mask must match")
  valid[1, 1, 1] <- NA
  expect_error(FibrTraQ:::ftq_plot_montage(values, "unused.png", "invalid", validity_mask = valid),
               "contain no missing values")
})

test_that("standard report forwards the saved fraction validity map to its plot", {
  skip_if_not_installed("RNifti")
  root <- tempfile("public-fraction-plot-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  template <- file.path(root, "template.nii.gz")
  RNifti::writeNifti(array(0, c(3, 3, 2)), template)
  result <- data.frame(fNode_Coords = c("2_2_1", "2_2_2"), tracks = c(4, 0),
                       local = 0, leave_no_target = c(4, 0), T = 0, none = c(0, 1))
  captured <- NULL
  testthat::local_mocked_bindings(
    ftq_plot_fate_summary = function(...) NULL,
    ftq_plot_rsc_summary = function(...) NULL,
    ftq_plot_tracks_distribution = function(...) NULL,
    ftq_plot_spatial_bundle = function(...) NULL,
    ftq_plot_montage = function(arr, path, title, ..., slices, validity_mask) {
      captured <<- list(values = as.array(arr), valid = validity_mask, slices = slices)
    }, .package = "FibrTraQ"
  )
  report <- ftq_report(result, rsc_summary(result, percent = TRUE, warn_legacy = FALSE),
    manifest = list(prefix = "fixture", fnode = template, endpoint_analysis = "none"),
    output_dir = file.path(root, "report"), report_slices = c(1L, 2L))
  expect_false(is.null(captured))
  expected <- array(FALSE, c(3, 3, 2))
  expected[2, 2, 1] <- TRUE
  expect_equal(captured$valid, expected, ignore_attr = TRUE)
  expect_identical(captured$slices, c(1L, 2L))
  expect_equal(as.numeric(captured$values), rep(0, 18L))
  expect_true(file.exists(report$html))
})

test_that("the fraction overlay draws defined zeros without coloring missing values", {
  layers <- list()
  draw_overlay <- ftq_test_function_copy(FibrTraQ:::ftq_draw_overlay_slice, exports = list(
    `graphics::image` = function(x, y, z, ...) {
      layers[[length(layers) + 1L]] <<- z
    }, `graphics::box` = function(...) NULL
  ))
  values <- matrix(c(0, NA, 0.5, -1), 2L)
  draw_overlay(values, 1:2, 1:2, c(0, 1), c("blue", "red"), "fraction", show_zeros = TRUE)
  expect_length(layers, 2L)
  expect_equal(as.numeric(layers[[2L]]), c(0, NA, 0.5, NA))
  layers <- list()
  draw_overlay(values, 1:2, 1:2, c(0, 1), c("blue", "red"), "counts")
  expect_equal(as.numeric(layers[[2L]]), c(NA, NA, 0.5, NA))
})

test_that("a returned report failure preserves completed scientific files and validation", {
  skip_if_not_installed("RNifti")
  root <- tempfile("public-report-failure-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  fixture <- ftq_example_data()
  captured_hashes <- NULL
  captured_slices <- NULL
  testthat::local_mocked_bindings(
    ftq_report = function(result_csv, report_slices, ...) {
      captured_hashes <<- tools::md5sum(list.files(dirname(result_csv), full.names = TRUE))
      captured_slices <<- report_slices
      stop("injected presentation failure")
    }, .package = "FibrTraQ"
  )
  manifest <- list(
    dataset_root = fixture$root, subject = "001", description = "report-failure",
    fnode = fixture$focal_mask, local = fixture$local_mask,
    tractogram = fixture$tractogram, ref_image = fixture$reference_image,
    targets = as.list(fixture$target_masks),
    endpoint_contexts = as.list(fixture$endpoint_contexts),
    coordinate = list(track_space = "voxel", discretization = "round", axis_flip = rep(FALSE, 3L)),
    analysis = list(focal_assignment_mode = "all_focal_voxels",
      streamline_voxelization_mode = "sample_points", target_intersection_mode = "full_streamline",
      endpoint_analysis = "exact_and_dilated1"),
    report = list(context_source = "none"),
    outputs = list(output_dir = root, write_maps = FALSE, write_report = TRUE, overwrite = FALSE)
  )
  run <- ftq_run(manifest, backend = "r", progress = FALSE,
                 report_slices = c(1L, 2L), report_error_action = "return")
  expect_null(run$report)
  expect_identical(run$report_status, "failed")
  expect_identical(run$report_error, "injected presentation failure")
  expect_identical(captured_slices, c(1L, 2L))
  expect_true(run$validation$passed)
  expect_true(length(captured_hashes) >= 7L)
  expect_identical(tools::md5sum(names(captured_hashes)), captured_hashes)
  expect_equal(FibrTraQ:::ftq_read_table(run$outputs$result_csv),
               FibrTraQ:::ftq_read_table(fixture$expected[["result"]]))
  saved_validation <- jsonlite::read_json(run$outputs$validation_json, simplifyVector = TRUE)
  expect_true(saved_validation$passed)
})
