tck_report_manifest <- function(paths, output_dir, ...) {
  modifyList(list(
    subject = "synthetic", description = "tck-reports",
    tractogram = paths$tck, fnode = paths$focal,
    targets = as.list(paths$targets),
    coordinate = list(track_space = "auto", discretization = "round"),
    execution = list(track_chunk_size = 2L),
    outputs = list(output_dir = output_dir)
  ), list(...))
}

tck_report_scene <- function(path) {
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  pattern <- '<script id="scene-data" type="application/json">(.*?)</script>'
  match <- regmatches(html, regexec(pattern, html, perl = TRUE))[[1L]]
  stopifnot(length(match) == 2L)
  jsonlite::fromJSON(match[[2L]], simplifyVector = FALSE)
}

test_that("native TCK creates every report and matches equivalent TrackVis outputs", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("freesurferformats")
  affine <- diag(c(-1.25, 2, 0.75, 1))
  affine[1:3, 4] <- c(15, -8, 2)
  fixture <- tck_fixture_case(affine)
  root <- tempfile("ftq-")
  paths <- tck_fixture_files(file.path(root, "inputs"), fixture)
  variants <- list(tck_stream = c("tck", "stream"),
                   tck_memory = c("tck", "memory"),
                   trk_stream = c("trk", "stream"))
  bundles <- lapply(names(variants), function(name) {
    variant <- variants[[name]]
    ftq_run_files(
      tractogram = paths[[variant[[1L]]]], focal_mask = paths$focal,
      target_masks = paths$targets, output_dir = file.path(root, name),
      subject = "synthetic", focal_label = "focal",
      # Omitting anatomy tests native TCK's real focal-image reference fallback.
      track_space = if (variant[[1L]] == "tck") "auto" else "voxelmm",
      axis_flip = if (variant[[1L]] == "tck") NULL else rep(FALSE, 3L),
      discretization = "round", track_reader = variant[[2L]],
      track_chunk_size = 2L, sample_per_target = 10L,
      streamline_voxelization_mode = "segment_interpolation",
      backend = "r", verbose = FALSE, progress = FALSE
    )
  })
  names(bundles) <- names(variants)
  for (bundle in bundles) {
    expect_identical(bundle$index$run_status, "completed")
    expect_true(bundle$index$validation_passed)
    expect_true(all(file.exists(unlist(bundle$index[c(
      "primary_report", "exploratory_report", "interactive_3d_qc",
      "endpoint_summary_csv", "metadata_json"
    )]))))
    expect_true(file.exists(bundle$report_index))
    endpoint <- utils::read.csv(bundle$index$endpoint_summary_csv)
    spacing_columns <- c("voxel_size_x_mm", "voxel_size_y_mm", "voxel_size_z_mm")
    expected_spacing <- c(1.25, 2, 0.75)
    for (i in seq_along(spacing_columns)) {
      expect_equal(endpoint[[spacing_columns[[i]]]], rep(expected_spacing[[i]], nrow(endpoint)))
    }
    provenance <- jsonlite::fromJSON(bundle$index$metadata_json)
    expect_equal(provenance$endpoint_analysis$reference_voxel_size_mm, expected_spacing)
  }
  reference <- bundles$tck_stream
  run_dir <- function(bundle) dirname(bundle$index$result_csv)
  for (bundle in bundles[-1L]) {
    # Compare every exported scientific CSV, including both report tables and
    # endpoint audits, not only the top-level count table.
    expected_csv <- list.files(run_dir(reference), "[.]csv$", recursive = TRUE)
    actual_csv <- list.files(run_dir(bundle), "[.]csv$", recursive = TRUE)
    expect_setequal(actual_csv, expected_csv)
    for (path in expected_csv) {
      actual <- utils::read.csv(file.path(run_dir(bundle), path), check.names = FALSE)
      expected <- utils::read.csv(file.path(run_dir(reference), path), check.names = FALSE)
      if (grepl("_roi-summary.csv$", path) &&
          identical(bundle$index$selected_track_space, "voxelmm")) {
        # Input-format-specific coordinate provenance must remain distinct;
        # every scientific count, denominator and metric must still match.
        expect_identical(actual$selected_track_space, "voxelmm")
        expect_identical(expected$selected_track_space, "world")
        expect_identical(actual$mapping_selection_rule, "pinned")
        expect_identical(expected$mapping_selection_rule,
                         "tck_world_ras_mm_from_nifti_affine")
        format_columns <- c("selected_track_space", "track_voxel_size_mm",
                            "mapping_selection_rule")
        actual[format_columns] <- expected[format_columns] <- NULL
      }
      expect_equal(actual, expected, info = path)
    }
    expected_maps <- list.files(run_dir(reference), "[.]nii[.]gz$", recursive = TRUE)
    actual_maps <- list.files(run_dir(bundle), "[.]nii[.]gz$", recursive = TRUE)
    expect_setequal(actual_maps, expected_maps)
    expect_true(any(grepl("tractogram-density", expected_maps, fixed = TRUE)))
    for (path in expected_maps) {
      expected <- RNifti::readNifti(file.path(run_dir(reference), path))
      actual <- RNifti::readNifti(file.path(run_dir(bundle), path))
      expect_equal(as.numeric(actual), as.numeric(expected), info = path)
      expect_equal(RNifti::xform(actual), RNifti::xform(expected), info = path)
    }
    expect_equal(tck_report_scene(bundle$index$interactive_3d_qc),
                 tck_report_scene(reference$index$interactive_3d_qc),
                 tolerance = 1e-6)
  }
  metadata <- jsonlite::fromJSON(reference$index$metadata_json)
  expect_identical(metadata$manifest$tractogram_format, "tck")
  expect_identical(metadata$manifest$track_space, "world")
  expect_identical(metadata$manifest$axis_flip, rep(FALSE, 3L))
  expect_identical(metadata$manifest$coordinate_mapping$input_coordinates,
                   "RAS+ world millimetres")
  expect_identical(metadata$files$context$metadata$tractogram_format, "tck")
  html <- paste(readLines(reference$index$primary_report), collapse = "\n")
  expect_match(html, "native TCK RAS+ world millimeters", fixed = TRUE)
  expect_match(html, "no axis flips or overlap-based orientation selection", fixed = TRUE)
})

test_that("standalone TCK density and viewer enforce the declared world geometry", {
  skip_if_not_installed("RNifti")
  root <- tempfile("tck-standalone-reports-")
  paths <- tck_fixture_files(file.path(root, "inputs"))
  context <- ftq_tractogram_context(
    paths$tck, paths$focal, file.path(root, "context"),
    track_chunk_size = 2L, density_discretization = "round"
  )
  expect_true(file.exists(context$density_image))
  expect_equal(context$summary$processed_streamlines, 6)
  expect_equal(context$summary$declared_streamlines, 6)
  expect_equal(context$summary$mapping_error_streamlines, 0)
  expect_identical(context$metadata$selected_track_space, "world")
  expect_identical(context$metadata$axis_flip, rep(FALSE, 3L))
  expect_error(ftq_tractogram_context(
    paths$tck, paths$focal, file.path(root, "bad-context"),
    axis_flip = c(TRUE, FALSE, FALSE)
  ), "axis_flip")
  manifest <- tck_report_manifest(paths, file.path(root, "viewer"))
  viewer <- ftq_interactive_qc_viewer(manifest, sample_per_target = 10L)
  expect_identical(viewer$track_space, "world")
  expect_true(file.exists(viewer$html))
  for (assignment in c("endpoint_voxels", "local_search")) {
    tck_manifest <- manifest
    tck_manifest$analysis <- list(target_intersection_mode = assignment)
    tck_viewer <- ftq_interactive_qc_viewer(
      tck_manifest, sample_per_target = 10L,
      output_html = file.path(root, paste0("tck-", assignment, ".html"))
    )
    trk_manifest <- tck_manifest
    trk_manifest$tractogram <- paths$trk
    trk_manifest$coordinate$track_space <- "voxelmm"
    trk_manifest$coordinate$axis_flip <- rep(FALSE, 3L)
    trk_viewer <- ftq_interactive_qc_viewer(
      trk_manifest, sample_per_target = 10L,
      output_html = file.path(root, paste0("trk-", assignment, ".html"))
    )
    expect_equal(tck_viewer$summary, trk_viewer$summary, info = assignment)
    expect_equal(tck_report_scene(tck_viewer$html), tck_report_scene(trk_viewer$html),
                 tolerance = 1e-6, info = assignment)
  }
  manifest$coordinate$track_space <- "voxel"
  expect_error(ftq_interactive_qc_viewer(
    manifest, output_html = file.path(root, "bad-viewer.html")
  ), "track_space must be 'auto' or 'world'", fixed = TRUE)
})

test_that("TCK preflight blocks invalid headers and ambiguous image geometry", {
  skip_if_not_installed("RNifti")
  root <- tempfile("tck-preflight-")
  paths <- tck_fixture_files(file.path(root, "inputs"))
  manifest <- tck_report_manifest(paths, file.path(root, "out"))
  preflight <- ftq_preflight(manifest)
  expect_true(preflight$passed)
  expect_true(all(preflight$checks$passed[grepl("^tck_", preflight$checks$check)]))
  expect_false("ref_image_recommended_for_auto_mapping" %in% preflight$checks$check)
  wrong_space <- manifest
  wrong_space$coordinate$track_space <- "voxelmm"
  expect_false(ftq_preflight(wrong_space)$passed)
  wrong_flip <- manifest
  wrong_flip$coordinate$axis_flip <- c(FALSE, TRUE, FALSE)
  expect_false(ftq_preflight(wrong_flip)$passed)
  bad_header <- manifest
  bad_header$tractogram <- file.path(root, "invalid.tck")
  writeLines("not a tractogram", bad_header$tractogram)
  failed <- ftq_preflight(bad_header)
  expect_false(failed$passed)
  expect_false(failed$checks$passed[failed$checks$check == "tck_header_readable"])
  bad_grid <- manifest
  unknown <- tck_fixture_case()
  RNifti::pixunits(unknown$targets$target_a) <- "unknown"
  unknown_path <- file.path(root, "unknown-units.nii.gz")
  RNifti::writeNifti(unknown$targets$target_a, unknown_path)
  bad_grid$targets$target_a <- unknown_path
  failed <- ftq_preflight(bad_grid)
  expect_false(failed$passed)
  expect_false(failed$checks$passed[
    failed$checks$check == "tck_spatial_geometry:target:target_a"
  ])
  expect_error(ftq_run(bad_grid), "preflight failed")
  expect_false(dir.exists(manifest$outputs$output_dir))
})


test_that("TrackVis files without a reference retain focal physical spacing in endpoint audits", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("freesurferformats")
  spacing_mm <- c(1.25, 2, 0.75)
  affine <- diag(c(spacing_mm, 1))
  for (units in c("mm", "m", "um")) {
    root <- tempfile(paste0("trk-no-reference-spacing-", units, "-"))
    paths <- tck_fixture_files(file.path(root, "inputs"),
                               tck_fixture_case(affine, units = units))
    manifest <- tck_report_manifest(
      paths, file.path(root, "out"), tractogram = paths$trk,
      coordinate = list(track_space = "voxelmm", discretization = "round",
                        axis_flip = rep(FALSE, 3L)),
      outputs = list(output_dir = file.path(root, "out"), write_maps = FALSE,
                     write_report = FALSE)
    )
    run <- ftq_run(manifest, backend = "r")
    expect_null(run$manifest$ref_image)
    expect_identical(attr(run$result, "fibrtraq_track_space"), "voxelmm")
    expect_identical(attr(run$result, "fibrtraq_axis_flip"), rep(FALSE, 3L))
    ftq_expect_core_equal(run$result, tck_fixture_expected())
    expect_equal(run$endpoint$metadata$reference_voxel_size_mm, spacing_mm,
                 tolerance = 1e-6, info = units)
    summary <- utils::read.csv(run$outputs$endpoint_summary_csv)
    for (axis in seq_along(spacing_mm)) {
      expect_equal(summary[[paste0("voxel_size_", c("x", "y", "z")[[axis]], "_mm")]],
                   rep(spacing_mm[[axis]], nrow(summary)), tolerance = 1e-6, info = units)
    }
    metadata <- jsonlite::fromJSON(run$outputs$metadata_json)
    expect_equal(metadata$endpoint_analysis$reference_voxel_size_mm, spacing_mm,
                 tolerance = 1e-6, info = units)
  }
})
