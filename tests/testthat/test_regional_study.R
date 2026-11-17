library(testthat)
library(FibrTraQ)

ftq_write_regional_study_test_trk <- function(path,
                                               tracks,
                                               dims = c(5L, 5L, 5L)) {
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

ftq_regional_study_fixture <- function(root) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  mask_a <- file.path(root, "A.nii.gz")
  mask_b <- file.path(root, "B.nii.gz")
  reference <- file.path(root, "dwi.nii.gz")
  bvals <- file.path(root, "dwi.bval")
  tractogram <- file.path(root, "tracks.trk")

  a <- array(0, dim = c(5, 5, 5))
  b <- a
  a[2, 2, 2] <- 1
  b[4, 2, 2] <- 1
  RNifti::writeNifti(a, mask_a)
  RNifti::writeNifti(b, mask_b, template = mask_a)

  dwi <- array(0, dim = c(5, 5, 5, 2))
  dwi[, , , 1] <- array(seq_len(125), dim = c(5, 5, 5))
  dwi[, , , 2] <- dwi[, , , 1] * 0.5
  RNifti::writeNifti(dwi, reference, template = mask_a)
  writeLines("0 1000", bvals)
  ftq_write_regional_study_test_trk(
    tractogram,
    list(
      rbind(c(2, 2, 2), c(3, 2, 2), c(4, 2, 2)),
      rbind(c(2, 2, 2), c(2, 3, 2))
    )
  )

  list(
    tractogram = tractogram,
    masks = c(A = mask_a, B = mask_b),
    reference = reference,
    bvals = bvals,
    expected_background = dwi[, , , 1]
  )
}

test_that("regional study writes one report with distinct results for two synthetic subjects", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("freesurferformats")

  fixture <- ftq_regional_study_fixture(
    file.path(tempdir(), "regional-study-complete-fixture")
  )
  output_dir <- file.path(tempdir(), "regional-output")
  result <- fibrtraq_study(
    tractograms = c(subject1 = fixture$tractogram, subject2 = fixture$tractogram),
    region_masks = list(subject1 = fixture$masks,
                        subject2 = setNames(rev(fixture$masks), names(fixture$masks))),
    output_dir = output_dir,
    focal_regions = "A",
    reference_images = c(subject1 = fixture$reference, subject2 = fixture$reference),
    bvals = c(subject1 = fixture$bvals, subject2 = fixture$bvals),
    local_masks = "union",
    description = "integration-test",
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    streamline_voxelization_mode = "segment_interpolation",
    backend = "r",
    sample_per_target = 1L,
    overwrite = TRUE,
    verbose = FALSE,
    progress = FALSE
  )

  expect_s3_class(result, "ftq_regional_study_result")
  expect_equal(nrow(result$index), 2L)
  expect_identical(result$index$subject, c("subject1", "subject2"))
  expect_identical(result$index$focal_region, c("A", "A"))
  # The second subject uses the other fixture mask as focal, so linking either
  # subject to the first subject's tables would return the wrong counts.
  hits <- vapply(result$index$roi_summary_csv, function(path) {
    utils::read.csv(path)$unique_focal_intersecting_streamlines
  }, numeric(1), USE.NAMES = FALSE)
  expect_equal(hits, c(2, 1))
  expect_true(all(result$index$run_status == "completed"))
  expect_true(all(result$index$validation_passed))
  expect_true(all(file.exists(result$index$primary_report)))
  expect_true(all(file.exists(result$index$exploratory_report)))
  expect_true(all(file.exists(result$index$interactive_3d_qc)))
  expect_true(all(file.exists(result$index$endpoint_summary_csv)))
  expect_true(file.exists(result$index_csv))
  expect_true(file.exists(result$prepared_inputs_csv))
  expect_identical(normalizePath(result$report_index), normalizePath(file.path(output_dir, "fibrtraq-report-index.html")))
  expect_length(list.files(output_dir, "[.]html$", recursive = TRUE), 1L)
  for (i in seq_len(nrow(result$index))) {
    focal_dir <- file.path(output_dir, result$index$subject[i], result$index$focal_region[i])
    expect_identical(result$index$primary_report[i], result$report_index)
    expect_identical(result$index$exploratory_report[i], result$report_index)
    expect_identical(result$index$interactive_3d_qc[i], result$report_index)
    expect_false(dir.exists(file.path(focal_dir, "reports")))
    expect_false(dir.exists(file.path(focal_dir, "exploratory-report")))
    expect_length(list.files(focal_dir, "[.]csv$"), 0L)
    expect_length(list.files(file.path(focal_dir, "tables"), "_exploratory-primary-target-lut[.]csv$"), 1L)
    expect_length(list.files(file.path(focal_dir, "figures"), "[.]png$"), 15L)
    expect_length(list.files(file.path(focal_dir, "figures"), "_qc-scene[.]js$"), 1L)
    expect_length(list.files(file.path(focal_dir, "maps"), "_local-fraction(-valid)?[.]nii[.]gz$"), 0L)
  }

  prepared <- result$prepared_inputs[1, ]
  expect_true(file.exists(prepared$reference_image))
  expect_true(file.exists(prepared$local_mask))
  background <- as.array(RNifti::readNifti(
    prepared$reference_image,
    internal = FALSE
  ))
  expect_identical(dim(background), dim(fixture$expected_background))
  expect_equal(
    as.numeric(background),
    as.numeric(fixture$expected_background),
    tolerance = 1e-6
  )

  index_html <- paste(readLines(result$report_index, warn = FALSE), collapse = "\n")
  expect_match(index_html, "subject1/A/tables/", fixed = TRUE)
  expect_match(index_html, "subject2/A/tables/", fixed = TRUE)
  expect_match(index_html, "id=['\"]ftq-run['\"]")
  expect_match(index_html, "id=['\"]ftq-primary['\"]")
  expect_match(index_html, "id=['\"]ftq-targets['\"]")
  expect_match(index_html, "id=['\"]ftq-open-viewer['\"]")
})

test_that("standard runs omit unused preparation and context folders but retain density provenance", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("freesurferformats")
  root <- tempfile("f-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  fixture <- ftq_regional_study_fixture(file.path(root, "inputs"))
  output_dir <- file.path(root, "results")
  result <- fibrtraq_single(
    tractogram = fixture$tractogram,
    focal_mask = fixture$masks[["A"]],
    target_masks = fixture$masks["B"],
    output_dir = output_dir,
    subject = "subject1", focal_label = "A",
    track_space = "voxel", discretization = "round", axis_flip = rep(FALSE, 3),
    backend = "r", write_exploratory_report = FALSE, write_qc_viewer = FALSE,
    verbose = FALSE, progress = FALSE
  )
  focal_dir <- file.path(output_dir, "subject1", "A")
  expect_false(dir.exists(file.path(output_dir, "subject1", "prepared")))
  expect_false(dir.exists(file.path(focal_dir, "context")))
  density_paths <- list.files(file.path(focal_dir, "maps"),
                              pattern = "_tractogram-density[.]nii[.]gz$", full.names = TRUE)
  expect_length(density_paths, 1L)
  expect_length(list.files(focal_dir, pattern = "tractogram-context", recursive = TRUE), 0L)
  metadata <- jsonlite::read_json(result$index$metadata_json, simplifyVector = TRUE)
  context <- metadata$files$context
  expect_identical(normalizePath(context$density_image), normalizePath(density_paths[[1L]]))
  expect_null(context$summary_csv)
  expect_null(context$metadata_json)
  expect_equal(context$summary$processed_streamlines, 2L)
  expect_identical(context$metadata$selected_track_space, "voxel")
  expect_identical(context$metadata$streamline_voxelization_mode, "sample_points")
  expect_identical(context$metadata$tractogram, normalizePath(fixture$tractogram))
})

test_that("aligned existing images and supplied local masks do not create prepared derivatives", {
  skip_if_not_installed("RNifti")
  root <- tempfile("regional-no-preparation-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  fixture <- ftq_regional_study_fixture(file.path(root, "inputs"))
  output_dir <- file.path(root, "results")
  result <- fibrtraq_single(
    tractogram = fixture$tractogram,
    focal_mask = fixture$masks[["A"]], target_masks = fixture$masks["B"],
    reference_image = fixture$masks[["A"]], local_mask = fixture$masks[["A"]],
    output_dir = output_dir, subject = "subject1", focal_label = "A",
    track_space = "voxel", discretization = "round", axis_flip = rep(FALSE, 3),
    backend = "r", write_report = FALSE, write_exploratory_report = FALSE,
    write_qc_viewer = FALSE, verbose = FALSE, progress = FALSE
  )
  expect_true(result$index$validation_passed)
  expect_false(dir.exists(file.path(output_dir, "subject1", "prepared")))
  expect_identical(result$prepared_inputs$reference_image, normalizePath(fixture$masks[["A"]]))
  expect_identical(result$prepared_inputs$local_mask, normalizePath(fixture$masks[["A"]]))
})

test_that("file runner provides the complete single-analysis workflow", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("freesurferformats")

  fixture <- ftq_regional_study_fixture(
    file.path(tempdir(), "file-run-complete-fixture")
  )
  output_dir <- file.path(tempdir(), "file-run-complete-output")
  result <- fibrtraq_single(
    tractogram = fixture$tractogram,
    focal_mask = fixture$masks[["A"]],
    target_masks = fixture$masks["B"],
    output_dir = output_dir,
    subject = "subject1",
    focal_label = "A",
    reference_image = fixture$reference,
    bvals = fixture$bvals,
    description = "file-run-test",
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    streamline_voxelization_mode = "segment_interpolation",
    backend = "r",
    sample_per_target = 1L,
    overwrite = TRUE,
    verbose = FALSE,
    progress = FALSE
  )

  expect_s3_class(result, "ftq_regional_study_result")
  expect_equal(nrow(result$index), 1L)
  expect_equal(result$index$subject, "subject1")
  expect_equal(result$index$focal_region, "A")
  expect_length(list.files(output_dir, "[.]html$", recursive = TRUE), 1L)
  expect_identical(normalizePath(dirname(result$report_index)),
                   normalizePath(file.path(output_dir, "subject1", "A")))
  expect_length(list.files(output_dir, "^OUTPUTS_README", recursive = TRUE), 0L)
  expect_equal(result$index$run_status, "completed")
  expect_true(result$index$validation_passed)
  expect_true(file.exists(result$index$primary_report))
  expect_true(file.exists(result$index$exploratory_report))
  expect_true(file.exists(result$index$interactive_3d_qc))
  expect_true(file.exists(result$index$endpoint_summary_csv))
  expect_true(file.exists(result$report_index))
})

test_that("file runner requires named target masks", {
  fixture <- ftq_regional_study_fixture(
    file.path(tempdir(), "file-run-target-label-fixture")
  )
  expect_error(
    fibrtraq_single(
      tractogram = fixture$tractogram,
      focal_mask = fixture$masks[["A"]],
      target_masks = unname(fixture$masks["B"]),
      output_dir = file.path(tempdir(), "file-run-target-label-output"),
      write_report = FALSE,
      write_exploratory_report = FALSE,
      write_qc_viewer = FALSE,
      verbose = FALSE,
      progress = FALSE
    ),
    "target_masks must be a non-empty named character vector",
    fixed = TRUE
  )
})

test_that("regional study runner requires b-values for a 4-D reference", {
  skip_if_not_installed("RNifti")

  fixture <- ftq_regional_study_fixture(
    file.path(tempdir(), "regional-study-missing-bvals-fixture")
  )
  expect_error(
    fibrtraq_study(
      tractograms = c(subject1 = fixture$tractogram),
      region_masks = list(subject1 = fixture$masks),
      output_dir = file.path(tempdir(), "regional-study-missing-bvals-output"),
      reference_images = c(subject1 = fixture$reference),
      local_masks = NULL,
      track_space = "voxel",
      discretization = "round",
      axis_flip = c(FALSE, FALSE, FALSE),
      write_maps = FALSE,
      write_report = FALSE,
      write_exploratory_report = FALSE,
      write_qc_viewer = FALSE,
      overwrite = TRUE,
      verbose = FALSE,
      progress = FALSE
    ),
    "four-dimensional; supply matching bvals",
    fixed = TRUE
  )
})

test_that("regional study runner makes affine reheadering explicit", {
  skip_if_not_installed("RNifti")

  fixture <- ftq_regional_study_fixture(
    file.path(tempdir(), "regional-study-affine-fixture")
  )
  reference <- RNifti::readNifti(fixture$reference, internal = TRUE)
  reflected <- diag(4)
  reflected[1, 1] <- -1
  reflected[1, 4] <- 4
  attr(reflected, "code") <- 1L
  RNifti::qform(reference) <- reflected
  RNifti::sform(reference) <- reflected
  RNifti::writeNifti(reference, fixture$reference)

  common_args <- list(
    tractograms = c(subject1 = fixture$tractogram),
    region_masks = list(subject1 = fixture$masks),
    reference_images = c(subject1 = fixture$reference),
    bvals = c(subject1 = fixture$bvals),
    focal_regions = "A",
    local_masks = NULL,
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    endpoint_analysis = "none",
    write_maps = FALSE,
    write_report = FALSE,
    write_exploratory_report = FALSE,
    write_qc_viewer = FALSE,
    backend = "r",
    overwrite = TRUE,
    verbose = FALSE,
    progress = FALSE
  )
  expect_error(
    do.call(
      fibrtraq_study,
      c(
        common_args,
        list(output_dir = file.path(tempdir(), "regional-study-affine-rejected"))
      )
    ),
    "set reheader_images_to_region_grid = TRUE",
    fixed = TRUE
  )

  accepted <- do.call(
    fibrtraq_study,
    c(
      common_args,
      list(
        output_dir = file.path(tempdir(), "regional-study-affine-accepted"),
        reheader_images_to_region_grid = TRUE
      )
    )
  )
  expect_true(accepted$index$validation_passed)
  prepared <- RNifti::readNifti(
    accepted$prepared_inputs$reference_image,
    internal = TRUE
  )
  mask <- RNifti::readNifti(fixture$masks[["A"]], internal = TRUE)
  expect_true(ftq_assert_same_grid(prepared, mask))
  provenance <- utils::read.csv(sub(
    "\\.nii\\.gz$",
    "-volumes.csv",
    accepted$prepared_inputs$reference_image
  ))
  expect_true(all(provenance$reheadered_to_region_grid))
  expect_equal(unique(provenance$preparation), "mean_low_b")
})

test_that("file runner reheaders 3-D images fully to the regional grid", {
  skip_if_not_installed("RNifti")

  root <- file.path(tempdir(), "file-runner-3d-reheader-fixture")
  fixture <- ftq_regional_study_fixture(root)
  reference_path <- file.path(root, "anatomy.nii.gz")
  values <- array(seq_len(125), dim = c(5, 5, 5))
  reference <- RNifti::asNifti(values)
  RNifti::pixunits(reference) <- "mm"
  source_affine <- diag(c(0.7, 0.8, 0.9, 1))
  source_affine[1:3, 4] <- c(2, 3, 4)
  attr(source_affine, "code") <- 1L
  RNifti::qform(reference) <- source_affine
  RNifti::sform(reference) <- source_affine
  RNifti::writeNifti(reference, reference_path)

  result <- fibrtraq_single(
    tractogram = fixture$tractogram,
    focal_mask = fixture$masks[["A"]],
    target_masks = fixture$masks["B"],
    output_dir = file.path(tempdir(), "file-runner-3d-reheader-output"),
    subject = "subject1",
    focal_label = "A",
    reference_image = reference_path,
    reheader_image_to_mask_grid = TRUE,
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    endpoint_analysis = "none",
    backend = "r",
    write_maps = FALSE,
    write_report = FALSE,
    write_exploratory_report = FALSE,
    write_qc_viewer = FALSE,
    overwrite = TRUE,
    verbose = FALSE,
    progress = FALSE
  )

  prepared <- RNifti::readNifti(
    result$prepared_inputs$reference_image,
    internal = TRUE
  )
  mask <- RNifti::readNifti(fixture$masks[["A"]], internal = TRUE)
  expect_true(ftq_assert_same_grid(prepared, mask))
  expect_equal(as.character(RNifti::pixunits(prepared))[[1L]], "mm")
  expect_identical(
    as.numeric(RNifti::readNifti(
      result$prepared_inputs$reference_image,
      internal = FALSE
    )),
    as.numeric(values)
  )
  provenance <- utils::read.csv(sub(
    "\\.nii\\.gz$",
    "-volumes.csv",
    result$prepared_inputs$reference_image
  ))
  expect_equal(unique(provenance$preparation), "reheader_without_resampling")
  expect_true(all(provenance$reheadered_to_region_grid))
})

ftq_public_regional_args <- function(fixture, output_dir, ...) {
  c(list(
    tractograms = c(`001` = fixture$tractogram),
    region_masks = list(`001` = fixture$masks),
    output_dir = output_dir, focal_regions = "A", local_masks = NULL,
    track_space = "voxel", discretization = "round", axis_flip = rep(FALSE, 3),
    backend = "r", write_maps = FALSE, write_report = TRUE,
    write_exploratory_report = TRUE, write_qc_viewer = TRUE,
    overwrite = TRUE, verbose = FALSE, progress = FALSE
  ), list(...))
}

test_that("presentation failures preserve completed data and independent stage status", {
  fixture <- ftq_regional_study_fixture(tempfile("public-stage-fixture-"))
  output_dir <- tempfile("public-stage-results-")
  analyses <- 0L
  testthat::local_mocked_bindings(
    ftq_run = function(manifest, ...) {
      analyses <<- analyses + 1L
      dir.create(manifest$outputs$output_dir, recursive = TRUE)
      csv <- file.path(manifest$outputs$output_dir, "complete-data.csv")
      writeLines("fNode_Coords,tracks\n2;2;2,2", csv)
      list(manifest = manifest, outputs = list(result_csv = csv),
           roi_summary = data.frame(processed_streamlines = 2,
             unique_focal_intersecting_streamlines = 2),
           validation = list(passed = TRUE, biological_qc_status = "not_reviewed"),
           report = NULL, report_status = "failed", report_error = "figure generation failed")
    },
    ftq_exploratory_report = function(...) stop("exploratory failed"),
    ftq_interactive_qc_viewer = function(run, ...) {
      scene <- file.path(dirname(run$outputs$result_csv), "figures", "test_qc-scene.js")
      dir.create(dirname(scene), showWarnings = FALSE)
      writeLines('document.currentScript.ftqScene = {};', scene)
      list(scene_file = scene)
    },
    .package = "FibrTraQ"
  )
  result <- do.call(fibrtraq_study,
                   ftq_public_regional_args(fixture, output_dir, continue_on_error = TRUE))
  row <- result$index[1L, ]
  expect_equal(analyses, 1L)
  expect_identical(row$subject, "001")
  expect_identical(row$run_status, "failed")
  expect_identical(row$analysis_status, "completed")
  expect_identical(row$validation_status, "passed")
  expect_true(row$validation_passed)
  expect_identical(row$biological_qc_status, "not_reviewed")
  expect_identical(row$primary_report_status, "failed")
  expect_identical(row$exploratory_report_status, "failed")
  expect_identical(row$interactive_3d_qc_status, "completed")
  expect_equal(row$processed_streamlines, 2)
  expect_true(file.exists(row$result_csv))
  expect_true(file.exists(row$interactive_3d_qc))
  expect_match(row$primary_report_error, "figure generation failed", fixed = TRUE)
  html <- paste(readLines(result$report_index), collapse = "\n")
  expect_match(html, "href='complete-data.csv'", fixed = TRUE)
  expect_match(html, "Per-voxel data", fixed = TRUE)
  expect_match(html, "figure generation failed", fixed = TRUE)
  expect_match(html, "exploratory failed", fixed = TRUE)
  expect_match(html, "3D QC", fixed = TRUE)
  expect_match(html, ">001 / A</h1>", fixed = TRUE)
  expect_false(grepl("Regional Study", html, fixed = TRUE))

  stopped_dir <- tempfile("public-report-default-stop-")
  expect_error(do.call(fibrtraq_study,
                      ftq_public_regional_args(fixture, stopped_dir)),
               "figure generation failed", fixed = TRUE)
  saved <- read.csv(file.path(stopped_dir, "fibrtraq-run-index.csv"),
                    colClasses = c(subject = "character"))
  expect_identical(saved$analysis_status, "completed")
  expect_true(saved$validation_passed)
  expect_true(file.exists(saved$result_csv))
  expect_length(list.files(file.path(stopped_dir, "001", "A"), "_report[.]html$"), 1L)
})

test_that("default stop persists failed run and earlier completed run in navigable index", {
  fixture <- ftq_regional_study_fixture(tempfile("public-stop-fixture-"))
  output_dir <- tempfile("public-stop-results-")
  calls <- 0L
  testthat::local_mocked_bindings(
    ftq_run = function(manifest, ...) {
      calls <<- calls + 1L
      if (calls == 2L) stop("counting failed")
      dir.create(manifest$outputs$output_dir, recursive = TRUE)
      csv <- file.path(manifest$outputs$output_dir, "complete-data.csv")
      writeLines("complete scientific output", csv)
      list(manifest = manifest, outputs = list(result_csv = csv),
           validation = list(passed = TRUE, biological_qc_status = "unavailable"))
    },
    .package = "FibrTraQ"
  )
  args <- ftq_public_regional_args(fixture, output_dir)
  args$focal_regions <- c("A", "B")
  args$write_report <- args$write_exploratory_report <- args$write_qc_viewer <- FALSE
  expect_error(do.call(fibrtraq_study, args), "counting failed", fixed = TRUE)
  csv <- file.path(output_dir, "fibrtraq-run-index.csv")
  html <- file.path(output_dir, "fibrtraq-report-index.html")
  expect_true(file.exists(csv))
  expect_true(file.exists(html))
  saved <- read.csv(csv, colClasses = c(subject = "character"))
  expect_identical(saved$subject, c("001", "001"))
  expect_identical(saved$analysis_status, c("completed", "failed"))
  expect_identical(saved$validation_status, c("passed", "not_run"))
  expect_identical(saved$primary_report_status, c("not_requested", "not_requested"))
  expect_true(file.exists(saved$result_csv[1L]))
  expect_true(is.na(saved$result_csv[2L]))
  text <- paste(readLines(html), collapse = "\n")
  expect_match(text, "complete-data.csv", fixed = TRUE)
  expect_match(text, "counting failed", fixed = TRUE)
})

test_that("a failed final report assembly retains validated results and clears report links", {
  fixture <- ftq_regional_study_fixture(tempfile("assembly-fixture-"))
  root <- tempfile("assembly-results-")
  args <- ftq_public_regional_args(fixture, root)
  args$write_exploratory_report <- FALSE
  args$write_qc_viewer <- FALSE
  testthat::local_mocked_bindings(
    ftq_write_unified_report = function(...) stop("simulated write failure"),
    .package = "FibrTraQ"
  )
  expect_error(do.call(fibrtraq_study, args), "simulated write failure")
  row <- utils::read.csv(file.path(root, "fibrtraq-run-index.csv"))
  expect_identical(row$run_status, "failed")
  expect_identical(row$analysis_status, "completed")
  expect_identical(row$validation_status, "passed")
  expect_true(file.exists(row$result_csv))
  expect_true(row$validation_passed)
  expect_identical(row$primary_report_status, "failed")
  expect_identical(row$exploratory_report_status, "not_requested")
  expect_identical(row$interactive_3d_qc_status, "not_requested")
  expect_match(row$primary_report_error, "simulated write failure")
  expect_true(is.na(row$primary_report))
  expect_length(list.files(root, pattern = "[.]html$", recursive = TRUE), 0L)
})

test_that("an overwrite cannot present old results while writing a new analysis", {
  fixture <- ftq_regional_study_fixture(tempfile("rf"))
  args <- ftq_public_regional_args(fixture, tempfile("rr"))
  args$write_exploratory_report <- args$write_qc_viewer <- FALSE
  prior <- do.call(fibrtraq_study, args)
  page <- prior$report_index
  expect_identical(prior$index$run_status, "completed")
  expect_match(paste(readLines(page), collapse = "\n"), "ftq-primary", fixed = TRUE)
  saved <- c(page, prior$index_csv, prior$index$result_csv)
  hashes <- tools::md5sum(saved)

  # If the canonical page is locked, the previous run must stay intact.
  testthat::with_mocked_bindings(
    expect_error(do.call(fibrtraq_study, args), "report locked"),
    ftq_write_report_page = function(...) stop("report locked"), .package = "FibrTraQ"
  )
  expect_identical(tools::md5sum(saved), hashes)

  replacement <- tempfile("rt", fileext = ".trk")
  ftq_write_regional_study_test_trk(replacement, list(rbind(c(2, 2, 2), c(4, 2, 2))))
  args$tractograms[[1L]] <- replacement
  original_run <- FibrTraQ:::ftq_run
  seen_during_run <- NULL
  testthat::local_mocked_bindings(
    ftq_run = function(...) {
      seen_during_run <<- paste(readLines(page), collapse = "\n")
      original_run(...)
    },
    ftq_write_unified_report = function(...) stop("assembly interrupted"),
    .package = "FibrTraQ"
  )
  expect_error(do.call(fibrtraq_study, args), "assembly interrupted")
  expect_match(seen_during_run, "analysis in progress", fixed = TRUE)
  expect_false(grepl("ftq-primary|href=", seen_during_run))
  failed <- paste(readLines(page), collapse = "\n")
  expect_match(failed, "report generation failed", fixed = TRUE)
  expect_match(failed, "assembly interrupted", fixed = TRUE)
  expect_false(grepl("ftq-primary|href=", failed))
  row <- utils::read.csv(prior$index_csv)
  expect_true(is.na(row$primary_report))
  expect_identical(row$run_status, "failed")
  expect_identical(row$validation_status, "passed")
  expect_equal(row$processed_streamlines, 1L)
  expect_equal(utils::read.csv(row$roi_summary_csv)$processed_streamlines, 1L)
  expect_equal(prior$index$processed_streamlines, 2L)
})

test_that("a cohort overwrite preserves excluded analyses and accepts a partial index", {
  fixture <- ftq_regional_study_fixture(tempfile("cf"))
  args <- ftq_public_regional_args(fixture, tempfile("cr"))
  args$write_report <- args$write_exploratory_report <- args$write_qc_viewer <- FALSE
  args$tractograms <- c(`001` = fixture$tractogram, `002` = fixture$tractogram)
  args$region_masks <- list(`001` = fixture$masks, `002` = fixture$masks)
  prior <- do.call(fibrtraq_study, args)
  saved <- c(prior$report_index, prior$index_csv, prior$index$result_csv)
  hashes <- tools::md5sum(saved)

  different <- args
  names(different$tractograms)[2] <- names(different$region_masks)[2] <- "003"
  expect_error(do.call(fibrtraq_study, different), "new output directory")
  different <- args
  different$focal_regions <- "B"
  expect_error(do.call(fibrtraq_study, different), "new output directory")
  expect_identical(tools::md5sum(saved), hashes)

  utils::write.csv(prior$index[1L, ], prior$index_csv, row.names = FALSE)
  resumed <- do.call(fibrtraq_study, args)
  expect_identical(resumed$index$subject, c("001", "002"))
  expect_true(all(resumed$index$run_status == "completed"))
  expect_length(list.files(args$output_dir, pattern = "[.]html$", recursive = TRUE), 1L)
})

test_that("overwriting a historical bundle preserves its reports and results", {
  fixture <- ftq_regional_study_fixture(tempfile("historical-fixture-"))
  root <- tempfile("historical-results-")
  report <- file.path(root, "001", "A", "reports", "subject", "subject_report.html")
  dir.create(dirname(report), recursive = TRUE)
  index <- file.path(root, "fibrtraq-run-index.csv")
  writeLines("historical index", index)
  writeLines("historical report", report)
  original <- tools::md5sum(c(index, report))
  args <- ftq_public_regional_args(fixture, root)
  expect_error(do.call(fibrtraq_study, args), "new output directory")
  expect_identical(tools::md5sum(c(index, report)), original)
  expect_setequal(list.files(root, recursive = TRUE, full.names = TRUE), c(index, report))
})

test_that("validation failure retains data without appearing to pass", {
  run <- list(outputs = list(result_csv = "saved.csv"),
              validation = list(passed = FALSE, biological_qc_status = "unavailable"))
  row <- FibrTraQ:::ftq_regional_study_run_row(
    "001", "002", run, NULL, NULL, "start", "validation failed",
    stage_errors = list(validation = "validation failed"),
    requested = c(primary_report = FALSE, exploratory_report = FALSE, interactive_3d_qc = FALSE)
  )
  expect_identical(row$analysis_status, "completed")
  expect_identical(row$validation_status, "failed")
  expect_identical(row$result_csv, "saved.csv")
  expect_false(row$validation_passed)
  expect_identical(row$validation_error, "validation failed")
})

test_that("subject defaults are explicit and sanitized folder collisions are rejected", {
  fixture <- ftq_regional_study_fixture(tempfile("public-subject-fixture-"))
  captured <- NULL
  testthat::local_mocked_bindings(
    fibrtraq_study = function(...) { captured <<- list(...); captured },
    .package = "FibrTraQ"
  )
  expect_message(fibrtraq_single(
    tractogram = fixture$tractogram, focal_mask = fixture$masks[["A"]],
    target_masks = fixture$masks["B"], output_dir = tempfile("subject-default-")
  ), "No subject identifier supplied", fixed = TRUE)
  expect_identical(names(captured$tractograms), "subject1")
  expect_message(fibrtraq_single(
    tractogram = fixture$tractogram, focal_mask = fixture$masks[["A"]],
    target_masks = fixture$masks["B"], output_dir = tempfile("subject-explicit-"),
    subject = "001", report_slices = c(2L, 4L)
  ), NA)
  expect_identical(names(captured$tractograms), "001")
  expect_identical(captured$report_slices, c(2L, 4L))
  check <- FibrTraQ:::ftq_regional_study_check_labels
  expect_error(check(c("S 1", "S-1"), "subject"), "distinct", fixed = TRUE)
  expect_error(check(c("S1", "s1"), "subject"), "distinct", fixed = TRUE)
  expect_error(check(c("..", "S1"), "subject"), "non-reserved", fixed = TRUE)
  expect_error(check(c("prepared", "A"), "region", "prepared"), "non-reserved", fixed = TRUE)
  for (label in c("CON", "con.txt", "PRN", "aux", "NUL", "COM1", "com9.nii",
                  "LPT1", "lpt9.csv", "S1.")) {
    expect_error(check(label, "subject"), "Windows", fixed = TRUE)
  }
  expect_silent(check(c("CON-1", "COM10", "LPT10", "S1.c"), "subject"))
  expect_silent(check(c("001", "S1c"), "subject"))
})
