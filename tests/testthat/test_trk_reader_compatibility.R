test_that("TRK memory import accepts the legacy reader signature", {
  path <- ftq_example_data()$tractogram
  records <- list(list(num_points = 2L, coords = rbind(c(1, 2, 3), c(4, 5, 6))))
  expected <- list(header = list(id_string = "TRACK"), tracks = records)
  called <- FALSE
  original_reader <- freesurferformats::read.dti.trk
  reader <- ftq_test_function_copy(ftq_read_tractogram, exports = list(
    `freesurferformats::read.dti.trk` = function(filepath, shift_origin = TRUE) {
      expect_identical(filepath, path)
      expect_true(shift_origin)
      called <<- TRUE
      expected
    }
  ))
  expect_no_warning(actual <- reader(path, "memory"))
  expect_true(called)
  expect_identical(actual, expected)
  expect_identical(freesurferformats::read.dti.trk, original_reader)
})

test_that("TRK memory import requests native coordinates without hiding warnings", {
  path <- ftq_example_data()$tractogram
  records <- list(list(num_points = 2L, coords = rbind(c(1, 2, 3), c(4, 5, 6))))
  expected <- list(header = list(id_string = "TRACK"), tracks = records)
  seen_coords <- NULL
  original_reader <- freesurferformats::read.dti.trk
  reader <- ftq_test_function_copy(ftq_read_tractogram, exports = list(
    `freesurferformats::read.dti.trk` = function(filepath, shift_origin = TRUE, coords = NULL) {
      seen_coords <<- coords
      warning("independent import warning", call. = FALSE)
      expected
    }
  ))
  expect_warning(actual <- reader(path, "memory"),
                 "independent import warning")
  expect_identical(seen_coords, "native")
  expect_identical(actual, expected)
  expect_identical(freesurferformats::read.dti.trk, original_reader)
})

test_that("TRK memory import preserves stored coordinates and normalizes containers", {
  affine <- diag(c(2, 3, 4, 1))
  affine[1:3, 4L] <- c(-5, 12, 9)
  fixture <- tck_fixture_case(affine)
  root <- tempfile("trk-reader-compatibility-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  path <- file.path(root, "translated-anisotropic.trk")
  tck_fixture_write_trk(path, fixture$voxel_tracks, fixture$dims, affine)
  expected_coords <- lapply(fixture$voxel_tracks, function(points) {
    sweep(points - 0.5, 2L, c(2, 3, 4), "*")
  })
  expect_no_warning(memory <- ftq_read_tractogram(path, "memory"))
  expect_true(is.list(memory$tracks))
  expect_false(is.object(memory$tracks))
  expect_length(memory$tracks, length(expected_coords))
  expect_equal(memory$header$vox2ras, affine, tolerance = 0)
  streamed <- ftq_trk_stream_sample(ftq_trk_stream(path), length(expected_coords))
  for (i in seq_along(expected_coords)) {
    expect_equal(memory$tracks[[i]]$coords, expected_coords[[i]], tolerance = 0)
    expect_equal(memory$tracks[[i]]$coords, streamed[[i]]$coords, tolerance = 0)
    expect_identical(memory$tracks[[i]]$num_points, nrow(expected_coords[[i]]))
  }
})

test_that("direct TRK containers and list forms retain known counts in both backends", {
  affine <- diag(c(2, 3, 4, 1))
  affine[1:3, 4L] <- c(-5, 12, 9)
  fixture <- tck_fixture_case(affine)
  root <- tempfile("trk-container-counts-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  path <- file.path(root, "translated-anisotropic.trk")
  tck_fixture_write_trk(path, fixture$voxel_tracks, fixture$dims, affine)
  read_trk <- freesurferformats::read.dti.trk
  read_args <- list(filepath = path)
  if ("coords" %in% names(formals(read_trk))) read_args$coords <- "native"
  expect_no_warning(direct <- do.call(read_trk, read_args))
  expect_null(ftq_track_format(direct$tracks))
  expect_identical(ftq_track_count(direct$tracks), length(fixture$voxel_tracks))
  records <- as.list(direct$tracks)
  representations <- list(
    dependency_object = direct,
    dependency_tracks = direct$tracks,
    memory_import = ftq_read_tractogram(path, "memory"),
    list_records = records,
    coordinate_matrices = lapply(records, function(track) track$coords)
  )
  results <- list()
  for (backend in c("r", "cpp")) {
    for (name in names(representations)) {
      result <- FibrTraQ(
        fNode = fixture$focal, cNode_list = fixture$targets,
        cNode_labels = names(fixture$targets), all_tracks = representations[[name]],
        ref_image = fixture$reference, track_space = "voxelmm",
        track_voxel_size = c(2, 3, 4), axis_flip = rep(FALSE, 3L),
        discretization = "round", assignment_mode = "full_streamline",
        streamline_voxelization_mode = "sample_points",
        endpoint_analysis = "exact_and_dilated1", backend = backend,
        verbose = FALSE, progress = FALSE
      )
      expect_equal(ftq_drop_result_attrs(result), tck_fixture_expected(), tolerance = 0,
                   info = paste(backend, name))
      results[[paste(backend, name)]] <- result
    }
  }
  baseline <- results[[1L]]
  for (result in results) {
    for (attribute in c("fibrtraq_endpoint_summary", "fibrtraq_endpoint_penetration",
                        "fibrtraq_endpoint_context")) {
      expect_equal(attr(result, attribute), attr(baseline, attribute), tolerance = 0)
    }
    actual_roi <- attr(result, "fibrtraq_roi_summary")
    expected_roi <- attr(baseline, "fibrtraq_roi_summary")
    fields <- grep("^(total_|processed_|unique_)", names(expected_roi), value = TRUE)
    expect_equal(actual_roi[fields], expected_roi[fields], tolerance = 0)
  }
})
