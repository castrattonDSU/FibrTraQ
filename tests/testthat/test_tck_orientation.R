tck_orientation_cases <- function() {
  translated <- diag(c(1.25, 2.5, 3.75, 1))
  translated[1:3, 4L] <- c(-19, 11, 5)
  permuted <- rbind(c(0, 0, 3, 9), c(1.25, 0, 0, -4),
                    c(0, 2.5, 0, 12), c(0, 0, 0, 1))
  reflected <- diag(c(-1.5, 2.5, 3, 1))
  reflected[1:3, 4L] <- c(30, -11, 7)
  reflected_two <- diag(c(-2, -3, 1.5, 1))
  reflected_two[1:3, 4L] <- c(33, 28, -14)
  angle <- pi / 9
  oblique <- diag(4L)
  oblique[1:3, 1:3] <- rbind(c(cos(angle), -sin(angle), 0),
                            c(sin(angle), cos(angle), 0), c(0, 0, 1)) %*%
    diag(c(1.25, 2.5, 3.75))
  oblique[1:3, 4L] <- c(-13.5, 8.25, 5.75)
  sheared <- rbind(c(1.5, 0.3, 0, -5), c(0, 2, 0.2, 8),
                   c(0, 0, 3, -7), c(0, 0, 0, 1))
  qform <- tck_fixture_case(oblique)
  qform_only <- function(image) {
    header <- RNifti::niftiHeader(image)
    header$qform_code <- 2L
    header$sform_code <- 0L
    # A positive rotation around z has quaternion (0, 0, sin(angle / 2)).
    header$quatern_b <- header$quatern_c <- 0
    header$quatern_d <- sin(angle / 2)
    header$qoffset_x <- oblique[1L, 4L]
    header$qoffset_y <- oblique[2L, 4L]
    header$qoffset_z <- oblique[3L, 4L]
    RNifti::asNifti(array(as.numeric(image), dim(image)), reference = header)
  }
  qform$focal <- qform_only(qform$focal)
  qform$targets <- lapply(qform$targets, qform_only)
  qform$reference <- qform_only(qform$reference)
  list(
    identity = tck_fixture_case(),
    translated_anisotropic = tck_fixture_case(translated),
    permuted_rotated = tck_fixture_case(permuted, voxel_order = "ASR"),
    negative_determinant = tck_fixture_case(reflected, voxel_order = "LAS"),
    two_negative_axes = tck_fixture_case(reflected_two, voxel_order = "LPS"),
    oblique = tck_fixture_case(oblique),
    sheared = tck_fixture_case(sheared),
    qform_only = qform,
    metres = tck_fixture_case(oblique, units = "m"),
    micrometres = tck_fixture_case(oblique, units = "um")
  )
}

tck_orientation_args <- function(fixture, tracks, backend = "r",
                                  assignment_mode = "full_streamline",
                                  voxelization = "sample_points",
                                  reverse = FALSE) {
  list(
    fNode = fixture$focal, cNode_list = fixture$targets,
    cNode_labels = names(fixture$targets), all_tracks = tracks,
    ref_image = fixture$reference, discretization = "round",
    assignment_mode = assignment_mode, search_radius_mm = 0.2,
    streamline_voxelization_mode = voxelization,
    reverse_streamlines = reverse, backend = backend,
    endpoint_analysis = "exact_and_dilated1", verbose = FALSE, progress = FALSE
  )
}

tck_expect_science_equal <- function(actual, expected, info = NULL) {
  expect_equal(ftq_drop_result_attrs(actual), ftq_drop_result_attrs(expected),
               tolerance = 0, info = info)
  a_roi <- attr(actual, "fibrtraq_roi_summary")
  e_roi <- attr(expected, "fibrtraq_roi_summary")
  count_fields <- grep("^(total_|processed_|unique_)", names(e_roi), value = TRUE)
  expect_equal(a_roi[count_fields], e_roi[count_fields], tolerance = 0, info = info)
  for (name in c("fibrtraq_endpoint_summary", "fibrtraq_endpoint_penetration",
                 "fibrtraq_endpoint_context")) {
    expect_equal(attr(actual, name), attr(expected, name),
                 tolerance = 1e-5, info = paste(info, name))
  }
  expect_equal(ftq_target_summary_tables(actual, a_roi),
               ftq_target_summary_tables(expected, e_roi),
               tolerance = 0, info = info)
}

for (tck_orientation_backend in c("r", "cpp")) {
  local({
    backend <- tck_orientation_backend
    test_that(paste("TCK and TRK recover independent affine known answers with", backend), {
      if (identical(backend, "cpp")) {
        skip_if_not(FibrTraQ:::ftq_compiled_backend_available())
      }
      root <- tempfile("tck-orientation-")
      dir.create(root)
      on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
      cases <- tck_orientation_cases()
      for (name in names(cases)) {
        fixture <- cases[[name]]
        paths <- tck_fixture_files(file.path(root, name), fixture)
        before <- tools::md5sum(unlist(paths, use.names = FALSE))
        fixture$focal <- RNifti::readNifti(paths$focal)
        fixture$targets <- lapply(paths$targets, RNifti::readNifti)
        fixture$reference <- RNifti::readNifti(paths$reference)
        for (reader in c("stream", "memory")) {
          tck <- FibrTraQ:::ftq_read_tractogram(paths$tck, reader, chunk_size = 2L)
          args <- tck_orientation_args(fixture, tck, backend)
          actual <- do.call(FibrTraQ, args)
          expect_equal(ftq_drop_result_attrs(actual), tck_fixture_expected(),
                       tolerance = 0, info = paste(name, reader))
          expect_identical(attr(actual, "fibrtraq_track_space"), "world")
          expect_identical(attr(actual, "fibrtraq_axis_flip"), rep(FALSE, 3L))
          expect_identical(attr(actual, "fibrtraq_mapping_selection_rule"),
                           "tck_world_ras_mm_from_nifti_affine")
          # RSC is coverage of the three focal voxels; this assertion does not
          # use either format's spatial converter to construct the oracle.
          expect_equal(colMeans(actual[c("target_a", "target_b")] > 0) * 100,
                       c(target_a = 100, target_b = 100 / 3))

          args$all_tracks <- FibrTraQ:::ftq_read_tractogram(paths$trk, reader,
                                                           chunk_size = 2L)
          args$track_space <- "trackvis"
          args$axis_flip <- rep(FALSE, 3L)
          trk_result <- do.call(FibrTraQ, args)
          tck_expect_science_equal(actual, trk_result, paste(name, reader))
        }
        expect_identical(tools::md5sum(names(before)), before)
      }
    })

    test_that(paste("TCK parity covers assignment modes, dense paths and reversal with", backend), {
      if (identical(backend, "cpp")) {
        skip_if_not(FibrTraQ:::ftq_compiled_backend_available())
      }
      root <- tempfile("tck-methods-")
      dir.create(root)
      on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
      fixture <- tck_orientation_cases()$oblique
      paths <- tck_fixture_files(root, fixture, datatype = "Float64BE")
      for (assignment_mode in c("full_streamline", "endpoint_voxels", "local_search")) {
        for (voxelization in c("sample_points", "segment_interpolation")) {
          # The baseline is explicit one-based voxel coordinates, independently
          # authored before serialization to either tractogram format.
          baseline_args <- tck_orientation_args(
            fixture, fixture$voxel_tracks, backend, assignment_mode, voxelization
          )
          baseline_args$track_space <- "voxel"
          baseline_args$axis_flip <- rep(FALSE, 3L)
          baseline <- do.call(FibrTraQ, baseline_args)
          expect_equal(ftq_drop_result_attrs(baseline),
                       tck_fixture_expected(assignment_mode, voxelization),
                       tolerance = 0, info = paste(assignment_mode, voxelization))
          for (reader in c("stream", "memory")) {
            trk_args <- baseline_args
            trk_args$all_tracks <- FibrTraQ:::ftq_read_tractogram(paths$trk, reader,
                                                                chunk_size = 1L)
            trk_args$track_space <- "trackvis"
            trk_result <- do.call(FibrTraQ, trk_args)
            tck_expect_science_equal(trk_result, baseline,
                                     paste(assignment_mode, voxelization, reader, "TRK"))
            for (reverse in c(FALSE, TRUE)) {
              args <- tck_orientation_args(
                fixture,
                FibrTraQ:::ftq_read_tractogram(paths$tck, reader, chunk_size = 1L),
                backend, assignment_mode, voxelization, reverse
              )
              actual <- do.call(FibrTraQ, args)
              tck_expect_science_equal(actual, baseline,
                                       paste(assignment_mode, voxelization, reader, reverse))
              if (identical(assignment_mode, "local_search")) {
                assignments <- attr(actual, "fibrtraq_local_search_assignments")
                expected <- attr(baseline, "fibrtraq_local_search_assignments")
                # Start/end columns exchange on reversal. Pair identity and
                # distances must retain the same scientific interpretation.
                fields <- c("streamline_index", "paired_target",
                            "paired_focal_distance_mm", "paired_target_distance_mm",
                            "paired_max_distance_mm", "paired_target_intersects_full_path")
                expect_equal(assignments[fields], expected[fields], tolerance = 1e-5)
              }
            }
          }
        }
      }
    })
  })
}

test_that("TCK uses its world convention even when reflections score more overlap", {
  fixture <- tck_fixture_case()
  true <- fixture$voxel_tracks[[1L]]
  reflected <- sweep(-true, 2L, fixture$dims + 1L, "+")
  voxel <- c(list(true), rep(list(reflected), 6L))
  world <- lapply(voxel, function(points) points - 1)
  path <- tempfile(fileext = ".tck")
  on.exit(unlink(path), add = TRUE)
  tck_fixture_write_tck(path, world)
  args <- tck_orientation_args(fixture,
    FibrTraQ:::ftq_read_tractogram(path, "stream", chunk_size = 2L))
  result <- do.call(FibrTraQ, args)
  expect_equal(result$tracks, rep(1, 3L))
  expect_equal(result$target_a, rep(1, 3L))
  expect_equal(attr(result, "fibrtraq_roi_summary")$processed_streamlines, 7)
  expect_identical(attr(result, "fibrtraq_axis_flip"), rep(FALSE, 3L))
})

test_that("TCK centres and TrackVis corners preserve independent boundary decisions", {
  root <- tempfile("tck-boundary-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  affine <- diag(c(2, 4, 8, 1))
  affine[1:3, 4L] <- c(-8, -12, -24)
  fixture <- tck_fixture_case(affine)
  # Exact half-voxel ties use round-to-even. Nearby values are separated by
  # 0.01 voxel, well beyond Float32 quantization at these world magnitudes.
  # Arbitrary oblique points within rounding error of a boundary do not imply
  # bitwise equality across independently quantized file formats.
  x <- c(1.49, 1.5, 1.51, 2.49, 2.5, 2.51,
         3.49, 3.5, 3.51, 4.49, 4.5, 4.51)
  fixture$voxel_tracks <- lapply(x, function(value) {
    rbind(c(value, 3, 2), c(8, 3, 2))
  })
  fixture$world_tracks <- lapply(fixture$voxel_tracks, function(voxel) {
    (cbind(voxel - 1, 1) %*% t(affine))[, 1:3, drop = FALSE]
  })
  paths <- tck_fixture_files(root, fixture)
  expected_counts <- list(round = c(4, 2, 4), floor = c(3, 3, 3),
                          ceil = c(3, 3, 3), any = c(5, 2, 5))
  backends <- "r"
  if (FibrTraQ:::ftq_compiled_backend_available()) backends <- c(backends, "cpp")
  for (backend in backends) {
    for (method in names(expected_counts)) {
      expected <- tck_fixture_expected()
      expected$tracks <- expected$target_a <- expected_counts[[method]]
      expected$local <- expected$leave_no_target <- expected$target_b <- rep(0, 3L)
      for (reader in c("stream", "memory")) {
        args <- tck_orientation_args(fixture,
          FibrTraQ:::ftq_read_tractogram(paths$tck, reader, chunk_size = 3L), backend)
        args$discretization <- method
        actual <- do.call(FibrTraQ, args)
        expect_equal(ftq_drop_result_attrs(actual), expected, tolerance = 0,
                     info = paste(backend, method, reader, "TCK"))
        args$all_tracks <- FibrTraQ:::ftq_read_tractogram(paths$trk, reader,
                                                         chunk_size = 3L)
        args$track_space <- "trackvis"
        args$axis_flip <- rep(FALSE, 3L)
        trk <- do.call(FibrTraQ, args)
        expect_equal(ftq_drop_result_attrs(trk), expected, tolerance = 0,
                     info = paste(backend, method, reader, "TRK"))
      }
    }
  }
})

test_that("TCK local search measures subvoxel distances in physical millimetres", {
  root <- tempfile("tck-physical-radius-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  affine <- diag(c(2, 4, 8, 1))
  affine[1:3, 4L] <- c(-8, -12, -24)
  backends <- "r"
  if (FibrTraQ:::ftq_compiled_backend_available()) backends <- c(backends, "cpp")
  for (units in c("mm", "m", "um")) {
    fixture <- tck_fixture_case(affine, units)
    # Each endpoint is 0.25 voxel from its nearest labelled centre along x,
    # whose spacing is 2 mm: the physical search distance must be 0.5 mm.
    fixture$voxel_tracks <- list(rbind(c(2.25, 3, 2), c(7.75, 3, 2)))
    fixture$world_tracks <- lapply(fixture$voxel_tracks, function(voxel) {
      (cbind(voxel - 1, 1) %*% t(affine))[, 1:3, drop = FALSE]
    })
    paths <- tck_fixture_files(file.path(root, units), fixture)
    for (backend in backends) {
      for (radius in c(0.49, 0.51)) {
        accepted <- radius > 0.5
        expected <- data.frame(
          fNode_Coords = c("2_3_2", "3_3_2", "4_3_2"),
          tracks = c(1, 0, 0), local = c(0, 0, 0),
          leave_no_target = c(as.integer(!accepted), 0, 0),
          target_a = c(as.integer(accepted), 0, 0), target_b = c(0, 0, 0),
          none = c(0, 1, 1), stringsAsFactors = FALSE
        )
        for (reader in c("stream", "memory")) {
          args <- tck_orientation_args(fixture,
            FibrTraQ:::ftq_read_tractogram(paths$tck, reader), backend,
            assignment_mode = "local_search")
          args$search_radius_mm <- radius
          actual <- do.call(FibrTraQ, args)
          expect_equal(ftq_drop_result_attrs(actual), expected, tolerance = 0,
                       info = paste(units, backend, radius, reader))
          if (accepted) {
            assignments <- attr(actual, "fibrtraq_local_search_assignments")
            expect_equal(assignments$paired_focal_distance_mm, 0.5, tolerance = 1e-5)
            expect_equal(assignments$paired_target_distance_mm, 0.5, tolerance = 1e-5)
          }
          args$all_tracks <- FibrTraQ:::ftq_read_tractogram(paths$trk, reader)
          args$track_space <- "trackvis"
          args$axis_flip <- rep(FALSE, 3L)
          trk <- do.call(FibrTraQ, args)
          tck_expect_science_equal(actual, trk, paste(units, backend, radius, reader))
        }
      }
    }
  }
})

test_that("TCK refuses explicit alternative coordinate interpretations", {
  fixture <- tck_fixture_case()
  path <- tempfile(fileext = ".tck")
  on.exit(unlink(path), add = TRUE)
  tck_fixture_write_tck(path, fixture$world_tracks)
  args <- tck_orientation_args(fixture, FibrTraQ:::ftq_read_tractogram(path, "memory"))
  for (space in c("voxel", "voxelmm", "vox10", "trackvis")) {
    bad <- args
    bad$track_space <- space
    expect_error(do.call(FibrTraQ, bad), "TCK|tck")
  }
  for (axis in seq_len(3L)) {
    bad <- args
    bad$axis_flip <- seq_len(3L) == axis
    expect_error(do.call(FibrTraQ, bad), "TCK|tck")
  }
  bad <- args
  bad$track_affine <- diag(4L)
  expect_error(do.call(FibrTraQ, bad), "TCK|tck")
  bad <- args
  bad$track_voxel_size <- c(1, 1, 1)
  expect_error(do.call(FibrTraQ, bad), "TCK|tck")
  args$track_space <- "world"
  args$axis_flip <- rep(FALSE, 3L)
  expect_equal(ftq_drop_result_attrs(do.call(FibrTraQ, args)), tck_fixture_expected())
})

test_that("TCK requires declared units and unambiguous invertible NIfTI geometry", {
  fixture <- tck_fixture_case()
  settings <- function(reference) {
    FibrTraQ:::ftq_tck_coordinate_settings("auto", NULL, reference)
  }
  plain_focal <- array(as.numeric(fixture$focal), fixture$dims)
  expect_error(settings(plain_focal), "NIfTI|affine|reference")
  undefined <- RNifti::asNifti(array(0, fixture$dims))
  RNifti::pixunits(undefined) <- "mm"
  expect_error(settings(undefined), "affine|qform|sform|transform")
  unknown_units <- tck_fixture_image(array(0, fixture$dims), units = "Unknown")
  expect_error(settings(unknown_units), "unit")

  singular_affine <- diag(4L)
  singular_affine[2L, ] <- singular_affine[1L, ]
  singular <- tck_fixture_image(array(0, fixture$dims), singular_affine)
  expect_error(settings(singular), "affine|singular|invert|transform")

  header <- RNifti::niftiHeader(fixture$reference)
  header$qform_code <- 2L
  header$qoffset_x <- 10
  header$pixdim <- rep(1, 8L)
  conflict <- RNifti::asNifti(array(0, fixture$dims), reference = header)
  expect_gt(RNifti::niftiHeader(conflict)$qform_code, 0)
  expect_gt(RNifti::niftiHeader(conflict)$sform_code, 0)
  expect_error(settings(conflict), "qform|sform|conflict|disagree")

  header$qoffset_x <- 0
  consistent <- RNifti::asNifti(array(0, fixture$dims), reference = header)
  expect_gt(RNifti::niftiHeader(consistent)$qform_code, 0)
  resolved <- settings(consistent)
  expect_identical(resolved$track_space, "world")
  expect_identical(resolved$axis_flip, rep(FALSE, 3L))
})

test_that("TCK uses a declared focal grid when no separate reference is supplied", {
  fixture <- tck_orientation_cases()$negative_determinant
  path <- tempfile(fileext = ".tck")
  on.exit(unlink(path), add = TRUE)
  tck_fixture_write_tck(path, fixture$world_tracks)
  args <- tck_orientation_args(fixture, FibrTraQ:::ftq_read_tractogram(path, "stream"))
  args$ref_image <- NULL
  expect_equal(ftq_drop_result_attrs(do.call(FibrTraQ, args)), tck_fixture_expected())
  args$fNode <- array(as.numeric(fixture$focal), fixture$dims)
  args$cNode_list <- lapply(fixture$targets, function(mask) {
    array(as.numeric(mask), fixture$dims)
  })
  expect_error(do.call(FibrTraQ, args), "NIfTI|affine|reference")
  args$ref_image <- fixture$reference
  expect_equal(ftq_drop_result_attrs(do.call(FibrTraQ, args)), tck_fixture_expected())
})

test_that("TCK validates every NIfTI ROI against the reference in mixed inputs", {
  fixture <- tck_fixture_case()
  path <- tempfile(fileext = ".tck")
  on.exit(unlink(path), add = TRUE)
  tck_fixture_write_tck(path, fixture$world_tracks)
  plain <- function(image) array(as.numeric(image), fixture$dims)
  for (reader in c("stream", "memory")) {
    args <- tck_orientation_args(fixture, FibrTraQ:::ftq_read_tractogram(path, reader))
    args$fNode <- plain(fixture$focal)
    args$cNode_list <- list(
      target_a = tck_fixture_image(plain(fixture$targets$target_a), units = "m"),
      target_b = plain(fixture$targets$target_b)
    )
    args$localROI <- tck_fixture_image(plain(fixture$focal), units = "um")
    args$endpoint_context_list <- list(anatomy = fixture$focal)
    expect_equal(ftq_drop_result_attrs(do.call(FibrTraQ, args)), tck_fixture_expected())

    shifted <- diag(4L)
    shifted[1L, 4L] <- 100
    bad_images <- list(
      shifted_affine = tck_fixture_image(plain(fixture$focal), shifted),
      different_dimensions = tck_fixture_image(array(1, fixture$dims + c(1L, 0L, 0L)))
    )
    for (problem in names(bad_images)) {
      for (role in c("target", "local", "context")) {
        bad <- args
        if (identical(role, "target")) bad$cNode_list$target_a <- bad_images[[problem]]
        if (identical(role, "local")) bad$localROI <- bad_images[[problem]]
        if (identical(role, "context")) {
          bad$endpoint_context_list$anatomy <- bad_images[[problem]]
        }
        expect_error(do.call(FibrTraQ, bad), "TCK ROI images.*same dimensions.*affine/grid",
                     info = paste(reader, problem, role))
      }
    }

    # Plain masks and coordinate matrices have no affine of their own. They
    # continue to use the declared reference grid and existing bounds checks.
    args$cNode_list <- list(target_a = matrix(c(8, 3, 2), ncol = 3L),
                           target_b = matrix(c(3, 7, 5), ncol = 3L))
    args$localROI <- rbind(c(2, 3, 2), c(3, 3, 2), c(4, 3, 2))
    args$endpoint_context_list$anatomy <- plain(fixture$focal)
    expect_equal(ftq_drop_result_attrs(do.call(FibrTraQ, args)), tck_fixture_expected())
    args$cNode_list$target_a <- matrix(c(13, 3, 2), ncol = 3L)
    expect_error(do.call(FibrTraQ, args), "bounds|outside|grid")
  }
})
