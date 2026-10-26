ftq_orientation_guard_fixture <- function() {
  focal <- array(0, dim = c(12L, 10L, 8L))
  target <- focal
  focal[2, 2, 2] <- 1
  target[5, 2, 2] <- 1

  # The known anatomical convention is the unreflected voxel grid. One track
  # intersects the focal and target masks. Three unrelated tracks occupy the
  # opposite corner; reflecting all axes makes those three look like hits.
  true_track <- rbind(c(2, 2, 2), c(5, 2, 2))
  reflected_track <- rbind(c(11, 9, 7), c(8, 9, 7))
  list(
    focal = focal,
    target = target,
    # TrackVis voxel-mm coordinates place voxel centres half a voxel from the
    # origin. Unit spacing makes this explicit and independent of the mapper.
    tracks = lapply(list(true_track, reflected_track, reflected_track,
                        reflected_track), function(x) x - 0.5)
  )
}

ftq_orientation_guard_arguments <- function(fixture, backend = "r") {
  list(
    fNode = fixture$focal,
    cNode_list = list(target = fixture$target),
    cNode_labels = "target",
    all_tracks = fixture$tracks,
    track_voxel_size = c(1, 1, 1),
    track_space = "auto",
    discretization = "round",
    assignment_mode = "full_streamline",
    backend = backend,
    verbose = FALSE,
    progress = FALSE
  )
}

ftq_orientation_guard_write_trk <- function(path, tracks, dims) {
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  write_fixed <- function(text, n) {
    bytes <- raw(n)
    value <- charToRaw(text)
    if (length(value)) bytes[seq_along(value)] <- value
    writeBin(bytes, con)
  }
  write_fixed("TRACK", 6L)
  writeBin(as.integer(dims), con, size = 2L, endian = "little")
  writeBin(c(1, 1, 1), con, size = 4L, endian = "little")
  writeBin(c(0, 0, 0), con, size = 4L, endian = "little")
  writeBin(0L, con, size = 2L, endian = "little")
  write_fixed("", 200L)
  writeBin(0L, con, size = 2L, endian = "little")
  write_fixed("", 200L)
  writeBin(as.numeric(t(diag(4))), con, size = 4L, endian = "little")
  write_fixed("", 444L)
  write_fixed("RAS", 4L)
  write_fixed("", 4L)
  writeBin(c(1, 0, 0, 0, 1, 0), con, size = 4L, endian = "little")
  write_fixed("", 2L)
  writeBin(rep(0L, 6L), con, size = 1L, endian = "little")
  writeBin(as.integer(length(tracks)), con, size = 4L, endian = "little")
  writeBin(c(2L, 1000L), con, size = 4L, endian = "little")
  for (track in tracks) {
    writeBin(as.integer(nrow(track)), con, size = 4L, endian = "little")
    writeBin(as.numeric(t(track)), con, size = 4L, endian = "little")
  }
  invisible(path)
}

for (orientation_backend in c("r", "cpp")) {
  local({
    backend <- orientation_backend
    test_that(paste("automatic reflections stop before counting with", backend), {
      if (identical(backend, "cpp")) {
        skip_if_not(FibrTraQ:::ftq_compiled_backend_available())
      }
      fixture <- ftq_orientation_guard_fixture()
      args <- ftq_orientation_guard_arguments(fixture, backend)
      for (entrypoint in list(FibrTraQ, fibrtraq)) {
        expect_error(do.call(entrypoint, args), "Automatic mapping would reflect")
        # Pinning only the units must not implicitly authorize reflection.
        args$track_space <- "voxelmm"
        expect_error(do.call(entrypoint, args), "Automatic mapping would reflect")
        args$track_space <- "auto"
      }
    })

    test_that(paste("explicit orientation preserves known counts with", backend), {
      if (identical(backend, "cpp")) {
        skip_if_not(FibrTraQ:::ftq_compiled_backend_available())
      }
      args <- ftq_orientation_guard_arguments(
        ftq_orientation_guard_fixture(), backend
      )
      for (entrypoint in list(FibrTraQ, fibrtraq)) {
        for (flip in c(FALSE, TRUE)) {
          args$axis_flip <- rep(flip, 3L)
          result <- do.call(entrypoint, args)
          expected_hits <- if (flip) 3L else 1L
          expect_identical(result$fNode_Coords, "2_2_2")
          expect_equal(result$tracks, expected_hits)
          expect_equal(result$target, expected_hits)
          expect_equal(result$local, 0L)
          expect_equal(result$leave_no_target, 0L)
          expect_identical(attr(result, "fibrtraq_axis_flip"), rep(flip, 3L))
          expect_identical(attr(result, "fibrtraq_track_space"), "voxelmm")
          roi <- attr(result, "fibrtraq_roi_summary")
          expect_equal(roi$processed_streamlines, 4)
          expect_equal(roi$unique_focal_intersecting_streamlines, expected_hits)
        }
      }
    })

    test_that(paste("automatic mapping retains an unreflected winner with", backend), {
      if (identical(backend, "cpp")) {
        skip_if_not(FibrTraQ:::ftq_compiled_backend_available())
      }
      fixture <- ftq_orientation_guard_fixture()
      fixture$tracks <- fixture$tracks[c(1L, 1L, 1L, 2L)]
      args <- ftq_orientation_guard_arguments(fixture, backend)
      for (entrypoint in list(FibrTraQ, fibrtraq)) {
        result <- do.call(entrypoint, args)
        expect_identical(attr(result, "fibrtraq_axis_flip"), rep(FALSE, 3L))
        expect_identical(attr(result, "fibrtraq_track_space"), "voxelmm")
        expect_identical(result$fNode_Coords, "2_2_2")
        expect_equal(result$tracks, 3L)
        expect_equal(result$target, 3L)
      }
    })
  })
}

test_that("file workflow blocks an unapproved reflection without scientific outputs", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("freesurferformats")
  root <- tempfile("orientation-guard-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  fixture <- ftq_orientation_guard_fixture()
  focal_path <- file.path(root, "focal.nii.gz")
  target_path <- file.path(root, "target.nii.gz")
  trk_path <- file.path(root, "tracks.trk")
  RNifti::writeNifti(fixture$focal, focal_path)
  RNifti::writeNifti(fixture$target, target_path, template = focal_path)
  ftq_orientation_guard_write_trk(trk_path, fixture$tracks, dim(fixture$focal))
  input_checksums <- tools::md5sum(c(focal_path, target_path, trk_path))

  rejected_output <- file.path(root, "automatic")
  args <- list(
    tractogram = trk_path,
    focal_mask = focal_path,
    target_masks = c(target = target_path),
    output_dir = rejected_output,
    subject = "synthetic",
    focal_label = "focal",
    assignment_mode = "full_streamline",
    backend = "r",
    verbose = FALSE,
    progress = FALSE
  )
  expect_error(do.call(ftq_run_files, args), "Automatic mapping would reflect")
  # A failed status/index is useful, but no count, map, or scientific report
  # may make the rejected orientation look like a completed analysis.
  expect_length(list.files(file.path(rejected_output, "synthetic", "focal"),
                           recursive = TRUE, all.files = TRUE), 0L)
  failed_index <- utils::read.csv(file.path(rejected_output, "fibrtraq-run-index.csv"))
  expect_identical(failed_index$run_status, "failed")
  expect_match(failed_index$error, "Automatic mapping would reflect")

  args$output_dir <- file.path(root, "verified-orientation")
  args$track_space <- "voxelmm"
  args$axis_flip <- c(FALSE, FALSE, FALSE)
  result <- do.call(ftq_run_files, args)
  run <- result$index[1L, , drop = FALSE]
  expect_identical(run$run_status, "completed")
  expect_true(run$validation_passed)
  counts <- utils::read.csv(run$result_csv, check.names = FALSE)
  expect_identical(counts$fNode_Coords, "2_2_2")
  expect_equal(counts$tracks, 1L)
  expect_equal(counts$target, 1L)
  roi <- utils::read.csv(run$roi_summary_csv)
  expect_equal(roi$processed_streamlines, 4)
  expect_equal(roi$unique_focal_intersecting_streamlines, 1)
  expect_identical(roi$selected_axis_flip, "FALSE,FALSE,FALSE")
  metadata <- jsonlite::read_json(run$metadata_json, simplifyVector = TRUE)
  tracks_map <- as.array(RNifti::readNifti(metadata$files$map_paths$tracks))
  expect_equal(tracks_map[2, 2, 2], 1)
  expect_equal(sum(tracks_map), 1)
  expect_true(all(file.exists(c(run$primary_report, run$exploratory_report,
                               run$interactive_3d_qc))))
  expect_identical(tools::md5sum(names(input_checksums)), input_checksums)
})
