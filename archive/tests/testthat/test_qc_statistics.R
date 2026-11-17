test_that("statistical helpers summarize, compare, and generate shifted null ROIs", {
  skip_if_not_installed("RNifti")

  result <- data.frame(
    fNode_Coords = paste(seq_len(6), 1, 1, sep = "_"),
    tracks = c(10, 12, 14, 16, 18, 20),
    local = c(1, 1, 1, 1, 1, 1),
    leave_no_target = c(1, 1, 1, 1, 1, 1),
    A = c(1, 2, 4, 8, 12, 16),
    B = c(8, 8, 7, 5, 4, 2),
    none = 0,
    stringsAsFactors = FALSE
  )
  summary <- ftq_stat_subject_summary(result, subject = "S1", group = "control", roi = "left-thalamus")
  expect_equal(summary$roi$focal_voxel_streamline_associations, sum(result$tracks))
  expect_equal(summary$targets$target[1], "A")
  expect_equal(nrow(summary$voxels), 6L)

  gradients <- ftq_stat_spatial_gradients(result, min_connected_voxels = 2L)
  expect_true(all(c("target", "axis", "slope", "q_value") %in% names(gradients)))
  expect_true(any(gradients$target == "A" & gradients$axis == "x"))

  correct <- data.frame(unit = c("u1", "u2", "u3"), assigned_tracks = c(10, 12, 14), target_hit_events = c(9, 10, 13))
  control <- data.frame(unit = c("u1", "u2", "u3"), assigned_tracks = c(10, 12, 14), target_hit_events = c(1, 2, 3))
  neg <- ftq_stat_negative_control(correct, control, by = "unit")
  expect_gt(neg$median_difference, 0)
  expect_equal(nrow(neg$paired), 3L)

  group_data <- data.frame(subject = paste0("S", 1:6), group = rep(c("A", "B"), each = 3), metric = c(3, 4, 5, 1, 2, 2))
  cmp <- ftq_stat_compare_groups(group_data, metric = "metric")
  expect_equal(cmp$metric, "metric")
  expect_true(is.finite(cmp$p_value))

  tmp <- file.path(tempdir(), "shifted-rois")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  mask <- file.path(tmp, "mask.nii.gz")
  arr <- array(0, dim = c(4, 4, 4))
  arr[2, 2, 2] <- 1
  RNifti::writeNifti(arr, mask)
  shifted <- ftq_null_shifted_rois(mask, tmp, data.frame(x = 1L, y = 0L, z = 0L), overwrite = TRUE)
  expect_true(file.exists(shifted$path))
  shifted_arr <- RNifti::readNifti(shifted$path, internal = FALSE)
  expect_equal(as.array(shifted_arr)[3, 2, 2], 1)
})
