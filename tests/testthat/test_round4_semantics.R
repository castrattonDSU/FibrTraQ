library(testthat)
library(FibrTraQ)

test_that("sample_points voxelization does not rasterize unsampled segment voxels", {
  fnode <- array(0, dim = c(5, 1, 1))
  fnode[1:5, 1, 1] <- 1
  target <- array(0, dim = dim(fnode))
  target[5, 1, 1] <- 1
  long_segment <- matrix(c(1, 1, 1, 5, 1, 1), ncol = 3, byrow = TRUE)

  r <- fibrtraq(
    fNode = fnode,
    localROI = array(0, dim = dim(fnode)),
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = list(long_segment),
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    focal_assignment_mode = "all_focal_voxels",
    streamline_voxelization_mode = "sample_points",
    backend = "r"
  )
  cpp <- fibrtraq(
    fNode = fnode,
    localROI = array(0, dim = dim(fnode)),
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = list(long_segment),
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    focal_assignment_mode = "all_focal_voxels",
    streamline_voxelization_mode = "sample_points",
    backend = "cpp"
  )

  expect_equal(r$tracks, c(1L, 0L, 0L, 0L, 1L))
  expect_equal(r$T, c(1L, 0L, 0L, 0L, 1L))
  expect_equal(attr(r, "fibrtraq_streamline_voxelization_mode"), "sample_points")
  ftq_expect_core_equal(cpp, r)
  ftq_expect_algorithm_metadata_equal(cpp, r)
})

test_that("exclusive ROI-level fates partition unique focal streamlines in R and C++", {
  fnode <- array(0, dim = c(4, 1, 1))
  fnode[1, 1, 1] <- 1
  local <- array(0, dim = dim(fnode))
  local[1:2, 1, 1] <- 1
  target <- array(0, dim = dim(fnode))
  target[2:3, 1, 1] <- 1
  tracks <- list(
    matrix(c(1, 1, 1, 2, 1, 1), ncol = 3, byrow = TRUE), # local and target-positive
    matrix(c(1, 1, 1, 1, 1, 1), ncol = 3, byrow = TRUE), # local and target-negative
    matrix(c(1, 1, 1, 3, 1, 1), ncol = 3, byrow = TRUE), # nonlocal and target-positive
    matrix(c(1, 1, 1, 4, 1, 1), ncol = 3, byrow = TRUE)  # nonlocal and target-negative
  )

  r <- fibrtraq(
    fNode = fnode,
    localROI = local,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = tracks,
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    focal_assignment_mode = "all_focal_voxels",
    backend = "r"
  )
  cpp <- fibrtraq(
    fNode = fnode,
    localROI = local,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = tracks,
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    focal_assignment_mode = "all_focal_voxels",
    backend = "cpp"
  )

  summary <- attr(r, "fibrtraq_roi_summary")
  expect_equal(summary$unique_focal_intersecting_streamlines, 4)
  expect_equal(summary$unique_local_target_positive_streamlines, 1)
  expect_equal(summary$unique_local_target_negative_streamlines, 1)
  expect_equal(summary$unique_nonlocal_target_positive_streamlines, 1)
  expect_equal(summary$unique_nonlocal_target_negative_streamlines, 1)
  expect_equal(summary$unique_leave_with_target_streamlines, 1)
  expect_equal(summary$unique_leave_no_target_streamlines, 1)
  expect_equal(summary$unique_any_target_positive_streamlines, 2)
  expect_equal(summary$unique_any_target_negative_streamlines, 2)
  expect_equal(summary$unique_exclusive_fate_sum, 4)
  expect_true(summary$unique_exclusive_fate_sum_matches_focal)
  ftq_expect_core_equal(cpp, r)
  ftq_expect_algorithm_metadata_equal(cpp, r)
  cpp_summary <- attr(cpp, "fibrtraq_roi_summary")
  summary$backend <- NULL
  cpp_summary$backend <- NULL
  expect_equal(cpp_summary, summary)
})

test_that("target summary tables separate unique-streamline and association denominators", {
  fnode <- array(0, dim = c(2, 1, 1))
  fnode[1:2, 1, 1] <- 1
  target <- array(0, dim = dim(fnode))
  target[2, 1, 1] <- 1
  tracks <- list(matrix(c(1, 1, 1, 2, 1, 1), ncol = 3, byrow = TRUE))

  res <- fibrtraq(
    fNode = fnode,
    localROI = array(0, dim = dim(fnode)),
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = tracks,
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    focal_assignment_mode = "all_focal_voxels",
    backend = "cpp"
  )
  tables <- ftq_target_summary_tables(
    res,
    attr(res, "fibrtraq_roi_summary"),
    target_family = "test",
    target_mask_voxels = c(T = 1),
    target_checksums = c(T = "abc")
  )

  expect_equal(tables$unique_streamlines$unique_target_positive_focal_intersecting_streamlines, 1)
  expect_equal(tables$unique_streamlines$percent_of_unique_focal_intersecting_streamlines, 100)
  expect_equal(tables$voxel_associations$focal_voxel_streamline_target_associations, 2)
  expect_equal(tables$voxel_associations$percent_of_all_focal_voxel_streamline_associations, 100)
  expect_false("percent_of_assigned_tracks" %in% names(tables$voxel_associations))
})

test_that("summary_connectivity labels target totals as associations, not unique streamlines", {
  fnode <- array(0, dim = c(2, 1, 1))
  fnode[1:2, 1, 1] <- 1
  target <- array(0, dim = dim(fnode))
  target[2, 1, 1] <- 1
  track <- matrix(c(1, 1, 1, 2, 1, 1), ncol = 3, byrow = TRUE)

  res <- fibrtraq(
    fNode = fnode,
    localROI = array(0, dim = dim(fnode)),
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = list(track),
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    focal_assignment_mode = "all_focal_voxels",
    backend = "cpp"
  )
  roi_summary <- attr(res, "fibrtraq_roi_summary")
  summary <- summary_connectivity(res, percent = FALSE, include_legacy_tracks = FALSE)

  expect_equal(roi_summary$`unique_target_streamlines:T`, 1)
  expect_equal(summary$target_associations, 2)
  expect_equal(summary$focal_voxel_streamline_target_associations, 2)
  expect_gt(summary$target_associations, roi_summary$`unique_target_streamlines:T`)
  expect_false("tracks" %in% names(summary))
  expect_equal(summary$rsc, 1)
  expect_equal(summary$total_focal_voxels, 2)

  legacy <- expect_warning(summary_connectivity(res, percent = FALSE), "deprecated")
  expect_equal(legacy$tracks, legacy$target_associations)
})

test_that("NULL localROI defaults to the focal ROI and exclusive fates still partition streamlines", {
  fnode <- array(0, dim = c(3, 1, 1))
  fnode[1, 1, 1] <- 1
  target <- array(0, dim = dim(fnode))
  target[2, 1, 1] <- 1
  tracks <- list(
    matrix(c(1, 1, 1, 1, 1, 1), ncol = 3, byrow = TRUE),
    matrix(c(1, 1, 1, 2, 1, 1), ncol = 3, byrow = TRUE),
    matrix(c(1, 1, 1, 3, 1, 1), ncol = 3, byrow = TRUE)
  )

  implicit <- fibrtraq(
    fNode = fnode,
    localROI = NULL,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = tracks,
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    backend = "cpp"
  )
  explicit <- fibrtraq(
    fNode = fnode,
    localROI = fnode,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = tracks,
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    backend = "r"
  )

  ftq_expect_core_equal(implicit, explicit)
  expect_true(attr(implicit, "fibrtraq_local_roi_defaulted_to_focal"))
  expect_false(attr(explicit, "fibrtraq_local_roi_defaulted_to_focal"))
  roi_summary <- attr(implicit, "fibrtraq_roi_summary")
  expect_equal(roi_summary$unique_local_target_positive_streamlines, 0)
  expect_equal(roi_summary$unique_local_target_negative_streamlines, 1)
  expect_equal(roi_summary$unique_nonlocal_target_positive_streamlines, 1)
  expect_equal(roi_summary$unique_nonlocal_target_negative_streamlines, 1)
  expect_equal(roi_summary$unique_exclusive_fate_sum, roi_summary$unique_focal_intersecting_streamlines)
})

test_that("unsupported streamline voxelization modes are rejected", {
  fnode <- array(1, dim = c(1, 1, 1))
  track <- matrix(c(1, 1, 1), ncol = 3)
  expect_error(
    fibrtraq(
      fNode = fnode,
      localROI = fnode,
      cNode_list = list(T = fnode),
      cNode_labels = "T",
      all_tracks = list(track),
      track_space = "voxel",
      discretization = "round",
      axis_flip = c(FALSE, FALSE, FALSE),
      streamline_voxelization_mode = "segment_supercover"
    ),
    "sample_points"
  )
})

test_that("segment interpolation detects voxels between sparse stored points", {
  fnode <- array(0, dim = c(4, 1, 1))
  fnode[2, 1, 1] <- 1
  target <- array(0, dim = dim(fnode))
  target[3, 1, 1] <- 1
  local <- array(1, dim = dim(fnode))
  track <- rbind(c(1, 1, 1), c(4, 1, 1))

  sample_points <- fibrtraq(
    fNode = fnode,
    localROI = local,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = list(track),
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    streamline_voxelization_mode = "sample_points",
    backend = "r"
  )
  expect_equal(sample_points$tracks, 0L)

  interpolated_r <- fibrtraq(
    fNode = fnode,
    localROI = local,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = list(track),
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    streamline_voxelization_mode = "segment_interpolation",
    backend = "r"
  )
  expect_equal(interpolated_r$tracks, 1L)
  expect_equal(interpolated_r$T, 1L)
  expect_equal(interpolated_r$local, 1L)

  interpolated_cpp <- fibrtraq(
    fNode = fnode,
    localROI = local,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = list(track),
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    streamline_voxelization_mode = "segment_interpolation",
    backend = "cpp"
  )
  ftq_expect_core_equal(interpolated_cpp, interpolated_r)
  ftq_expect_algorithm_metadata_equal(interpolated_cpp, interpolated_r)
})
