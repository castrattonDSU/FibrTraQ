library(testthat)
library(FibrTraQ)

test_that("iso_coords extracts coordinates", {
  arr <- array(0, dim = c(2,2,2))
  arr[1,1,1] <- 1
  arr[2,2,2] <- 1
  coords <- iso_coords(arr)
  expect_equal(nrow(coords), 2)
  expect_true(all(coords[1,] == c(1,1,1)))
  expect_true(all(coords[2,] == c(2,2,2)))
})

test_that("fibrtraq counts target, local, and leave-no-target tracks", {
  fnode <- array(0, dim = c(3, 1, 1))
  fnode[1, 1, 1] <- 1
  target <- array(0, dim = c(3, 1, 1))
  target[2, 1, 1] <- 1

  tracks <- list(
    matrix(c(1,1,1, 1,1,1), ncol = 3, byrow = TRUE), # local-only
    matrix(c(1,1,1, 2,1,1), ncol = 3, byrow = TRUE), # target hit
    matrix(c(1,1,1, 3,1,1), ncol = 3, byrow = TRUE)  # leaves fNode, no target
  )

  res <- fibrtraq(
    fNode = fnode,
    localROI = fnode,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = tracks,
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    focal_assignment_mode = "first_intersection"
  )

  expect_equal(res$tracks, 3L)
  expect_equal(res$local, 1L)
  expect_equal(res$leave_no_target, 1L)
  expect_equal(res$T, 1L)
  expect_equal(res$none, 0L)

  expect_warning(summary <- rsc_summary(res, percent = TRUE), "deprecated")
  expect_equal(summary$node, "T")
  expect_equal(summary$tracks, 1)
  expect_equal(summary$target_associations, 1)
})

test_that("fibrtraq synthetic benchmark gives exact known-answer counts", {
  fnode <- array(0, dim = c(6, 2, 1))
  fnode[1, 1, 1] <- 1
  fnode[2, 1, 1] <- 1
  fnode[3, 1, 1] <- 1

  local <- array(0, dim = dim(fnode))
  local[1, 1, 1] <- 1

  target_a <- array(0, dim = dim(fnode))
  target_a[5, 1, 1] <- 1
  target_b <- array(0, dim = dim(fnode))
  target_b[5, 2, 1] <- 1

  tracks <- list(
    matrix(c(1,1,1, 1,1,1), ncol = 3, byrow = TRUE),
    matrix(c(1,1,1, 5,1,1), ncol = 3, byrow = TRUE),
    matrix(c(1,1,1, 5,2,1), ncol = 3, byrow = TRUE),
    matrix(c(2,1,1, 5,1,1, 5,2,1), ncol = 3, byrow = TRUE),
    matrix(c(2,1,1, 4,2,1), ncol = 3, byrow = TRUE),
    matrix(c(2,1,1, 3,1,1, 5,1,1), ncol = 3, byrow = TRUE),
    matrix(c(3,1,1, 2,1,1, 5,2,1), ncol = 3, byrow = TRUE)
  )

  res <- fibrtraq(
    fNode = fnode,
    localROI = local,
    cNode_list = list(A = target_a, B = target_b),
    cNode_labels = c("A", "B"),
    all_tracks = tracks,
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    focal_assignment_mode = "first_intersection"
  )

  expected <- data.frame(
    fNode_Coords = c("1_1_1", "2_1_1", "3_1_1"),
    tracks = c(3L, 3L, 1L),
    local = c(1L, 0L, 0L),
    leave_no_target = c(0L, 1L, 0L),
    A = c(1L, 2L, 0L),
    B = c(1L, 1L, 1L),
    none = c(0L, 0L, 0L),
    stringsAsFactors = FALSE
  )
  ftq_expect_core_equal(res, expected)

  expect_warning(summary <- rsc_summary(res, percent = TRUE), "deprecated")
  expect_equal(summary$node, c("A", "B"))
  expect_equal(summary$tracks, c(3, 3))
  expect_equal(summary$target_associations, c(3, 3))
  expect_equal(summary$focal_voxel_streamline_target_associations, c(3, 3))
  expect_equal(summary$voxels_connected, c(2L, 3L))
  expect_equal(summary$rsc_percent, c(200 / 3, 100))

  profiled <- fibrtraq(
    fNode = fnode,
    localROI = local,
    cNode_list = list(A = target_a, B = target_b),
    cNode_labels = c("A", "B"),
    all_tracks = tracks,
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    focal_assignment_mode = "first_intersection",
    profile = TRUE
  )
  profile <- attr(profiled, "fibrtraq_profile")
  attr(profiled, "fibrtraq_profile") <- NULL
  attr(profiled, "fibrtraq_backend") <- NULL
  attr(profiled, "fibrtraq_focal_assignment_mode") <- NULL
  attr(profiled, "fibrtraq_target_intersection_mode") <- NULL
  attr(profiled, "fibrtraq_track_space") <- NULL
  attr(profiled, "fibrtraq_axis_flip") <- NULL
  attr(profiled, "fibrtraq_discretization") <- NULL
  attr(profiled, "fibrtraq_roi_summary") <- NULL
  attr(profiled, "fibrtraq_mapping_selection_rule") <- NULL
  ftq_expect_core_equal(profiled, expected)
  expect_true(all(c("prepare_rois_and_luts", "mapping_detection", "count_tracks", "assemble_output") %in% profile$stage))
  expect_true(all(profile$seconds >= 0))

  cpp <- fibrtraq(
    fNode = fnode,
    localROI = local,
    cNode_list = list(A = target_a, B = target_b),
    cNode_labels = c("A", "B"),
    all_tracks = tracks,
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    focal_assignment_mode = "first_intersection",
    backend = "cpp"
  )
  ftq_expect_core_equal(cpp, expected)
  ftq_expect_core_equal(cpp, res)
  ftq_expect_algorithm_metadata_equal(cpp, res)

  auto <- fibrtraq(
    fNode = fnode,
    localROI = local,
    cNode_list = list(A = target_a, B = target_b),
    cNode_labels = c("A", "B"),
    all_tracks = tracks,
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    focal_assignment_mode = "first_intersection",
    backend = "auto",
    profile = TRUE
  )
  expect_equal(attr(auto, "fibrtraq_backend"), "cpp")
  attr(auto, "fibrtraq_profile") <- NULL
  attr(auto, "fibrtraq_backend") <- NULL
  ftq_expect_core_equal(auto, res)
})

test_that("fibrtraq applies TrackVis header affine for raw .trk coordinates", {
  skip_if_not_installed("RNifti")

  tmp <- tempdir()
  ref <- file.path(tmp, "trackvis-ref.nii.gz")
  RNifti::writeNifti(array(0, dim = c(3, 3, 3)), ref)

  fnode <- array(0, dim = c(3, 3, 3))
  fnode[1, 1, 1] <- 1
  target <- array(0, dim = dim(fnode))
  target[2, 1, 1] <- 1

  trk <- list(
    header = list(vox2ras = diag(4)),
    tracks = list(
      list(
        num_points = 2L,
        coords = matrix(c(0, 0, 0, 1, 0, 0), ncol = 3, byrow = TRUE)
      )
    )
  )

  res <- fibrtraq(
    fNode = fnode,
    localROI = fnode,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = trk,
    track_space = "trackvis",
    ref_image = ref,
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE)
  )

  expect_equal(res$tracks, 1L)
  expect_equal(res$local, 0L)
  expect_equal(res$leave_no_target, 0L)
  expect_equal(res$T, 1L)

  cpp <- fibrtraq(
    fNode = fnode,
    localROI = fnode,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = trk,
    track_space = "trackvis",
    ref_image = ref,
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    backend = "cpp"
  )
  ftq_expect_core_equal(cpp, res)
  ftq_expect_algorithm_metadata_equal(cpp, res)
})

test_that("compiled backend matches R backend for rounding and axis flips", {
  fnode <- array(0, dim = c(4, 2, 1))
  fnode[2, 1, 1] <- 1
  target <- array(0, dim = dim(fnode))
  target[1, 2, 1] <- 1
  local <- array(0, dim = dim(fnode))
  local[2, 1, 1] <- 1

  tracks <- list(
    matrix(c(2.5, 1, 1, 4, 2, 1), ncol = 3, byrow = TRUE),
    matrix(c(2.5, 1, 1, 3.2, 1, 1), ncol = 3, byrow = TRUE)
  )

  r <- fibrtraq(
    fNode = fnode,
    localROI = local,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = tracks,
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(TRUE, FALSE, FALSE),
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
    axis_flip = c(TRUE, FALSE, FALSE),
    backend = "cpp"
  )
  ftq_expect_core_equal(cpp, r)
  ftq_expect_algorithm_metadata_equal(cpp, r)
})

test_that("streamlines crossing multiple focal voxels are assigned to first focal hit", {
  fnode <- array(0, dim = c(4, 1, 1))
  fnode[1, 1, 1] <- 1
  fnode[2, 1, 1] <- 1
  target <- array(0, dim = dim(fnode))
  target[4, 1, 1] <- 1

  res <- fibrtraq(
    fNode = fnode,
    localROI = array(0, dim = dim(fnode)),
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = list(matrix(c(2,1,1, 1,1,1, 4,1,1), ncol = 3, byrow = TRUE)),
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    focal_assignment_mode = "first_intersection"
  )

  expect_equal(res$tracks, c(0L, 1L))
  expect_equal(res$T, c(0L, 1L))
  expect_equal(res$none, c(1L, 0L))
})
