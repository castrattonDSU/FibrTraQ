library(testthat)
library(FibrTraQ)

ftq_round3_run <- function(fnode, local, targets, tracks, mode, backend = "r",
                           discretization = "round", axis_flip = c(FALSE, FALSE, FALSE),
                           reverse_streamlines = FALSE) {
  fibrtraq(
    fNode = fnode,
    localROI = local,
    cNode_list = targets,
    cNode_labels = names(targets),
    all_tracks = tracks,
    track_space = "voxel",
    discretization = discretization,
    axis_flip = axis_flip,
    focal_assignment_mode = mode,
    reverse_streamlines = reverse_streamlines,
    backend = backend
  )
}

test_that("all_focal_voxels is reversal invariant and suppresses duplicate focal visits", {
  fnode <- array(0, dim = c(5, 1, 1))
  fnode[1:3, 1, 1] <- 1
  local <- array(0, dim = dim(fnode))
  target_a <- array(0, dim = dim(fnode))
  target_a[5, 1, 1] <- 1
  target_b <- array(0, dim = dim(fnode))
  target_b[4, 1, 1] <- 1

  track <- matrix(
    c(1, 1, 1,
      2, 1, 1,
      2, 1, 1,
      3, 1, 1,
      2, 1, 1,
      5, 1, 1,
      4, 1, 1),
    ncol = 3,
    byrow = TRUE
  )
  reversed_track <- track[nrow(track):1, , drop = FALSE]

  r <- ftq_round3_run(
    fnode, local, list(A = target_a, B = target_b),
    list(track), "all_focal_voxels", backend = "r"
  )
  cpp <- ftq_round3_run(
    fnode, local, list(A = target_a, B = target_b),
    list(track), "all_focal_voxels", backend = "cpp"
  )
  reversed <- ftq_round3_run(
    fnode, local, list(A = target_a, B = target_b),
    list(reversed_track), "all_focal_voxels", backend = "cpp"
  )

  expected <- data.frame(
    fNode_Coords = c("1_1_1", "2_1_1", "3_1_1"),
    tracks = c(1L, 1L, 1L),
    local = c(0L, 0L, 0L),
    leave_no_target = c(0L, 0L, 0L),
    A = c(1L, 1L, 1L),
    B = c(1L, 1L, 1L),
    none = c(0L, 0L, 0L),
    stringsAsFactors = FALSE
  )
  ftq_expect_core_equal(r, expected)
  ftq_expect_core_equal(cpp, r)
  ftq_expect_core_equal(reversed, r)
  ftq_expect_algorithm_metadata_equal(cpp, r)
  ftq_expect_algorithm_metadata_equal(reversed, r)

  summary <- attr(r, "fibrtraq_roi_summary")
  expect_equal(summary$unique_focal_intersecting_streamlines, 1)
  expect_equal(summary$`unique_target_streamlines:A`, 1)
  expect_equal(summary$`unique_target_streamlines:B`, 1)
  expect_equal(summary$total_focal_voxel_streamline_associations, 3)
	  expect_equal(summary$total_focal_voxel_streamline_target_associations, 6)
	  expect_equal(summary$focal_assignment_mode, "all_focal_voxels")
	  expect_equal(summary$streamline_voxelization_mode, "sample_points")
	})

test_that("first_intersection preserves point-order dependent legacy assignment", {
  fnode <- array(0, dim = c(4, 1, 1))
  fnode[1:2, 1, 1] <- 1
  local <- array(0, dim = dim(fnode))
  target <- array(0, dim = dim(fnode))
  target[4, 1, 1] <- 1

  forward <- matrix(c(1, 1, 1, 2, 1, 1, 4, 1, 1), ncol = 3, byrow = TRUE)
  reversed <- forward[nrow(forward):1, , drop = FALSE]

  fwd <- ftq_round3_run(fnode, local, list(T = target), list(forward),
                        "first_intersection", backend = "cpp")
  rev <- ftq_round3_run(fnode, local, list(T = target), list(reversed),
                        "first_intersection", backend = "r")
  fwd_r <- ftq_round3_run(fnode, local, list(T = target), list(forward),
                          "first_intersection", backend = "r")
  fwd_reversed_in_memory <- ftq_round3_run(
    fnode, local, list(T = target), list(forward),
    "first_intersection", backend = "cpp", reverse_streamlines = TRUE
  )

  expect_equal(fwd$tracks, c(1L, 0L))
  expect_equal(rev$tracks, c(0L, 1L))
  expect_equal(fwd_reversed_in_memory$tracks, c(0L, 1L))
  expect_equal(fwd$T, c(1L, 0L))
  expect_equal(rev$T, c(0L, 1L))
  expect_equal(fwd_reversed_in_memory$T, c(0L, 1L))
  ftq_expect_core_equal(fwd, fwd_r)
  ftq_expect_algorithm_metadata_equal(fwd, fwd_r)
  expect_true(attr(fwd_reversed_in_memory, "fibrtraq_reverse_streamlines"))

  expect_equal(attr(fwd, "fibrtraq_roi_summary")$total_focal_voxel_streamline_associations, 1)
})

test_that("fate fields handle local-only target-only overlap and leave-no-target cases", {
  fnode <- array(0, dim = c(4, 1, 1))
  fnode[1, 1, 1] <- 1
  local <- array(0, dim = dim(fnode))
  local[1:2, 1, 1] <- 1
  target_a <- array(0, dim = dim(fnode))
  target_a[2, 1, 1] <- 1
  target_b <- array(0, dim = dim(fnode))
  target_b[3, 1, 1] <- 1

  tracks <- list(
    matrix(c(1, 1, 1, 1, 1, 1), ncol = 3, byrow = TRUE),
    matrix(c(1, 1, 1, 2, 1, 1), ncol = 3, byrow = TRUE),
    matrix(c(1, 1, 1, 3, 1, 1), ncol = 3, byrow = TRUE),
    matrix(c(1, 1, 1, 4, 1, 1), ncol = 3, byrow = TRUE)
  )

  r <- ftq_round3_run(fnode, local, list(A = target_a, B = target_b), tracks,
                      "all_focal_voxels", backend = "r")
  cpp <- ftq_round3_run(fnode, local, list(A = target_a, B = target_b), tracks,
                        "all_focal_voxels", backend = "cpp")

  expect_equal(r$tracks, 4L)
  expect_equal(r$local, 2L)
  expect_equal(r$leave_no_target, 1L)
  expect_equal(r$A, 1L)
  expect_equal(r$B, 1L)
  ftq_expect_core_equal(cpp, r)
  ftq_expect_algorithm_metadata_equal(cpp, r)

  summary <- attr(r, "fibrtraq_roi_summary")
  expect_equal(summary$unique_focal_intersecting_streamlines, 4)
  expect_equal(summary$unique_local_streamlines, 2)
	  expect_equal(summary$unique_leave_no_target_streamlines, 1)
	  expect_equal(summary$unique_focal_intersecting_with_target_streamlines, 2)
	  expect_equal(summary$unique_local_target_positive_streamlines, 1)
	  expect_equal(summary$unique_local_target_negative_streamlines, 1)
	  expect_equal(summary$unique_nonlocal_target_positive_streamlines, 1)
	  expect_equal(summary$unique_nonlocal_target_negative_streamlines, 1)
	  expect_equal(summary$unique_any_target_positive_streamlines, 2)
	  expect_equal(summary$unique_any_target_negative_streamlines, 2)
	  expect_equal(summary$unique_exclusive_fate_sum, summary$unique_focal_intersecting_streamlines)
	  expect_true(summary$unique_exclusive_fate_sum_matches_focal)
	  expect_equal(summary$`unique_target_streamlines:A`, 1)
	  expect_equal(summary$`unique_target_streamlines:B`, 1)
	})

test_that("all_focal_voxels keeps empty focal voxels and zero-hit runs explicit", {
  fnode <- array(0, dim = c(3, 1, 1))
  fnode[1:2, 1, 1] <- 1
  local <- array(0, dim = dim(fnode))
  target <- array(0, dim = dim(fnode))
  target[3, 1, 1] <- 1

  res <- ftq_round3_run(
    fnode, local, list(T = target),
    list(matrix(c(3, 1, 1, 3, 1, 1), ncol = 3, byrow = TRUE)),
    "all_focal_voxels",
    backend = "cpp"
  )

  expect_equal(res$tracks, c(0L, 0L))
  expect_equal(res$none, c(1L, 1L))
  summary <- attr(res, "fibrtraq_roi_summary")
  expect_equal(summary$processed_streamlines, 1)
  expect_equal(summary$unique_focal_intersecting_streamlines, 0)
  expect_equal(summary$total_focal_voxel_streamline_associations, 0)
})

test_that("invalid focal assignment modes are rejected", {
  fnode <- array(0, dim = c(1, 1, 1))
  fnode[1, 1, 1] <- 1
  target <- fnode
  track <- matrix(c(1, 1, 1, 1, 1, 1), ncol = 3, byrow = TRUE)

  expect_error(
    fibrtraq(
      fNode = fnode,
      localROI = fnode,
      cNode_list = list(T = target),
      cNode_labels = "T",
      all_tracks = list(track),
      track_space = "voxel",
      discretization = "round",
      axis_flip = c(FALSE, FALSE, FALSE),
      focal_assignment_mode = "bad_mode"
    ),
    "should be one of"
  )
})
