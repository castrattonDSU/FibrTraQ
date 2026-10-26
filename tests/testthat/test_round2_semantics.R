library(testthat)
library(FibrTraQ)

test_that("first-intersection focal assignment follows stored point order in both backends", {
  fnode <- array(0, dim = c(4, 1, 1))
  fnode[1, 1, 1] <- 1
  fnode[2, 1, 1] <- 1
  target <- array(0, dim = dim(fnode))
  target[4, 1, 1] <- 1
  empty_local <- array(0, dim = dim(fnode))

  tracks <- list(
    matrix(c(1, 1, 1, 2, 1, 1, 4, 1, 1), ncol = 3, byrow = TRUE),
    matrix(c(2, 1, 1, 1, 1, 1, 4, 1, 1), ncol = 3, byrow = TRUE)
  )

  r <- fibrtraq(
    fNode = fnode,
    localROI = empty_local,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = tracks,
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    focal_assignment_mode = "first_intersection",
    backend = "r"
  )

  cpp <- fibrtraq(
    fNode = fnode,
    localROI = empty_local,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = tracks,
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    focal_assignment_mode = "first_intersection",
    backend = "cpp"
  )

  expect_equal(r$tracks, c(1L, 1L))
  expect_equal(r$T, c(1L, 1L))
  ftq_expect_core_equal(cpp, r)
  ftq_expect_algorithm_metadata_equal(cpp, r)
})

test_that("local and target flags are non-exclusive when masks overlap", {
  fnode <- array(0, dim = c(2, 1, 1))
  fnode[1, 1, 1] <- 1
  local <- fnode
  target <- fnode
  track <- matrix(c(1, 1, 1, 1, 1, 1), ncol = 3, byrow = TRUE)

  res <- fibrtraq(
    fNode = fnode,
    localROI = local,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = list(track),
    track_space = "voxel",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    backend = "cpp"
  )

  expect_equal(res$tracks, 1L)
  expect_equal(res$local, 1L)
  expect_equal(res$T, 1L)
  expect_equal(res$leave_no_target, 0L)
})

test_that("discretization any uses one selected method for focal target and local status", {
  fnode <- array(0, dim = c(2, 1, 1))
  fnode[1, 1, 1] <- 1
  local <- array(0, dim = dim(fnode))
  local[2, 1, 1] <- 1
  target <- fnode

  track <- matrix(c(1.6, 1, 1, 1.6, 1, 1), ncol = 3, byrow = TRUE)

  r <- fibrtraq(
    fNode = fnode,
    localROI = local,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = list(track),
    track_space = "voxel",
    discretization = "any",
    axis_flip = c(FALSE, FALSE, FALSE),
    backend = "r"
  )

  cpp <- fibrtraq(
    fNode = fnode,
    localROI = local,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = list(track),
    track_space = "voxel",
    discretization = "any",
    axis_flip = c(FALSE, FALSE, FALSE),
    backend = "cpp"
  )

  expect_equal(nrow(r), 1L)
  expect_equal(r$tracks, 1L)
  expect_equal(r$local, 0L)
  expect_equal(r$T, 1L)
  expect_equal(r$leave_no_target, 0L)
  ftq_expect_core_equal(cpp, r)
  ftq_expect_algorithm_metadata_equal(cpp, r)
})
