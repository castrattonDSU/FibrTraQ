test_that("iso_tracks identifies tracks touching ROI", {
  track <- list(num_points = 2,
                coords = matrix(c(0.1,0.1,0.1, 0.1,0.1,0.1), ncol = 3, byrow = TRUE))
  fNode <- matrix(c(1,1,1), ncol = 3)
  localROI <- matrix(c(1,1,1), ncol = 3)
  res <- iso_tracks(track, fNode, localROI, verbose = FALSE)
  expect_type(res, "list")
  expect_true(res$local == 1)
})

test_that("cNode_tracks exposes leave_no_target for composed workflows", {
  fnode <- matrix(c(1,1,1), ncol = 3)
  local <- fnode
  target <- matrix(c(2,1,1), ncol = 3)
  tracks <- list(
    iso_tracks(matrix(c(1,1,1, 2,1,1), ncol = 3, byrow = TRUE),
               fnode, local, track_space = "voxel", discretization = "round"),
    iso_tracks(matrix(c(1,1,1, 3,1,1), ncol = 3, byrow = TRUE),
               fnode, local, track_space = "voxel", discretization = "round")
  )

  res <- cNode_tracks(list(T = target), "T", track_relister(tracks, fnode))

  expect_equal(names(res), c("fNode_Coords", "tracks", "local",
                             "leave_no_target", "T", "none"))
  expect_equal(res$tracks, 2L)
  expect_equal(res$leave_no_target, 1L)
  expect_equal(res$T, 1L)
})

test_that("voxelmm converts TrackVis voxel-mm coordinates to 1-based voxel centers", {
  voxel_size <- c(0.1, 0.2, 0.5)
  raw <- rbind(
    c(0.05, 0.10, 0.25),
    c(0.25, 0.30, 0.25)
  )
  track <- list(num_points = 2L, coords = raw)

  fnode <- array(0, dim = c(4, 3, 2))
  fnode[1, 1, 1] <- 1
  target <- array(0, dim = dim(fnode))
  target[3, 2, 1] <- 1

  through <- iso_tracks(
    track = track,
    fNode = fnode,
    localROI = fnode,
    track_space = "voxelmm",
    track_voxel_size = voxel_size,
    discretization = "round"
  )
  expect_equal(unname(as.matrix(through$tracks)), rbind(c(1L, 1L, 1L), c(3L, 2L, 1L)))

  through_auto <- iso_tracks(
    track = track,
    fNode = fnode,
    localROI = fnode,
    track_space = "auto",
    track_voxel_size = voxel_size,
    discretization = "round"
  )
  expect_equal(unname(as.matrix(through_auto$tracks)), unname(as.matrix(through$tracks)))

  trk <- list(header = list(voxel_size = voxel_size), tracks = list(track))
  r_result <- fibrtraq(
    fNode = fnode,
    localROI = fnode,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = trk,
    track_space = "voxelmm",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    backend = "r"
  )
  expect_equal(r_result$tracks, 1L)
  expect_equal(r_result$T, 1L)
  expect_equal(attr(r_result, "fibrtraq_track_voxel_size"), voxel_size)

  cpp_result <- fibrtraq(
    fNode = fnode,
    localROI = fnode,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = trk,
    track_space = "voxelmm",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    backend = "cpp"
  )
  ftq_expect_core_equal(cpp_result, r_result)

  auto_result <- fibrtraq(
    fNode = fnode,
    localROI = fnode,
    cNode_list = list(T = target),
    cNode_labels = "T",
    all_tracks = trk,
    track_space = "auto",
    discretization = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    backend = "r"
  )
  expect_identical(attr(auto_result, "fibrtraq_track_space"), "voxelmm")
  ftq_expect_core_equal(auto_result, r_result)

  expect_error(
    iso_tracks(track, fnode, track_space = "voxelmm"),
    "requires track_voxel_size"
  )
  expect_error(
    fibrtraq(
      fNode = fnode,
      localROI = fnode,
      cNode_list = list(T = target),
      cNode_labels = "T",
      all_tracks = list(track),
      track_space = "voxelmm",
      axis_flip = c(FALSE, FALSE, FALSE)
    ),
    "requires TrackVis header voxel sizes"
  )
})
