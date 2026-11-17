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
