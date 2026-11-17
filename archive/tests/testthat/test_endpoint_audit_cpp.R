test_that("endpoint audit separates intersections from endpoint pairs", {
  vol_dim <- c(7L, 3L, 3L)
  n_voxels <- prod(vol_dim)

  region_a <- array(FALSE, dim = vol_dim)
  region_b <- array(FALSE, dim = vol_dim)
  region_a[2, 2, 2] <- TRUE
  region_b[5, 2, 2] <- TRUE

  region_a_dilated <- array(FALSE, dim = vol_dim)
  region_b_dilated <- array(FALSE, dim = vol_dim)
  region_a_dilated[1:3, 1:3, 1:3] <- TRUE
  region_b_dilated[4:6, 1:3, 1:3] <- TRUE

  context <- array(FALSE, dim = vol_dim)
  context[1:6, , ] <- TRUE

  tracks <- list(
    rbind(c(2, 2, 2), c(3, 2, 2), c(4, 2, 2), c(5, 2, 2)),
    rbind(c(1, 2, 2), c(2, 2, 2), c(3, 2, 2), c(4, 2, 2), c(5, 2, 2)),
    rbind(c(2, 2, 2), c(3, 2, 2), c(4, 2, 2), c(5, 2, 2), c(6, 2, 2)),
    rbind(c(1, 2, 2), c(7, 2, 2))
  )

  result <- ftq_endpoint_audit_chunk_cpp(
    tracks = tracks,
    raw_to_voxel = diag(4),
    method = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    vol_dim = vol_dim,
    region_lut = matrix(
      c(as.vector(region_a), as.vector(region_b)),
      nrow = n_voxels
    ),
    region_lut_dilated = matrix(
      c(as.vector(region_a_dilated), as.vector(region_b_dilated)),
      nrow = n_voxels
    ),
    context_lut = matrix(as.vector(context), nrow = n_voxels),
    streamline_voxelization_mode = "segment_interpolation",
    penetration_hist_max = 10L
  )

  expect_equal(result$processed, 4L)
  expect_equal(result$valid_streamlines, 4L)
  expect_equal(result$endpoint_count, 8L)
  expect_equal(result$endpoint_out_of_bounds, 0L)
  expect_equal(result$endpoint_out_of_bounds_low, c(0L, 0L, 0L))
  expect_equal(result$endpoint_out_of_bounds_high, c(0L, 0L, 0L))
  expect_equal(result$endpoint_out_of_bounds_low_distance, matrix(0L, 3, 2))
  expect_equal(result$endpoint_out_of_bounds_high_distance, matrix(0L, 3, 2))
  expect_equal(result$region_intersection, c(4L, 4L))
  expect_equal(result$endpoint_region_exact, c(2L, 2L))
  expect_equal(result$endpoint_region_dilated, c(4L, 3L))

  expect_equal(result$intersection_pair[1, 2], 4L)
  expect_equal(result$intersection_pair_target_endpoint_exact[1, 2], 2L)
  expect_equal(result$intersection_pair_focal_endpoint_exact[1, 2], 2L)
  expect_equal(result$intersection_pair_any_relevant_endpoint_exact[1, 2], 3L)
  expect_equal(result$intersection_pair_any_relevant_endpoint_dilated[1, 2], 4L)
  expect_equal(result$focal_endpoint_target_intersection_exact[1, 2], 2L)
  expect_equal(result$endpoint_pair_exact[1, 2], 1L)
  expect_equal(result$endpoint_pair_dilated[1, 2], 3L)
  expect_equal(result$intersection_pair_endpoint_pair_dilated[1, 2], 3L)
  expect_equal(result$endpoint_class_exact[2, ], c(0L, 1L, 1L, 1L, 1L))

  expect_equal(result$endpoint_context, 7L)
  expect_equal(result$streamline_any_endpoint_context, 4L)
  expect_equal(result$streamline_both_endpoints_context, 3L)
  expect_equal(result$penetration_hist[2, 2], 4L)
})

test_that("endpoint audit reports the image boundary crossed by endpoints", {
  vol_dim <- c(7L, 3L, 3L)
  n_voxels <- prod(vol_dim)
  result <- ftq_endpoint_audit_chunk_cpp(
    tracks = list(
      rbind(c(0, 2, 2), c(8, 4, 2)),
      matrix(c(-2, 2, 2), nrow = 1)
    ),
    raw_to_voxel = diag(4),
    method = "round",
    axis_flip = c(FALSE, FALSE, FALSE),
    vol_dim = vol_dim,
    region_lut = matrix(FALSE, nrow = n_voxels, ncol = 1),
    region_lut_dilated = matrix(FALSE, nrow = n_voxels, ncol = 1),
    context_lut = matrix(FALSE, nrow = n_voxels, ncol = 1),
    streamline_voxelization_mode = "segment_interpolation",
    penetration_hist_max = 10L
  )

  expect_equal(result$endpoint_count, 3L)
  expect_equal(result$endpoint_out_of_bounds, 3L)
  expect_equal(result$endpoint_out_of_bounds_low, c(2L, 0L, 0L))
  expect_equal(result$endpoint_out_of_bounds_high, c(1L, 1L, 0L))
  expect_equal(
    result$endpoint_out_of_bounds_low_distance,
    rbind(c(1L, 1L), c(0L, 0L), c(0L, 0L))
  )
  expect_equal(
    result$endpoint_out_of_bounds_high_distance,
    rbind(c(1L, 0L), c(1L, 0L), c(0L, 0L))
  )
})

test_that("endpoint audit validates lookup dimensions and modes", {
  expect_error(
    ftq_endpoint_audit_chunk_cpp(
      tracks = list(matrix(c(1, 1, 1), nrow = 1)),
      raw_to_voxel = diag(4),
      method = "round",
      axis_flip = c(FALSE, FALSE, FALSE),
      vol_dim = c(2L, 2L, 2L),
      region_lut = matrix(FALSE, nrow = 7, ncol = 1),
      region_lut_dilated = matrix(FALSE, nrow = 7, ncol = 1),
      context_lut = matrix(FALSE, nrow = 7, ncol = 1),
      streamline_voxelization_mode = "sample_points"
    ),
    "one row per reference voxel"
  )

  expect_error(
    ftq_endpoint_audit_chunk_cpp(
      tracks = list(matrix(c(1, 1, 1), nrow = 1)),
      raw_to_voxel = diag(4),
      method = "round",
      axis_flip = c(FALSE, FALSE, FALSE),
      vol_dim = c(2L, 2L, 2L),
      region_lut = matrix(FALSE, nrow = 8, ncol = 1),
      region_lut_dilated = matrix(FALSE, nrow = 8, ncol = 1),
      context_lut = matrix(FALSE, nrow = 8, ncol = 1),
      streamline_voxelization_mode = "unsupported"
    ),
    "Unsupported streamline_voxelization_mode"
  )
})
