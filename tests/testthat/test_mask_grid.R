library(testthat)
library(FibrTraQ)

ftq_test_nifti_image <- function(arr,
                                 shift = c(0, 0, 0),
                                 scale = c(1, 1, 1),
                                 spatial_units = NULL) {
  scale <- rep_len(as.numeric(scale), 3L)
  shift <- rep_len(as.numeric(shift), 3L)
  image <- RNifti::asNifti(
    arr,
    reference = list(
      sform_code = 2L,
      qform_code = 0L,
      srow_x = c(scale[[1L]], 0, 0, shift[[1L]]),
      srow_y = c(0, scale[[2L]], 0, shift[[2L]]),
      srow_z = c(0, 0, scale[[3L]], shift[[3L]])
    )
  )
  if (!is.null(spatial_units)) RNifti::pixunits(image) <- spatial_units
  image
}

ftq_test_write_nifti <- function(path, arr, shift = c(0, 0, 0), datatype = "uint8") {
  RNifti::writeNifti(
    ftq_test_nifti_image(arr, shift = shift),
    path,
    datatype = datatype
  )
  invisible(path)
}

test_that("mask grid report distinguishes dimension and affine mismatch", {
  skip_if_not_installed("RNifti")

  tmp <- tempfile("ftq-mask-grid-")
  dir.create(tmp)
  template <- file.path(tmp, "template.nii.gz")
  mask <- file.path(tmp, "mask.nii.gz")
  shifted <- file.path(tmp, "shifted.nii.gz")
  mismatched <- file.path(tmp, "mismatched.nii.gz")
  reference <- file.path(tmp, "reference.nii.gz")
  report_csv <- file.path(tmp, "grid-report.csv")

  arr <- array(0, dim = c(3, 3, 3))
  arr[2, 2, 2] <- 1
  ref_arr <- array(0, dim = c(3, 3, 3))
  ref_arr[2, 2, 2] <- 1
  ref_arr[3, 3, 3] <- 1

  ftq_test_write_nifti(template, array(0, dim = c(3, 3, 3)))
  ftq_test_write_nifti(mask, arr)
  ftq_test_write_nifti(shifted, arr, shift = c(1, 0, 0))
  ftq_test_write_nifti(mismatched, array(0, dim = c(4, 3, 3)))
  ftq_test_write_nifti(reference, ref_arr)

  ok <- ftq_mask_grid_report(mask, template, reference_masks = c(ref = reference), output_csv = report_csv)
  expect_true(file.exists(report_csv))
  expect_true(ok$dims_match)
  expect_true(ok$affine_match)
  expect_true(ok$grid_match)
  expect_equal(ok$nonzero_voxels, 1L)
  expect_equal(ok$voxel_bbox_x_min, 2L)
  expect_equal(ok$voxel_bbox_z_max, 2L)
  expect_equal(ok$overlap_voxels, 1L)
  expect_equal(ok$union_voxels, 2L)
  expect_equal(ok$jaccard, 0.5)

  shifted_report <- ftq_mask_grid_report(shifted, template)
  expect_true(shifted_report$dims_match)
  expect_false(shifted_report$affine_match)
  expect_false(shifted_report$grid_match)

  mismatched_report <- ftq_mask_grid_report(mismatched, template)
  expect_false(mismatched_report$dims_match)
  expect_false(mismatched_report$grid_match)
})

test_that("union mask requires same grid and records overlap metadata", {
  skip_if_not_installed("RNifti")

  tmp <- tempfile("ftq-union-")
  dir.create(tmp)
  template <- file.path(tmp, "template.nii.gz")
  a <- file.path(tmp, "a.nii.gz")
  b <- file.path(tmp, "b.nii.gz")
  shifted <- file.path(tmp, "shifted.nii.gz")
  output <- file.path(tmp, "union.nii.gz")

  arr_a <- array(0, dim = c(3, 3, 3))
  arr_a[1, 1, 1] <- 1
  arr_a[2, 2, 2] <- 1
  arr_b <- array(0, dim = c(3, 3, 3))
  arr_b[2, 2, 2] <- 1
  arr_b[3, 3, 3] <- 1

  ftq_test_write_nifti(template, array(0, dim = c(3, 3, 3)))
  ftq_test_write_nifti(a, arr_a)
  ftq_test_write_nifti(b, arr_b)
  ftq_test_write_nifti(shifted, arr_b, shift = c(1, 0, 0))

  union <- ftq_make_union_mask(c(A = a, B = b), template, output, label = "test-union")
  expect_true(file.exists(union$path))
  expect_equal(union$summary$union_voxels, 3L)
  expect_equal(union$inputs$nonzero_voxels, c(2L, 2L))
  expect_equal(union$pairwise_overlaps$overlap_voxels, 1L)
  expect_equal(sum(as.array(RNifti::readNifti(output, internal = FALSE)) != 0), 3L)
  expect_error(
    ftq_make_union_mask(c(A = a, shifted = shifted), template, file.path(tmp, "bad-union.nii.gz")),
    "affine/grid differs"
  )
  expect_error(
    ftq_make_union_mask(c(A = a), template, output),
    "already exists"
  )
})

test_that("regrid mask preserves identity grids and handles trusted affine translations", {
  skip_if_not_installed("RNifti")

  tmp <- tempfile("ftq-regrid-")
  dir.create(tmp)
  template <- file.path(tmp, "template.nii.gz")
  source <- file.path(tmp, "source.nii.gz")
  shifted_source <- file.path(tmp, "shifted-source.nii.gz")
  identity_out <- file.path(tmp, "identity-out.nii.gz")
  shifted_out <- file.path(tmp, "shifted-out.nii.gz")

  arr <- array(0, dim = c(3, 3, 3))
  arr[1, 1, 1] <- 1
  arr[2, 2, 2] <- 1
  ftq_test_write_nifti(template, array(0, dim = c(3, 3, 3)))
  ftq_test_write_nifti(source, arr)
  ftq_test_write_nifti(shifted_source, arr, shift = c(1, 0, 0))

  identity <- ftq_regrid_mask(source, template, identity_out)
  expect_true(file.exists(identity$path))
  expect_equal(as.array(RNifti::readNifti(identity_out, internal = FALSE)) != 0, arr != 0)
  expect_equal(identity$summary$output_source_voxel_ratio, 1)

  shifted <- ftq_regrid_mask(shifted_source, template, shifted_out)
  shifted_arr <- as.array(RNifti::readNifti(shifted_out, internal = FALSE)) != 0
  expected <- array(FALSE, dim = c(3, 3, 3))
  expected[2, 1, 1] <- TRUE
  expected[3, 2, 2] <- TRUE
  expect_equal(shifted_arr, expected)
})

test_that("regrid mask warns on suspicious voxel loss", {
  skip_if_not_installed("RNifti")

  tmp <- tempfile("ftq-regrid-warning-")
  dir.create(tmp)
  source <- file.path(tmp, "source.nii.gz")
  template <- file.path(tmp, "tiny-template.nii.gz")
  output <- file.path(tmp, "out.nii.gz")

  ftq_test_write_nifti(source, array(1, dim = c(3, 3, 3)))
  ftq_test_write_nifti(template, array(0, dim = c(2, 2, 2)))

  expect_warning(
    out <- ftq_regrid_mask(source, template, output),
    "Suspicious voxel-count change"
  )
  expect_true(out$summary$output_source_voxel_ratio < 0.5)
})

test_that("preflight reports harmonization guidance for same-dimension affine mismatch", {
  skip_if_not_installed("RNifti")

  tmp <- tempfile("ftq-preflight-grid-")
  dir.create(tmp)
  fnode <- file.path(tmp, "fnode.nii.gz")
  target <- file.path(tmp, "target.nii.gz")
  shifted_target <- file.path(tmp, "shifted-target.nii.gz")
  unit_mismatched_target <- file.path(tmp, "unit-mismatched-target.nii.gz")
  region_a <- file.path(tmp, "region-a.nii.gz")
  region_b <- file.path(tmp, "region-b.nii.gz")
  local <- file.path(tmp, "local-union.nii.gz")
  tractogram <- file.path(tmp, "tracts.trk")

  fnode_arr <- array(0, dim = c(3, 3, 3))
  fnode_arr[1, 1, 1] <- 1
  target_arr <- array(0, dim = c(3, 3, 3))
  target_arr[3, 3, 3] <- 1
  region_a_arr <- fnode_arr
  region_b_arr <- target_arr

  ftq_test_write_nifti(fnode, fnode_arr)
  ftq_test_write_nifti(target, target_arr)
  ftq_test_write_nifti(shifted_target, target_arr, shift = c(1, 0, 0))
  unit_mismatched_image <- RNifti::asNifti(
    target_arr,
    reference = list(
      sform_code = 2L,
      qform_code = 0L,
      srow_x = c(1, 0, 0, 0),
      srow_y = c(0, 1, 0, 0),
      srow_z = c(0, 0, 1, 0)
    )
  )
  RNifti::pixunits(unit_mismatched_image) <- "m"
  RNifti::writeNifti(unit_mismatched_image, unit_mismatched_target)
  ftq_test_write_nifti(region_a, region_a_arr)
  ftq_test_write_nifti(region_b, region_b_arr)
  file.create(tractogram)

  ftq_make_union_mask(c(focal = region_a, target = region_b), fnode, local)

  manifest <- list(
    dataset_root = tmp,
    subject = "sub-001",
    description = "preflight-grid",
    fnode = "fnode.nii.gz",
    local = "local-union.nii.gz",
    tractogram = "tracts.trk",
    targets = list(T = "target.nii.gz"),
    outputs = list(output_dir = file.path(tmp, "out"), overwrite = FALSE)
  )
  preflight <- ftq_preflight(manifest)
  expect_true(preflight$passed)
  expect_true(preflight$checks$passed[preflight$checks$check == "local_affine_matches_fnode"])

  bad_manifest <- manifest
  bad_manifest$targets <- list(T = "shifted-target.nii.gz")
  bad <- ftq_preflight(bad_manifest)
  expect_false(bad$passed)
  affine_row <- bad$checks[bad$checks$check == "target_affine_matches_fnode:T", , drop = FALSE]
  expect_equal(nrow(affine_row), 1L)
  expect_false(affine_row$passed)
  expect_match(affine_row$detail, "dimensions match but affine/grid differs")
  expect_match(affine_row$detail, "ftq_make_union_mask", fixed = TRUE)
  expect_match(affine_row$detail, "ftq_regrid_mask", fixed = TRUE)

  unit_manifest <- manifest
  unit_manifest$targets <- list(T = "unit-mismatched-target.nii.gz")
  unit_mismatch <- ftq_preflight(unit_manifest)
  unit_affine_row <- unit_mismatch$checks[
    unit_mismatch$checks$check == "target_affine_matches_fnode:T",
    ,
    drop = FALSE
  ]
  expect_false(unit_mismatch$passed)
  expect_equal(nrow(unit_affine_row), 1L)
  expect_false(unit_affine_row$passed)
  expect_match(
    unit_affine_row$detail,
    "dimensions match but affine/grid differs"
  )

  local_search_manifest <- manifest
  local_search_manifest$ref_image <- "shifted-target.nii.gz"
  local_search_manifest$analysis <- list(
    target_intersection_mode = "local_search",
    local_search_radius_mm = 2
  )
  local_search <- ftq_preflight(local_search_manifest)
  expect_false(local_search$passed)
  expect_true(local_search$checks$passed[
    local_search$checks$check == "local_search_geometry_available"
  ])
  unit_row <- local_search$checks[
    local_search$checks$check == "local_search_spatial_units_known",
    ,
    drop = FALSE
  ]
  expect_equal(nrow(unit_row), 1L)
  expect_equal(unit_row$severity, "warning")
  expect_false(unit_row$passed)
  expect_match(unit_row$detail, "Unknown", fixed = TRUE)
  expect_match(unit_row$detail, "scale_to_mm = 1", fixed = TRUE)
  expect_match(unit_row$detail, "unknown_assumed_mm", fixed = TRUE)
  expect_true(
    "local_search_spatial_units_known" %in%
      FibrTraQ:::ftq_compact_preflight(local_search)$checks$check
  )
  reference_affine <- local_search$checks[
    local_search$checks$check ==
      "local_search_ref_image_affine_matches_fnode",
    ,
    drop = FALSE
  ]
  expect_equal(nrow(reference_affine), 1L)
  expect_false(reference_affine$passed)
  expect_match(
    reference_affine$detail,
    "dimensions match but affine/grid differs"
  )
  expect_error(
    ftq_run(local_search_manifest, verbose = FALSE),
    "preflight failed",
    fixed = TRUE
  )
})

test_that("direct FibrTraQ validates ROI dimensions, geometry, and coordinates", {
  skip_if_not_installed("RNifti")

  dims <- c(3L, 3L, 3L)
  focal_array <- array(0, dim = dims)
  focal_array[1, 1, 1] <- 1
  target_array <- array(0, dim = dims)
  target_array[3, 1, 1] <- 1
  tracks <- list(rbind(c(1, 1, 1), c(3, 1, 1)))
  run_direct <- function(focal = focal_array,
                         local = NULL,
                         target = target_array,
                         ...) {
    FibrTraQ(
      fNode = focal,
      localROI = local,
      cNode_list = list(T = target),
      cNode_labels = "T",
      all_tracks = tracks,
      discretization = "round",
      focal_assignment_mode = "all_focal_voxels",
      streamline_voxelization_mode = "sample_points",
      track_space = "voxel",
      axis_flip = c(FALSE, FALSE, FALSE),
      backend = "r",
      endpoint_analysis = "none",
      ...
    )
  }

  expect_error(
    run_direct(local = array(0, dim = c(4, 3, 3))),
    "localROI dimensions must match the fNode grid",
    fixed = TRUE
  )
  expect_error(
    run_direct(target = array(0, dim = c(3, 4, 3))),
    "cNode_list[[1]] ('T') dimensions must match the fNode grid",
    fixed = TRUE
  )
  expect_error(
    run_direct(focal = array(0, dim = c(3, 3, 3, 2))),
    "fNode must be a 3-D array",
    fixed = TRUE
  )
  expect_error(
    run_direct(local = array(0, dim = c(3, 3, 3, 2))),
    "localROI must be a 3-D array",
    fixed = TRUE
  )
  expect_error(
    run_direct(target = array(0, dim = c(3, 3, 3, 2))),
    "cNode_list[[1]] ('T') must be a 3-D array",
    fixed = TRUE
  )

  focal_image <- ftq_test_nifti_image(focal_array)
  target_image <- ftq_test_nifti_image(target_array)
  expect_error(
    run_direct(
      focal = focal_image,
      local = ftq_test_nifti_image(focal_array, shift = c(1, 0, 0)),
      target = target_image
    ),
    "localROI must use the same unit-normalized affine/grid as fNode",
    fixed = TRUE
  )
  expect_error(
    run_direct(
      focal = focal_image,
      target = ftq_test_nifti_image(target_array, shift = c(1, 0, 0))
    ),
    "cNode_list[[1]] ('T') must use the same unit-normalized affine/grid as fNode",
    fixed = TRUE
  )

  # A NIfTI object with transform code 0 has no declared affine provenance.
  # RNifti supplies an identity fallback for calculations, but validation must
  # not present that fallback as a measured grid mismatch.
  undefined_focal_grid <- RNifti::asNifti(focal_array)
  expect_equal(attr(RNifti::xform(undefined_focal_grid), "code"), 0)
  undefined_grid_result <- run_direct(
    focal = undefined_focal_grid,
    target = ftq_test_nifti_image(target_array, shift = c(1, 0, 0))
  )
  expect_equal(undefined_grid_result$T, 1L)

  expect_error(
    run_direct(local = matrix(c(1.5, 1, 1), nrow = 1L)),
    "localROI coordinates must be integer-compatible voxel indices",
    fixed = TRUE
  )
  expect_error(
    run_direct(target = matrix(c(Inf, 1, 1), nrow = 1L)),
    "coordinates must all be finite numeric values",
    fixed = TRUE
  )
  expect_error(
    run_direct(target = matrix(c(4, 1, 1), nrow = 1L)),
    "coordinates must be within the fNode voxel grid",
    fixed = TRUE
  )

  coordinate_result <- run_direct(
    focal = focal_image,
    local = rbind(c(1, 1, 1), c(3, 1, 1)),
    target = matrix(c(3, 1, 1), nrow = 1L)
  )
  expect_equal(coordinate_result$local, 1L)
  expect_equal(coordinate_result$T, 1L)

  focal_mm <- ftq_test_nifti_image(
    focal_array,
    scale = 1000,
    spatial_units = "mm"
  )
  target_m <- ftq_test_nifti_image(
    target_array,
    scale = 1,
    spatial_units = "m"
  )
  expect_no_error(run_direct(focal = focal_mm, target = target_m))
})

test_that("preflight rejects four-dimensional masks without array errors", {
  skip_if_not_installed("RNifti")

  tmp <- tempfile("ftq-preflight-exact-3d-")
  dir.create(tmp)
  fnode_3d <- file.path(tmp, "fnode-3d.nii.gz")
  fnode_4d <- file.path(tmp, "fnode-4d.nii.gz")
  target_3d <- file.path(tmp, "target-3d.nii.gz")
  target_4d <- file.path(tmp, "target-4d.nii.gz")
  local_4d <- file.path(tmp, "local-4d.nii.gz")
  context_4d <- file.path(tmp, "context-4d.nii.gz")
  tractogram <- file.path(tmp, "tracts.trk")

  mask_3d <- array(0, dim = c(3, 3, 3))
  mask_3d[1, 1, 1] <- 1
  mask_4d <- array(0, dim = c(3, 3, 3, 2))
  mask_4d[1, 1, 1, 1] <- 1
  ftq_test_write_nifti(fnode_3d, mask_3d)
  ftq_test_write_nifti(fnode_4d, mask_4d)
  ftq_test_write_nifti(target_3d, mask_3d)
  ftq_test_write_nifti(target_4d, mask_4d)
  ftq_test_write_nifti(local_4d, mask_4d)
  ftq_test_write_nifti(context_4d, mask_4d)
  file.create(tractogram)

  manifest <- list(
    dataset_root = tmp,
    subject = "sub-exact-3d",
    description = "exact-3d",
    fnode = basename(fnode_3d),
    tractogram = basename(tractogram),
    targets = list(T = basename(target_3d)),
    outputs = list(
      output_dir = file.path(tmp, "out"),
      overwrite = FALSE
    )
  )

  all_4d_manifest <- manifest
  all_4d_manifest$fnode <- basename(fnode_4d)
  all_4d_manifest$targets <- list(T = basename(target_4d))
  expect_no_error(all_4d <- ftq_preflight(all_4d_manifest))
  expect_false(all_4d$passed)
  expect_false(all_4d$checks$passed[all_4d$checks$check == "fnode_is_3d"])
  expect_false(all_4d$checks$passed[
    all_4d$checks$check == "target_is_3d:T"
  ])
  expect_match(
    all_4d$checks$detail[all_4d$checks$check == "fnode_is_3d"],
    "3x3x3x2",
    fixed = TRUE
  )

  target_4d_manifest <- manifest
  target_4d_manifest$targets <- list(T = basename(target_4d))
  expect_no_error(target_4d_preflight <- ftq_preflight(target_4d_manifest))
  expect_false(target_4d_preflight$passed)
  target_grid_row <- target_4d_preflight$checks[
    target_4d_preflight$checks$check == "target_dims_match_fnode:T",
    ,
    drop = FALSE
  ]
  expect_equal(nrow(target_grid_row), 1L)
  expect_false(target_grid_row$passed)
  expect_match(target_grid_row$detail, "must be exactly 3-D", fixed = TRUE)
  expect_match(target_grid_row$detail, "3x3x3x2", fixed = TRUE)

  local_4d_manifest <- manifest
  local_4d_manifest$local <- basename(local_4d)
  expect_no_error(local_4d_preflight <- ftq_preflight(local_4d_manifest))
  expect_false(local_4d_preflight$passed)
  expect_false(local_4d_preflight$checks$passed[
    local_4d_preflight$checks$check == "local_is_3d"
  ])

  context_4d_manifest <- manifest
  context_4d_manifest$endpoint_contexts <- list(
    anatomy = basename(context_4d)
  )
  expect_no_error(context_4d_preflight <- ftq_preflight(context_4d_manifest))
  expect_false(context_4d_preflight$passed)
  expect_false(context_4d_preflight$checks$passed[
    context_4d_preflight$checks$check == "endpoint_context_is_3d:anatomy"
  ])
})

test_that("direct FibrTraQ validates every used reference image against fNode", {
  skip_if_not_installed("RNifti")

  dims <- c(3L, 3L, 3L)
  focal_array <- array(0, dim = dims)
  focal_array[1, 1, 1] <- 1
  target_array <- array(0, dim = dims)
  target_array[3, 1, 1] <- 1
  focal <- ftq_test_nifti_image(focal_array, spatial_units = "mm")
  target <- ftq_test_nifti_image(target_array, spatial_units = "mm")
  mismatched_dims <- ftq_test_nifti_image(
    array(0, dim = c(4, 3, 3)),
    spatial_units = "mm"
  )
  four_dimensional <- ftq_test_nifti_image(
    array(0, dim = c(3, 3, 3, 2)),
    spatial_units = "mm"
  )
  mismatched_affine <- ftq_test_nifti_image(
    array(0, dim = dims),
    shift = c(1, 0, 0),
    spatial_units = "mm"
  )
  tracks <- list(rbind(c(1, 1, 1), c(3, 1, 1)))

  scenarios <- list(
    full_world = list(
      assignment_mode = "full_streamline",
      track_space = "world"
    ),
    endpoint_trackvis = list(
      assignment_mode = "endpoint_voxels",
      track_space = "trackvis",
      track_affine = diag(4L)
    ),
    full_auto = list(
      assignment_mode = "full_streamline",
      track_space = "auto"
    ),
    local_voxel = list(
      assignment_mode = "local_search",
      track_space = "voxel"
    )
  )
  run_scenario <- function(scenario, reference) {
    do.call(FibrTraQ, c(
      list(
        fNode = focal,
        localROI = NULL,
        cNode_list = list(T = target),
        cNode_labels = "T",
        all_tracks = tracks,
        ref_image = reference,
        discretization = "round",
        focal_assignment_mode = "all_focal_voxels",
        streamline_voxelization_mode = "sample_points",
        axis_flip = c(FALSE, FALSE, FALSE),
        backend = "r",
        endpoint_analysis = "none"
      ),
      scenario
    ))
  }

  for (scenario_name in names(scenarios)) {
    scenario <- scenarios[[scenario_name]]
    expect_error(
      run_scenario(scenario, mismatched_dims),
      "ref_image dimensions must match the fNode grid",
      fixed = TRUE,
      info = scenario_name
    )
    expect_error(
      run_scenario(scenario, four_dimensional),
      "ref_image must provide a three-dimensional voxel grid",
      fixed = TRUE,
      info = scenario_name
    )
    expect_error(
      run_scenario(scenario, mismatched_affine),
      "ref_image must use the same unit-normalized affine/grid as fNode",
      fixed = TRUE,
      info = scenario_name
    )
  }

  focal_mm <- ftq_test_nifti_image(
    focal_array,
    scale = 1000,
    spatial_units = "mm"
  )
  equivalent_m <- ftq_test_nifti_image(
    array(0, dim = dims),
    scale = 1,
    spatial_units = "m"
  )
  expect_no_error(FibrTraQ(
    fNode = focal_mm,
    localROI = NULL,
    cNode_list = list(T = target_array),
    cNode_labels = "T",
    all_tracks = tracks,
    discretization = "round",
    focal_assignment_mode = "all_focal_voxels",
    streamline_voxelization_mode = "sample_points",
    track_space = "auto",
    ref_image = equivalent_m,
    axis_flip = c(FALSE, FALSE, FALSE),
    backend = "r",
    endpoint_analysis = "none"
  ))
})

test_that("preflight checks explicit mapping references with normalized units", {
  skip_if_not_installed("RNifti")

  tmp <- tempfile("ftq-preflight-reference-grid-")
  dir.create(tmp)
  fnode <- file.path(tmp, "fnode-mm.nii.gz")
  target <- file.path(tmp, "target-mm.nii.gz")
  equivalent_ref <- file.path(tmp, "reference-m.nii.gz")
  shifted_ref <- file.path(tmp, "reference-shifted.nii.gz")
  mismatched_ref <- file.path(tmp, "reference-dims.nii.gz")
  undefined_fnode <- file.path(tmp, "fnode-undefined-affine.nii.gz")
  undefined_target <- file.path(tmp, "target-undefined-affine.nii.gz")
  tractogram <- file.path(tmp, "tracts.trk")

  focal_array <- array(0, dim = c(3, 3, 3))
  focal_array[1, 1, 1] <- 1
  target_array <- array(0, dim = c(3, 3, 3))
  target_array[3, 1, 1] <- 1
  RNifti::writeNifti(
    ftq_test_nifti_image(focal_array, scale = 1000, spatial_units = "mm"),
    fnode
  )
  RNifti::writeNifti(
    ftq_test_nifti_image(target_array, scale = 1000, spatial_units = "mm"),
    target
  )
  RNifti::writeNifti(
    ftq_test_nifti_image(
      array(0, dim = c(3, 3, 3)),
      scale = 1,
      spatial_units = "m"
    ),
    equivalent_ref
  )
  RNifti::writeNifti(
    ftq_test_nifti_image(
      array(0, dim = c(3, 3, 3)),
      shift = c(1, 0, 0),
      scale = 1000,
      spatial_units = "mm"
    ),
    shifted_ref
  )
  RNifti::writeNifti(
    ftq_test_nifti_image(
      array(0, dim = c(4, 3, 3)),
      scale = 1000,
      spatial_units = "mm"
    ),
    mismatched_ref
  )
  RNifti::writeNifti(focal_array, undefined_fnode)
  RNifti::writeNifti(target_array, undefined_target)
  file.create(tractogram)

  manifest <- list(
    dataset_root = tmp,
    subject = "sub-reference",
    description = "mapping-grid",
    fnode = basename(fnode),
    tractogram = basename(tractogram),
    ref_image = basename(equivalent_ref),
    targets = list(T = basename(target)),
    coordinate = list(track_space = "auto"),
    outputs = list(output_dir = file.path(tmp, "out"), overwrite = FALSE)
  )
  equivalent <- ftq_preflight(manifest)
  expect_true(equivalent$passed)
  expect_true(equivalent$checks$passed[
    equivalent$checks$check == "mapping_ref_image_affine_matches_fnode"
  ])

  shifted_manifest <- manifest
  shifted_manifest$ref_image <- basename(shifted_ref)
  shifted <- ftq_preflight(shifted_manifest)
  expect_false(shifted$passed)
  expect_false(shifted$checks$passed[
    shifted$checks$check == "mapping_ref_image_affine_matches_fnode"
  ])

  mismatched_manifest <- manifest
  mismatched_manifest$ref_image <- basename(mismatched_ref)
  mismatched <- ftq_preflight(mismatched_manifest)
  expect_false(mismatched$passed)
  expect_false(mismatched$checks$passed[
    mismatched$checks$check == "mapping_ref_image_dims_match_fnode"
  ])

  for (mapping in c("world", "trackvis")) {
    mapping_manifest <- shifted_manifest
    mapping_manifest$coordinate <- list(track_space = mapping)
    mapping_preflight <- ftq_preflight(mapping_manifest)
    expect_false(mapping_preflight$passed, info = mapping)
    expect_false(mapping_preflight$checks$passed[
      mapping_preflight$checks$check ==
        "mapping_ref_image_affine_matches_fnode"
    ], info = mapping)
  }

  undefined_manifest <- manifest
  undefined_manifest$fnode <- basename(undefined_fnode)
  undefined_manifest$targets <- list(T = basename(undefined_target))
  undefined_manifest$ref_image <- basename(shifted_ref)
  undefined_manifest$outputs$output_dir <- file.path(tmp, "undefined-out")
  undefined <- ftq_preflight(undefined_manifest)
  expect_true(undefined$passed)
  undefined_ref_row <- undefined$checks[
    undefined$checks$check == "mapping_ref_image_affine_matches_fnode",
    ,
    drop = FALSE
  ]
  expect_equal(nrow(undefined_ref_row), 1L)
  expect_equal(undefined_ref_row$severity, "warning")
  expect_false(undefined_ref_row$passed)
  expect_match(undefined_ref_row$detail, "affine/grid comparison is unavailable")
})
