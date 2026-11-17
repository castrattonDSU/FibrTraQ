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
