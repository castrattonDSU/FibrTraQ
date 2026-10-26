test_that("compact context preserves research density and metadata without sidecar files", {
  skip_if_not_installed("RNifti")
  root <- tempfile("tractogram-context-layout-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  inputs <- system.file("extdata", "fibrtraq-known-answer", "inputs", package = "FibrTraQ")
  args <- list(
    tractogram = file.path(inputs, "sub-example_tractogram.trk"),
    template_image = file.path(inputs, "sub-example_focal.nii.gz"),
    prefix = "fixture", track_space = "voxel", track_chunk_size = 2L,
    streamline_voxelization_mode = "segment_interpolation"
  )
  research <- do.call(ftq_tractogram_context, c(args, list(output_dir = file.path(root, "research"))))
  compact <- do.call(FibrTraQ:::ftq_build_tractogram_context, c(args, list(
    output_dir = file.path(root, "public"), density_subdir = "maps", write_sidecars = FALSE
  )))
  expect_true(all(file.exists(unlist(research[c("density_image", "summary_csv", "metadata_json")]))))
  expect_identical(basename(dirname(research$density_image)), "context")
  expect_identical(basename(dirname(compact$density_image)), "maps")
  expect_null(compact$summary_csv)
  expect_null(compact$metadata_json)
  expect_length(list.files(file.path(root, "public"), recursive = TRUE), 1L)
  compact_density <- RNifti::readNifti(compact$density_image)
  research_density <- RNifti::readNifti(research$density_image)
  expect_true(FibrTraQ:::ftq_same_grid(compact_density, research_density))
  expect_equal(as.numeric(compact_density), as.numeric(research_density))
  expect_equal(compact$summary, research$summary)
  metadata_fields <- setdiff(names(research$metadata), "generated_at")
  expect_equal(compact$metadata[metadata_fields], research$metadata[metadata_fields])
  expect_error(do.call(FibrTraQ:::ftq_build_tractogram_context, c(args, list(
    output_dir = file.path(root, "public"), density_subdir = "maps", write_sidecars = FALSE
  ))), "already exist")
})
