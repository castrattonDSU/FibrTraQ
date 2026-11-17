ftq_drop_result_attrs <- function(x) {
  keep <- attributes(x)[c("names", "row.names", "class")]
  attributes(x) <- keep
  x
}

ftq_expect_core_equal <- function(actual, expected) {
  testthat::expect_equal(ftq_drop_result_attrs(actual), ftq_drop_result_attrs(expected))
}

ftq_expect_algorithm_metadata_equal <- function(actual, expected) {
  for (nm in c(
    "fibrtraq_focal_assignment_mode",
    "fibrtraq_streamline_voxelization_mode",
    "fibrtraq_target_intersection_mode",
    "fibrtraq_track_space",
    "fibrtraq_axis_flip",
    "fibrtraq_discretization",
    "fibrtraq_reverse_streamlines",
    "fibrtraq_mapping_selection_rule"
  )) {
    testthat::expect_equal(attr(actual, nm), attr(expected, nm), info = nm)
  }
}

ftq_test_path <- function(path) {
  gsub("\\\\", "/", path)
}

ftq_test_example_files <- function() {
  root <- system.file("extdata", "fibrtraq-known-answer", package = "FibrTraQ",
                      mustWork = TRUE)
  input <- function(name) file.path(root, "inputs", paste0("sub-example_", name))
  expected <- list.files(file.path(root, "expected"), pattern = "[.]csv$", full.names = TRUE)
  list(
    root = root, parameters = jsonlite::fromJSON(file.path(root, "analysis-parameters.json")),
    tractogram = input("tractogram.trk"), focal_mask = input("focal.nii.gz"),
    target_masks = c(cortex = input("target-cortex.nii.gz"),
                     brainstem = input("target-brainstem.nii.gz")),
    local_mask = input("local.nii.gz"), reference_image = input("reference.nii.gz"),
    endpoint_contexts = c(analysis_domain = input("context.nii.gz")),
    expected = stats::setNames(expected, tools::file_path_sans_ext(basename(expected)))
  )
}
