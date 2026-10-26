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
