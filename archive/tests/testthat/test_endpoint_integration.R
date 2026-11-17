test_that("archived statistics accept endpoint results", {
  result <- ftq_run_endpoint_fixture("cpp")
  endpoints <- ftq_endpoint_tables(result)
  old_subject <- ftq_stat_subject_summary(result)
  expect_identical(names(old_subject), c("roi", "targets", "voxels"))
  endpoint_subject <- ftq_stat_subject_summary(
    result,
    endpoint_summary = endpoints$summary,
    endpoint_penetration = endpoints$penetration
  )
  expect_true(all(c("endpoints", "penetration") %in% names(endpoint_subject)))

  old_targets <- ftq_target_summary_tables(result, attr(result, "fibrtraq_roi_summary"))
  expect_identical(names(old_targets), c("unique_streamlines", "voxel_associations"))
  endpoint_targets <- ftq_target_summary_tables(
    result,
    attr(result, "fibrtraq_roi_summary"),
    endpoint_summary = endpoints$summary
  )
  expect_true("endpoint_streamlines" %in% names(endpoint_targets))
})
