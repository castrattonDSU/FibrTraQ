test_that("bundled known-answer inputs are complete and checksummed", {
  example <- ftq_example_data()

  expect_equal(example$id, "fibrtraq-known-answer-v1")
  expect_true(file.exists(example$tractogram))
  expect_true(file.exists(example$focal_mask))
  expect_true(all(file.exists(example$target_masks)))
  expect_true(file.exists(example$local_mask))
  expect_true(file.exists(example$reference_image))
  expect_true(all(file.exists(example$endpoint_contexts)))
  expect_true(all(file.exists(example$expected)))
  expect_equal(length(example$expected), 7L)
  expected_claims <- data.frame(
    claim = c(
      "processed_streamlines",
      "unique_focal_intersecting_streamlines",
      "total_focal_voxel_streamline_associations",
      "unique_any_target_positive_streamlines",
      "unique_local_streamlines",
      "unique_leave_no_target_streamlines",
      "cortex_target_intersections",
      "brainstem_target_intersections",
      "cortex_focal_voxel_target_associations",
      "brainstem_focal_voxel_target_associations",
      "cortex_exact_paired_endpoints",
      "cortex_dilated1_paired_endpoints",
      "brainstem_exact_paired_endpoints",
      "reference_grid_endpoint_count",
      "reference_grid_endpoints_inside",
      "fate_local_associations",
      "fate_leave_no_target_associations",
      "fate_leave_with_target_associations"
    ),
    source_table = c(
      rep("roi_summary", 6L),
      rep("endpoint_summary", 2L),
      rep("rsc", 2L),
      rep("endpoint_summary", 3L),
      rep("endpoint_context", 2L),
      rep("fate", 3L)
    ),
    value_column = c(
      "processed_streamlines",
      "unique_focal_intersecting_streamlines",
      "total_focal_voxel_streamline_associations",
      "unique_any_target_positive_streamlines",
      "unique_local_streamlines",
      "unique_leave_no_target_streamlines",
      "target_intersection_streamlines",
      "target_intersection_streamlines",
      "target_associations",
      "target_associations",
      "paired_focal_target_endpoints",
      "paired_focal_target_endpoints",
      "paired_focal_target_endpoints",
      "endpoint_count",
      "endpoints_inside",
      "tracks",
      "tracks",
      "tracks"
    ),
    filter_column_1 = c(
      rep(NA_character_, 6L),
      "target", "target", "node", "node", "target", "target", "target",
      "context", "context", "fate", "fate", "fate"
    ),
    filter_value_1 = c(
      rep(NA_character_, 6L),
      "cortex", "brainstem", "cortex", "brainstem", "cortex", "cortex",
      "brainstem", "reference_grid", "reference_grid", "local",
      "leave_no_target", "leave_with_target"
    ),
    filter_column_2 = c(
      rep(NA_character_, 6L),
      "tolerance", "tolerance", NA_character_, NA_character_,
      "tolerance", "tolerance", "tolerance", rep(NA_character_, 5L)
    ),
    filter_value_2 = c(
      rep(NA_character_, 6L),
      "exact", "exact", NA_character_, NA_character_,
      "exact", "dilated1", "exact", rep(NA_character_, 5L)
    ),
    expected_value = c(
      14L, 12L, 15L, 10L, 1L, 1L, 6L, 4L, 6L, 7L, 2L, 3L, 2L,
      27L, 27L, 1L, 1L, 13L
    ),
    stringsAsFactors = FALSE
  )
  expect_identical(example$claims, expected_claims)
  expect_equal(example$parameters$biological_interpretation, FALSE)
  expect_equal(
    example$parameters$target_intersection_mode,
    "full_streamline"
  )
  expect_equal(example$parameters$example_execution$sample_per_target, 6L)
  expect_equal(nrow(example$inventory), 17L)

  text_paths <- c(
    file.path(example$root, "checksums.csv"),
    example$inventory$path[
      grepl("\\.(csv|json|md)$", example$inventory$relative_path)
    ]
  )
  contains_cr <- vapply(text_paths, function(path) {
    bytes <- readBin(path, what = "raw", n = file.info(path)$size)
    any(bytes == as.raw(13L))
  }, logical(1))
  expect_false(any(contains_cr), info = paste(text_paths[contains_cr],
                                              collapse = ", "))
  expect_error(ftq_example_data(verify = NA), "verify must be TRUE or FALSE")
})

test_that("bundled example reproduces every scientific table and report", {
  output_dir <- tempfile("f-")
  result <- ftq_run_example(
    output_dir = output_dir,
    backend = "auto",
    overwrite = TRUE,
    verbose = FALSE
  )

  expect_s3_class(result, "ftq_example_result")
  expect_true(result$validation$passed)
  expect_true(result$validation$package_validation_passed)
  expect_equal(nrow(result$validation$tables), 7L)
  expect_true(all(result$validation$tables$passed))
  expect_equal(nrow(result$validation$claims), 18L)
  expect_true(all(result$validation$claims$stored_passed))
  expect_true(all(result$validation$claims$actual_passed))
  expect_true(all(result$validation$claims$passed))
  expect_true(all(result$validation$reports$passed))
  expect_true(file.exists(result$report_index))
  expect_true(file.exists(result$primary_report))
  expect_true(file.exists(result$exploratory_report))
  expect_true(file.exists(result$interactive_3d_qc))
  params <- result$data$parameters
  planned <- ftq_regional_study_output_paths(output_dir, params$subject,
    c(params$focal_label, names(result$data$target_masks)), params$focal_label,
    params$description)
  actual_files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  # A new report artifact must also enter preflight, so future filename changes
  # cannot silently escape the Windows budget check.
  expect_length(setdiff(normalizePath(actual_files, winslash = "/"), planned), 0L)
  expect_output(print(result), "Exact tables: 7/7", fixed = TRUE)
  expect_output(print(result), "Known answers: 18/18", fixed = TRUE)

  viewer <- paste(
    readLines(result$interactive_3d_qc, warn = FALSE),
    collapse = "\n"
  )
  expect_match(viewer, "intersection-positive QC sample", fixed = TRUE)
  expect_match(viewer, "Exact", fixed = TRUE)
  expect_match(viewer, "1-voxel sensitivity", fixed = TRUE)
})

test_that("known-answer comparison detects scientific drift", {
  example <- ftq_example_data()
  changed <- tempfile(fileext = ".csv")
  table <- utils::read.csv(
    example$expected[["result"]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  table$tracks[1] <- table$tracks[1] + 1
  utils::write.csv(table, changed, row.names = FALSE)

  actual <- example$expected
  actual[["result"]] <- changed
  comparison <- FibrTraQ:::ftq_example_compare_tables(
    actual,
    example$expected
  )

  expect_false(comparison$passed[comparison$artifact == "result"])
  expect_true(all(comparison$passed[comparison$artifact != "result"]))
})

test_that("known-answer claim validation detects declaration drift", {
  example <- ftq_example_data()
  changed_claims <- example$claims
  changed_claims$expected_value[1L] <-
    changed_claims$expected_value[1L] + 1

  comparison <- FibrTraQ:::ftq_example_check_claims(
    actual = example$expected,
    expected = example$expected,
    claims = changed_claims
  )

  expect_false(comparison$stored_passed[1L])
  expect_false(comparison$actual_passed[1L])
  expect_false(comparison$passed[1L])
  expect_true(all(comparison$passed[-1L]))
})
