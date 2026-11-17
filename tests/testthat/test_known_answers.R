test_that("the single-subject command reproduces the bundled reference tables and claims", {
  fixture <- ftq_test_example_files()
  p <- fixture$parameters
  root <- tempfile("known-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  analysis <- fibrtraq_single(
    tractogram = fixture$tractogram, focal_mask = fixture$focal_mask,
    target_masks = fixture$target_masks, local_mask = fixture$local_mask,
    reference_image = fixture$reference_image, endpoint_contexts = fixture$endpoint_contexts,
    output_dir = root, subject = p$subject, focal_label = p$focal_label,
    track_space = p$track_space, discretization = p$discretization,
    axis_flip = as.logical(p$axis_flip), description = p$description,
    focal_assignment_mode = p$focal_assignment_mode,
    streamline_voxelization_mode = p$streamline_voxelization_mode,
    assignment_mode = p$target_intersection_mode, endpoint_analysis = p$endpoint_analysis,
    track_reader = p$track_reader, track_chunk_size = p$track_chunk_size,
    backend = "r", verbose = FALSE, progress = FALSE
  )
  expect_true(analysis$index$validation_passed)

  read_science <- function(path) {
    tbl <- utils::read.csv(path, check.names = FALSE)
    if ("source_path" %in% names(tbl)) tbl$source_path <- basename(tbl$source_path)
    tbl
  }
  expected <- lapply(fixture$expected, read_science)
  actual <- lapply(names(expected), function(name) {
    read_science(analysis$index[[paste0(name, "_csv")]])
  })
  names(actual) <- names(expected)
  expect_length(actual, 7L)
  for (name in names(actual)) expect_identical(actual[[name]], expected[[name]], label = name)

  claims <- utils::read.csv(file.path(fixture$root, "known-answer-claims.csv"))
  expect_equal(nrow(claims), 18L)
  for (i in seq_len(nrow(claims))) {
    claim <- claims[i, ]
    tbl <- actual[[claim$source_table]]
    for (j in 1:2) {
      field <- claim[[paste0("filter_column_", j)]]
      if (!is.na(field) && nzchar(field)) {
        tbl <- tbl[tbl[[field]] == claim[[paste0("filter_value_", j)]], , drop = FALSE]
      }
    }
    expect_equal(tbl[[claim$value_column]], claim$expected_value, info = claim$claim)
  }
})
