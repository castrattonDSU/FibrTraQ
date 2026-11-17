test_that("the archived study preserves the 0.3.5 regional-union default", {
  expect_identical(formals(research$ftq_run_regional_study)$local_masks, "union")
  expect_null(formals(FibrTraQ::fibrtraq_study)$local_masks)
})

test_that("index shares method context and distinguishes example file checks", {
  root <- tempfile("public-index-")
  dir.create(root)
  report <- file.path(root, "index.html")
  row <- FibrTraQ:::ftq_regional_study_run_row(
    "001", "A", list(validation = list(passed = TRUE)), NULL, NULL, "start", NULL,
    manifest = list()
  )
  row$target_intersection_mode <- "local_search"
  row$local_search_radius_mm <- 2
  row$focal_assignment_mode <- "all_focal_voxels"
  index <- rbind(row, transform(row, focal_region = "B"))
  checks <- data.frame(artifact = "example", passed = TRUE, details = "matched")
  validation <- list(package_validation_passed = TRUE, tables = checks,
                     claims = transform(checks, claim = "claim"), reports = checks)
  ftq_write_regional_study_index(index, report, root, "test", validation)
  html <- paste(readLines(report), collapse = "\n")
  expect_length(regmatches(html, gregexpr("radius 2 mm", html, fixed = TRUE))[[1L]], 1L)
  expect_match(html, "Known-answer comparisons: 1/1 tables and 1/1 claims", fixed = TRUE)
  expect_match(html, "Generated HTML files present and nonempty: 1/1", fixed = TRUE)
  expect_match(html, "browser rendering and interactive behavior require separate review", fixed = TRUE)
})
