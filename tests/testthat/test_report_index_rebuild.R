test_that("complete report links and figures survive reserved output-folder names", {
  root <- file.path(tempdir(), " #%20\u00e9")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  example <- ftq_run_example(output_dir = root)
  expect_true(example$validation$passed)
  reports <- list.files(root, pattern = "[.]html$", recursive = TRUE, full.names = TRUE)
  expect_gte(length(reports), 4L)
  for (report in reports) {
    html <- paste(readLines(report, warn = FALSE), collapse = "\n")
    attributes <- regmatches(html, gregexpr("(?:href|src)=['\"][^'\"]+['\"]", html,
                                            perl = TRUE))[[1L]]
    for (attribute in attributes) {
      href <- sub("^[^=]+=['\"](.*)['\"]$", "\\1", attribute)
      if (grepl("^(#|https?:|data:)", href)) next
      expect_true(file.exists(file.path(dirname(report), utils::URLdecode(href))),
                  info = paste(basename(report), href))
    }
  }
})

test_that("a Windows summary rebuilds after relocation and retains leading-zero IDs", {
  root <- tempfile("moved reports # 100%20 ")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  paths <- file.path(root, "001", "002", c(
    "reports/standard.html", "exploratory-report/exploratory.html",
    "reports/figures/viewer.html"
  ))
  for (path in paths) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines("<!doctype html><title>Existing report</title>", path)
  }
  index <- data.frame(
    subject = "001", focal_region = "002", run_status = "completed",
    validation_passed = TRUE, biological_qc_status = "review",
    primary_report = "C:\\old root\\001\\002\\reports\\standard.html",
    exploratory_report = "C:\\old root\\001\\002\\exploratory-report\\exploratory.html",
    interactive_3d_qc = "C:\\old root\\001\\002\\reports\\figures\\viewer.html"
  )
  csv <- file.path(root, "fibrtraq-run-index.csv")
  original <- file.path(root, "fibrtraq-report-index.html")
  write.csv(index, csv, row.names = FALSE)
  writeLines("Original summary retained", original)
  hashes <- tools::md5sum(c(csv, original, paths))

  repaired <- FibrTraQ:::ftq_rebuild_report_index(root)
  expect_identical(basename(repaired), "fibrtraq-report-index-portable.html")
  html <- paste(readLines(repaired), collapse = "\n")
  hrefs <- regmatches(html, gregexpr("href='[^']+'", html))[[1L]]
  hrefs <- sub("^href='(.*)'$", "\\1", hrefs)
  expect_length(hrefs, 5L)
  expect_setequal(
    normalizePath(file.path(root, utils::URLdecode(hrefs)), winslash = "/"),
    normalizePath(c(paths, csv, file.path(root, "OUTPUTS_README.html")), winslash = "/")
  )
  expect_identical(tools::md5sum(c(csv, original, paths)), hashes)
  expect_error(FibrTraQ:::ftq_rebuild_report_index(root), "overwrite|exist", ignore.case = TRUE)

  unlink(repaired)
  duplicate <- file.path(dirname(paths[1L]), "archive", basename(paths[1L]))
  dir.create(dirname(duplicate))
  file.copy(paths[1L], duplicate)
  expect_error(FibrTraQ:::ftq_rebuild_report_index(root), "found 2")
  expect_false(file.exists(repaired))
  unlink(duplicate)
  unlink(paths[1L])
  expect_error(FibrTraQ:::ftq_rebuild_report_index(root), "found 0")
  expect_false(file.exists(repaired))
})
