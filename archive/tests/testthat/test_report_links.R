ftq_report_link_test_index <- function(reports) {
  data.frame(
    subject = "subject1", focal_region = "region1", run_status = "completed",
    validation_passed = TRUE, biological_qc_status = "not_reviewed",
    primary_report = reports[[1L]], exploratory_report = reports[[2L]],
    interactive_3d_qc = reports[[3L]], stringsAsFactors = FALSE
  )
}

ftq_report_link_test_hrefs <- function(html) {
  links <- regmatches(html, gregexpr("href=['\"][^'\"]+['\"]", html))[[1L]]
  sub("['\"]$", "", sub("^href=['\"]", "", links))
}

test_that("regional index handles native Windows strings and the original schema", {
  # Capture only the write boundary; all path processing uses the real writer.
  html <- NULL
  writer <- ftq_write_regional_study_index
  environment(writer) <- list2env(list(
    writeLines = function(text, con, ...) {
      html <<- paste(text, collapse = "\n")
    }
  ), parent = environment(writer))
  reports <- c("C:\\study #\\output\\subject 1\\standard %20.html",
               "C:\\study #\\output\\subject 1\\exploratory.html",
               "C:\\study #\\output\\subject 1\\QC.html")
  for (native_windows in c(TRUE, FALSE)) {
    index_path <- "C:\\study #\\output\\indexes\\index.html"
    output_dir <- "C:\\study #\\output"
    report_paths <- reports
    if (!native_windows) {
      index_path <- gsub("\\", "/", index_path, fixed = TRUE)
      output_dir <- gsub("\\", "/", output_dir, fixed = TRUE)
      report_paths <- gsub("\\", "/", report_paths, fixed = TRUE)
    }
    writer(ftq_report_link_test_index(report_paths), path = index_path,
           output_dir = output_dir, description = "link-regression")
    expect_identical(ftq_report_link_test_hrefs(html),
                     c("../subject%201/standard%20%2520.html",
                       "../subject%201/exploratory.html", "../subject%201/QC.html"))
    expect_match(html, "not recorded", fixed = TRUE)
  }
})

test_that("raw regional hrefs survive relocation of the complete report tree", {
  parent <- tempfile("report-relocation-")
  original <- file.path(parent, "original tree # %20")
  relocated <- file.path(parent, "relocated tree # %20")
  dir.create(file.path(original, "indexes"), recursive = TRUE)
  dir.create(file.path(original, "subject 1 #"))
  on.exit(unlink(parent, recursive = TRUE), add = TRUE)
  names <- c("standard %20.html", "exploratory #.html", "QC caf\u00e9.html")
  reports <- file.path(original, "subject 1 #", names)
  for (i in seq_along(reports)) writeLines(paste("report", i), reports[[i]])
  historical_guide <- file.path(original, "OUTPUTS_README.html")
  writeLines("Original field definitions.", historical_guide)
  index <- file.path(original, "indexes", "index.html")
  ftq_write_regional_study_index(
    ftq_report_link_test_index(reports), path = index, output_dir = original,
    description = "relocation-regression"
  )
  html <- paste(readLines(index, warn = FALSE), collapse = "\n")
  hrefs <- ftq_report_link_test_hrefs(html)
  expect_identical(hrefs, c("../subject%201%20%23/standard%20%2520.html",
                           "../subject%201%20%23/exploratory%20%23.html",
                           "../subject%201%20%23/QC%20caf%C3%A9.html",
                           "../OUTPUTS_README.html"))
  guide_content <- readLines(file.path(original, "OUTPUTS_README.html"), warn = FALSE)
  expect_identical(guide_content, "Original field definitions.")
  expect_true(file.rename(original, relocated))
  expect_false(dir.exists(original))
  moved_index <- file.path(relocated, "indexes", "index.html")
  expect_identical(paste(readLines(moved_index, warn = FALSE), collapse = "\n"), html)
  for (i in seq_along(hrefs)) {
    target <- file.path(dirname(moved_index), utils::URLdecode(hrefs[[i]]))
    expect_true(file.exists(target))
    expect_identical(readLines(target, warn = FALSE),
                     if (i <= length(reports)) paste("report", i) else guide_content)
    expect_true(startsWith(normalizePath(target, winslash = "/"),
                           paste0(normalizePath(relocated, winslash = "/"), "/")))
  }
})
